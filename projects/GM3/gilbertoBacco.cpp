#include <perens.hpp>

#include "highPrec.hpp"

// https://arxiv.org/pdf/1903.06476.pdf

using PrecVect=
  Matrix<PrecFloat,Dynamic,1>;

using PrecMatr=
  Matrix<PrecFloat,Dynamic,Dynamic>;

// template <typename F>
// PrecFloat integrateUpToInfinity(const F& f,
// 				const PrecFloat& min)
// {
//   PrecFloat s=0;
  
//   PrecFloat dX=1e-4;
//   for(PrecFloat x=1;x>0;x-=dX)
//     s+=f(min+(PrecFloat(1.0)-x)/x)/x/x;
  
//   s*=dX;
  
//   cout<<s.get()<<endl;
  
//   return
//     s;
// }

void perens_t::Gilberto() const
{
  PrecFloat::setDefaultPrecision(512);
  
  constexpr bool useBw=1;
  
  const double eu=2.0/3,ed=-1.0/3;
   djvec_t corr=
     getAveForRego(0,nSources,1,REGO_TM)*sqr(Z[regoZId[REGO_TM]])*(sqr(eu)+sqr(ed));
   
  // djack_t Z2(gauss_filler_t{0.1,0.01,234324});
  // djack_t M=djack_t(gauss_filler_t{0.4,0.01,235253});
  
  // for(size_t t=0;t<=TH;t++)
  //   corr[t]=two_pts_corr_fun(Z2,M,TH,t,+1);
  
  auto bT=
    [this](const PrecFloat& t,const PrecFloat& E)
    {
      return
	exp(-t*E)+useBw*
	exp(-((PrecFloat)T-t)*E);
    };
  
  const int nT=48,tMin=1;
  
  PrecVect R(nT);
  for(int iT=0;iT<nT;iT++)
    R[iT]=
      1.0/PrecFloat(iT+tMin)+useBw*
      1.0/PrecFloat(T-iT-tMin);
  
  grace_file_t RFile("/tmp/R.xmg");
  RFile.new_data_set();
  for(int iT=0;iT<nT;iT++)
    RFile.write_xy(iT,R[iT].get());
  
  /////////////////////////////////////////////////////////////////
  
  PrecMatr A(nT,nT);
  const PrecFloat E0=0.1;
  
  for(int iR=0;iR<nT;iR++)
    for(int iT=0;iT<nT;iT++)
      {
	const PrecFloat a=(PrecFloat)iR+iT+2*tMin;
	const PrecFloat b=(PrecFloat)T-iR+iT;
	const PrecFloat c=(PrecFloat)T+iR-iT;
	const PrecFloat d=(PrecFloat)2*T-iR-iT-2*tMin;
	A(iR,iT)=
	  exp(-a*E0)/a+
	  useBw*(
	     exp(-b*E0)/b+
	     exp(-c*E0)/c+
	     exp(-d*E0)/d);
      }
  
  /////////////////////////////////////////////////////////////////
  
  const PrecFloat lambda=0;//0.5;
  
  const auto Z=
    [](const PrecFloat& Estar,
       const PrecFloat& sigma)
    {
      return
	(PrecFloat(1)+erf(Estar/(M_SQRT2*sigma)))/2;
    };
  
  const PrecFloat sigma=0.1;
  const PrecFloat alpha=0.0;
  
  const auto N=
    [&sigma,&alpha,&lambda,&Z](const PrecFloat& Estar,
			       const PrecFloat& k)
    {
      return
	(1-lambda)/(2*Z(Estar,sigma))*exp((alpha-k)*((alpha-k)*sigma*sigma+2*Estar)/2);
    };
  
  const auto F=
    [&sigma,&alpha,&E0](const PrecFloat& Estar,
			const PrecFloat& k)
    {
      return
	1+erf(((alpha-k)*sigma*sigma+Estar-E0)/(M_SQRT2*sigma));
    };
  
  grace_file_t out("/tmp/E.xmg");
  for(double Estar=1.5;Estar<1.500001;Estar+=0.01)
    {
      const auto DeltaFun=
	[&sigma,&Estar,&Z](const PrecFloat& E)
	{
	  return
	    exp(-sqr(E-Estar)/(2*sigma*sigma))/(M_SQRT2*sqrt(M_PI)*sigma*Z(Estar,sigma));
	};
      
      PrecVect f(nT);
      for(int iT=0;iT<nT;iT++)
	f(iT)=
	  N(Estar,iT+tMin)*F(Estar,iT+tMin)+useBw*
	  N(Estar,T-iT-tMin)*F(Estar,T-iT-tMin);
      
      PrecMatr W(nT,nT);
      
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  W(iR,iT)=
	    (1-lambda)*A(iR,iT)+
	     lambda*(iR==iT)*sqr(corr[iR+tMin].err()/corr[0].ave());
      
      grace_file_t fFile("/tmp/f"+to_string(Estar)+".xmg");
      for(int i=0;i<nT;i++)
	fFile.write_xy(i,f(i).get());
      
      const auto Winv=
	W.inverse();
      
      const PrecFloat num=
	1-R.transpose()*Winv*f;
      
      const PrecFloat den=
	R.transpose()*Winv*R;
      
      const PrecVect g=
	Winv*f+
	Winv*R*num/den;
      
      // const MatrixXd Sdiag=
      //   svd.singularValues().asDiagonal();
      
      // const auto AinvBis=
      //   svd.matrixV()*Sdiag.inverse()*svd.matrixU().transpose();
      
      // cout<<(Ainv-AinvBis).eval()<<endl;
      
      // const auto num=
      // 	(Ainv*R).eval();
      
      // const auto den=
      // 	(R.transpose()*Ainv*R)(0);
      
      grace_file_t gFile("/tmp/g"+to_string(Estar)+".xmg");
      // const auto g=
      // 	(num/den).eval();
      for(int iT=0;iT<nT;iT++)
	gFile.write_xy(iT,g[iT].get());
      
      grace_file_t RecDelta("/tmp/RecDelta"+to_string(Estar)+".xmg");
      grace_file_t TargDelta("/tmp/TargDelta"+to_string(Estar)+".xmg");
      grace_file_t ErrDelta("/tmp/ErrDelta"+to_string(Estar)+".xmg");
      PrecFloat x=0,x2=0,in=0;
      for(PrecFloat E=E0;E<2.0;E+=0.01)
	{
	  PrecFloat s=0,sTarg=DeltaFun(E);
	  for(int iT=0;iT<nT;iT++)
	    s+=g(iT)*bT(iT+tMin,E);
	  RecDelta.write_xy(E.get(),s.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=s-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	  
	  x+=E*s*0.01;
	  x2+=E*E*s*0.01;
	  in+=s*0.01;
	}
      x/=in;
      x2/=in;
      x2-=x*x;
      // Cent.write_xy(Estar,x.get());
      // Res.write_xy(Estar,sqrt(x2).get());
      
      djack_t s{};
      for(int iT=0;iT<nT;iT++)
	for(size_t ijack=0;ijack<=njacks;ijack++)
	  s[ijack]+=(g[iT]*corr[iT+tMin][ijack]).get();
      
      out.write_ave_err(Estar,s.ave_err());
    }
  
  cout<<"Gilbert done"<<endl;
}
