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
  PrecFloat::setDefaultPrecision(128);
  
  const double eu=2.0/3,ed=-1.0/3;
   djvec_t corr=
     getAveForRego(0,nSources,1,REGO_TM)*sqr(Z[regoZId[REGO_TM]])*(sqr(eu)+sqr(ed));
   
  // djack_t Z2(gauss_filler_t{0.1,0.01,234324});
  // djack_t M=djack_t(gauss_filler_t{0.4,0.01,235253});
  
  // for(size_t t=0;t<=TH;t++)
  //   corr[t]=two_pts_corr_fun(Z2,M,TH,t,+1);
  
  auto bT=
    [this](const int& t,const PrecFloat& E)
    {
      return
	(exp(-(PrecFloat)t*E)+exp(-(PrecFloat)(T-t)*E));
    };
  
  const int nT=10,tMin=1;
  PrecVect R(nT);
  for(int it=0;it<nT;it++)
    {
      const int t=
	it+tMin;
      R[it]=
	PrecFloat(1.0)/t+PrecFloat(1.0)/(T-t);
    }
  
  {
    grace_file_t RFile("/tmp/R.xmg");
    RFile.new_data_set();
        for(int it=0;it<nT;it++)
      {
	const int t=
	  it+tMin;
	RFile.write_xy(t,R[it].get());
      }
  }
  
  const auto Afun=
    [&bT](const double& t,const double& r,const PrecFloat& Estar)
  {
    return
      gslIntegrateUpToInfinity([t,r,Estar,&bT](const PrecFloat& E)
      {
	return
	  (sqr(E-Estar)*bT(t,E)*bT(r,E)).get();
      },0);
  };
  
  const auto AfunAn=
    [this](const double& t,const double& r)
  {
	return
	  PrecFloat(1)/(r+t)+
	  PrecFloat(1)/(T+r-t)+
	  PrecFloat(1)/(T-r+t)+
	  PrecFloat(1)/(2*T-r-t);
  };
  
  grace_file_t Cent("/tmp/Cent.xmg");
  grace_file_t Res("/tmp/Res.xmg");
  grace_file_t out("/tmp/E.xmg");
  for(double Estar=0.1;Estar<=5;Estar+=0.1)
    {
      PrecMatr A(nT,nT);
      PrecMatr Aan(nT,nT);
      for(int r=0;r<nT;r++)
	for(int c=0;c<nT;c++)
	  {
	    A(r,c)=Afun(r+tMin,c+tMin,Estar);
	    Aan(r,c)=AfunAn(r+tMin,c+tMin,Estar);
	  }
      
      const auto Ainv=
	A.inverse();
      
      // const MatrixXd Sdiag=
      //   svd.singularValues().asDiagonal();
      
      // const auto AinvBis=
      //   svd.matrixV()*Sdiag.inverse()*svd.matrixU().transpose();
      
      // cout<<(Ainv-AinvBis).eval()<<endl;
      
      const auto num=
	(Ainv*R).eval();
      
      const auto den=
	(R.transpose()*Ainv*R)(0);
      
      grace_file_t gFile("/tmp/g"+to_string(Estar)+".xmg");
      const auto g=
	(num/den).eval();
      for(int iT=0;iT<nT;iT++)
	gFile.write_xy(iT,g[iT].get());
      
      grace_file_t Delta("/tmp/Delta"+to_string(Estar)+".xmg");
      PrecFloat x=0,x2=0,in=0;
      for(double E=0;E<5.0;E+=0.01)
	{
	  PrecFloat s=0;
	  for(int iT=0;iT<nT;iT++)
	    s+=g(iT)*bT(iT,E);
	  Delta.write_xy(E,s.get());
	  
	  x+=E*s*0.01;
	  x2+=E*E*s*0.01;
	  in+=s*0.01;
	}
      x/=in;
      x2/=in;
      x2-=x*x;
      Cent.write_xy(Estar,x.get());
      Res.write_xy(Estar,sqrt(x2).get());
      
      djack_t s{};
      for(int iT=0;iT<nT;iT++)
	for(size_t ijack=0;ijack<=njacks;ijack++)
	  s[ijack]+=(g[iT]*corr[iT+tMin][ijack]).get();
      
      out.write_ave_err(Estar,s.ave_err());
    }
}
