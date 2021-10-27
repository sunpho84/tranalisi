#include <tantaloBacco.hpp>

#include <perens.hpp>


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
  
  // Becchi's delta f(x)=(e2+x2)/(2*(e2-x2)^2)/cosh^2(x/(e2-x2)) eq 270 pag.90
  
  const double eu=2.0/3,ed=-1.0/3;
   djvec_t corr=
     getAveForRego(0,nSources,1,REGO_TM)*sqr(Z[regoZId[REGO_TM]])*(sqr(eu)+sqr(ed));
   
  // djack_t Z2(gauss_filler_t{0.1,0.01,234324});
  // djack_t M=djack_t(gauss_filler_t{1.4,0.01,235253});
  
  // for(size_t t=0;t<=TH;t++)
  //   corr[t]=two_pts_corr_fun(Z2,M,TH,t,+1);
  
  
  const int nT=30,tMin=1;
  
  /////////////////////////////////////////////////////////////////
  
  const PrecFloat E0=0.1;
  const PrecFloat lambda=0; 
  const PrecFloat sigma=0.1;
  const PrecFloat alpha=0.0;
  const int useBw=1;
  
  /////////////////////////////////////////////////////////////////
  
  TantaloBaccoPars pars(T,tMin,nT,E0,lambda,sigma,useBw);
  TantaloBaccoRecoEngine recoEngine(pars);
  
  grace_file_t out("/tmp/E.xmg");
  for(double Estar=E0.get();Estar<4;Estar+=sigma.get()/2)
    {
      TantaloBaccoReco reco(recoEngine,Estar,corr);
      
      // const MatrixXd Sdiag=
      //   svd.singularValues().asDiagonal();
      
      // const auto AinvBis=
      //   svd.matrixV()*Sdiag.inverse()*svd.matrixU().transpose();
      
      // cout<<(Ainv-AinvBis).eval()<<endl;
      
      // const auto num=
      // 	(Ainv*R).eval();
      
      // const auto den=
      // 	(R.transpose()*Ainv*R)(0);
      
      grace_file_t RecDelta("/tmp/RecDelta"+to_string(Estar)+".xmg");
      grace_file_t TargDelta("/tmp/TargDelta"+to_string(Estar)+".xmg");
      grace_file_t ErrDelta("/tmp/ErrDelta"+to_string(Estar)+".xmg");
      PrecFloat x=0,x2=0,in=0;
      for(PrecFloat E=E0;E<2.0;E+=0.01)
	{
	  
	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
	  const PrecFloat sReco=reco.recoDelta(E);
	  RecDelta.write_xy(E.get(),sReco.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=sReco-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	  
	  x+=E*sReco*0.01;
	  x2+=E*E*sReco*0.01;
	  in+=sReco*0.01;
	}
      x/=in;
      x2/=in;
      x2-=x*x;
      // Cent.write_xy(Estar,x.get());
      // Res.write_xy(Estar,sqrt(x2).get());
      
      const djack_t d=reco.recoDensity();
      out.write_ave_err(Estar,d.ave_err());
    }
  
  cout<<"Gilbert done"<<endl;
}
