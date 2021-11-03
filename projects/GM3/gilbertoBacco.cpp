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
  PrecFloat::setDefaultPrecision(1024);
  
  // Becchi's delta f(x)=(e2+x2)/(2*(e2-x2)^2)/cosh^2(x/(e2-x2)) eq 270 pag.90
  
  const double eu=2.0/3,ed=-1.0/3;
  vector<jack_t<PrecFloat>> corr(THp1);
  
  {
    djvec_t _corr=
      getAveForRego(0,nSources,1,REGO_TM)*sqr(Z[regoZId[REGO_TM]])*(sqr(eu)+sqr(ed));
    for(size_t it=0;it<=TH;it++)
      corr[it]=_corr[it];
  }
   
   // djack_t Z2(0.01);
   // djack_t M(1.4);
   // djack_t e(gauss_filler_t{1,0.1,235253});
   // for(size_t t=0;t<=TH;t++)
   //   corr[t]=two_pts_corr_fun(Z2,M,TH,t,+1*0)*e;
  
  const int nT=35,tMin=1;
  
  /////////////////////////////////////////////////////////////////
  
  const PrecFloat E0=0.05;
  grace_file_t out("/tmp/dens.xmg");
  PrecFloat sigma;
  for(PrecFloat Estar=E0.get();Estar<4;Estar+=sigma/2)
    {
      sigma=0.05*sqrt(Estar/E0);
      
      const string Etag=to_string(Estar.get());
      const PrecFloat lambda=0.0;
      const PrecFloat alpha=0.0;
      const int useBw=1;
      constexpr bool useTantalo=true;
      TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0,lambda,sigma,useBw);
      TantaloBaccoRecoEngine recoEngine(pars,Estar);
      
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
      
      grace_file_t RecDelta("/tmp/RecDelta"+Etag+".xmg");
      grace_file_t TargDelta("/tmp/TargDelta"+Etag+".xmg");
      grace_file_t ErrDelta("/tmp/ErrDelta"+Etag+".xmg");
      // PrecFloat x=0,x2=0,in=0;
      for(PrecFloat E=E0;E<5.0;E+=0.01)
	{
	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
	  const PrecFloat sReco=reco.recoDelta(E);
	  RecDelta.write_xy(E.get(),sReco.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=sReco-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	  
	  // x+=E*sReco*0.01;
	  // x2+=E*E*sReco*0.01;
	  // in+=sReco*0.01;
	}
      // x/=in;
      // x2/=in;
      // x2-=x*x;
      // Cent.write_xy(Estar,x.get());
      // Res.write_xy(Estar,sqrt(x2).get());
      
      djack_t d;
      {
	const jack_t<PrecFloat> _d=reco.recoDensity();
	for(size_t ijack=0;ijack<=njacks;ijack++)
	  d[ijack]=_d[ijack].get();
      }
      out.write_ave_err(Estar.get(),d.ave_err());
    }
  
  cout<<"Gilbert done"<<endl;
}
