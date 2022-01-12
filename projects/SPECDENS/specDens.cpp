// #define FAKE_HP
#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

// constexpr bool useTantalo=false;

int main()
{
  input_file_t input("input.txt");
  const size_t precision=input.read<size_t>("Precision");
  const size_t T=input.read<size_t>("T");
  const size_t tMin=input.read<size_t>("tMin");
  const size_t tMax=input.read<size_t>("tMax");
  const size_t nT=tMax-tMin;
  cout<<"nT: "<<nT<<endl;
  const int useBw=input.read<int>("useBw");
  const string corrPath=input.read<string>("CorrPath");
  const PrecFloat lambda=input.read<double>("Lambda");
  const double a=input.read<double>("a");
  
  set_njacks(30);
  PrecFloat::setDefaultPrecision(precision);
  
  PrecFloat sigma=0.1;
  PrecFloat Estar=0.05;
  grace_file_t Dens("plots/Dens.xmg");
  grace_file_t SigmaPlot("plots/Sigma.xmg");
  grace_file_t RsFile("plots/Rs.xmg");
  RsFile.set_title("R(s)");
  RsFile.set_xaxis_label("E [GeV]");
  
  const PrecFloat E0=0.0;
  
  const PrecFloat alpha=0.0;
  const djvec_t corr=read_djvec(corrPath,T/2+1);
  // raw_file_t inputCorr(corrPath,"r");
  // for(size_t t=0;t<=T/2;t++)
  //   {
  //     const double f=inputCorr.read<double>();
  //     const double e=inputCorr.read<double>();
  //     corr[t].fill_gauss({f,e},t+33);
  //   }
  
  corr.ave_err().write("plots/Corr.xmg");
  effective_mass(corr,T/2,useBw?+1:0).ave_err().write("plots/EffMass.xmg");
  
  // // gen_t g(31241);
  // do
  //   {
  //     cout<<"Estar: "<<Estar.get()<<endl;
      
  //     // double lambda=1/(sqr(corr[1].ave())+1);
  //     // for(int i=0;i<1;i++)
  //     // 	{
  //     // 	  const double eps=1e-3;
  //     // 	  djack_t d[2];
  //     // 	  TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0,lambda,sigma,useBw);
  //     // 	  PrecFloat width2;
	  
  //     // 	  PrecFloat Err2=0.0;
	  
  //     // 	  for(int iEps=0;iEps<=1;iEps++)
  //     // 	    {
  //     // 	      TantaloBaccoRecoEngine recoEngine(pars,Estar+eps*iEps);
  //     // 	      TantaloBaccoReco reco(recoEngine,Estar+eps*iEps,corr);
  //     // 	      d[iEps]=reco.recoDensity();
  //     // 	      cout<<iEps<<" "<<d[iEps]<<endl;
  //     // 	      if(iEps==0)
  //     // 		{
  //     // 		  width2=reco.g.transpose()*reco.A*reco.g;
  //     // 		  for(size_t iT=0;iT<nT;iT++)
  //     // 		    Err2+=sqr(reco.g(iT)*corr[iT+tMin].err());
  //     // 		}
  //     // 	    }
  //     // 	  const djack_t bd=(d[1]-d[0])/eps;
  //     // 	  // const djack_t fd=(d[2]-d[1])/eps;
  //     // 	  cout<<"Der: "<<bd.ave_err()<<endl;
  //     // 	  // const PrecFloat A=width2*sqr(bd).ave();
  //     // 	  lambda=1/(sqr(bd)+1).ave();
  //     // 	  cout<<"Lambda: "<<lambda<<endl;
  //     // 	}
  //     //lambda=0.5;
      
  //     TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0,lambda,sigma,useBw);
  //     TantaloBaccoRecoEngine recoEngine(pars,Estar);
  //     TantaloBaccoReco reco(recoEngine,Estar,corr);
      
  //     /// Verify the norm
  //     // PrecFloat norm1=0.0;
  //     // for(size_t iT=0;iT<nT;iT++)
  //     // 	norm1+=reco.R(iT)*reco.g[iT];
  //     // const PrecFloat norm1m1=norm1-1;
  //     // cout<<"norm-1: "<<norm1m1.get()<<endl;
  //     /// Find last coeff
      
  //     /// Compute the width
  //     PrecFloat width2=reco.A.formWith(reco.g,reco.g);
  //     PrecFloat ave=0.0;
  //     for(size_t i=0;i<nT;i++)
  // 	ave+=reco.g[i]*reco.M[i];
  //     const PrecFloat width=sqrt(width2);
  //     cout<<"width: "<<width.get()<<endl;
  //     cout<<"ave: "<<ave.get()<<endl;
      
  //     // PrecFloat var=0.0,ave=0.0;
  //     // const PrecFloat dE=Estar/160;
  //     // for(PrecFloat E=0;E<Estar*8;E+=dE)
  //     // 	{
  //     // 	  const PrecFloat y=reco.recoDelta(E);
  //     // 	  ave+=y*E;
  //     // 	  var+=y*y*sqr(E-Estar);
  //     // 	}
  //     // var*=dE;
  //     // ave*=dE;
  //     // // var-=ave*ave;
  //     // const PrecFloat sigmaEff=sqrt(var);
  //     // cout<<"Ave: "<<ave.get()<<" "<<Estar.get()<<endl;
  //     // cout<<"sigma: "<<sigmaEff.get()<<endl;
      
  //     SigmaPlot.write_xy(Estar.get(),width.get());
      
  //     PrecFloat Err2=0.0;
  //     for(size_t iT=0;iT<nT;iT++)
  // 	Err2+=sqr(reco.g[iT]*corr[iT+tMin].err());
  //     // const PrecFloat F=(1-lambda)*width2+lambda*Err2;
  //     // cout<<"Functionals: "<<F.get()<<endl;
      
  //     const djack_t rd=reco.recoDensity();
  //     Dens.write_ave_err(Estar.get(),rd.ave_err());
  //     const djack_t rs=rd/sqr(Estar.get());
  //     RsFile.write_ave_err(Estar.get()/a,rs.ave_err());
      
  //     const string Etag=to_string(Estar.get());
  //     grace_file_t RecDelta("plots/RecDelta"+Etag+".xmg");
  //     grace_file_t TargDelta("plots/TargDelta"+Etag+".xmg");
  //     grace_file_t ErrDelta("plots/ErrDelta"+Etag+".xmg");
      
  //     for(PrecFloat E=0;E<Estar*3;E+=Estar/40)
  // 	{
  // 	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
  // 	  const PrecFloat sReco=reco.recoDelta(E);
  // 	  RecDelta.write_xy(E.get(),sReco.get());
  // 	  TargDelta.write_xy(E.get(),sTarg.get());
  // 	  const PrecFloat dS=sReco-sTarg;
  // 	  ErrDelta.write_xy(E.get(),dS.get());
  // 	}
      
  //     grace_file_t gFile("plots/g"+Etag+".xmg");
  //     for(size_t it=0;it<nT;it++)
  // 	gFile.write_xy(it,reco.g[it].get());
      
  //     // const PrecFloat EstarPrime=sqr(Estar+sigma)/Estar;
  //     // const PrecFloat sigmaPrime=sigma*sqrt(EstarPrime/Estar);
      
  //     //cout<<(sigma+sigmaPrime).get()<<" "<<(EstarPrime-Estar).get()<<endl;
      
  //     // minimizer_pars_t fitPars;
  //     // vector<double> fitG(nT-1);
  //     // for(size_t iT=0;iT<nT-1;iT++)
  //     // 	fitPars.add("G"+to_string(iT),reco.g(iT).get(),fabs(reco.g(iT).get()/100));

  //     // fun_minuit_wraper_t FCN([&reco,&corr,nT,lambda](const vector<double>& p)
  //     // {
  //     // 	for(size_t iT=0;iT<nT-1;iT++)
  //     // 	  reco.g[iT]=p[iT];
	
  //     // 	PrecFloat lastCoeff=1.0;
  //     // 	for(size_t iT=0;iT<nT-1;iT++)
  //     // 	  lastCoeff-=reco.R(iT)*reco.g[iT];
  //     // 	lastCoeff/=reco.R(nT-1);
  //     // 	cout<<"lastCoeff: "<<lastCoeff.get()<<endl;
	
  //     // 	/// Compute the width
  //     // 	PrecFloat width2=reco.g.transpose()*reco.B*reco.g;
  //     // 	const PrecFloat width=sqrt(width2);
  //     // 	cout<<"width: "<<width.get()<<endl;
	
  //     // 	const djack_t rd=reco.recoDensity();
  //     // 	const ave_err_t ae=rd.ave_err();
  //     // 	//(Estar.get(),rd.ave_err());
	
  //     // 	// PrecFloat Err2=0.0;
  //     // 	// for(size_t iT=0;iT<nT;iT++)
  //     // 	//   Err2+=sqr(reco.g(iT)*corr[iT+tMin].err());
  //     // 	const PrecFloat F=(1-lambda)*width*sqr(ae.ave())+lambda*ae.err();
  //     // 	cout<<"Functionals: "<<F.get()<<endl;
	
  //     // 	return F.get();
  //     // });
      
  //     // minimizer_t migrad(FCN,fitPars);
  //     // migrad.minimize();
      
  //     Estar+=width;///4;//10;
  //     // Estar+=0.05;//sigmaEff;//EstarPrime;
  //     //sigma=sigmaPrime;
      
  //     // cout<<Estar.get()<<" "<<width.get()<<endl;
  //     cout<<"================="<<endl;
  //   }
  // while(Estar<1);
  
  return 0;
}
