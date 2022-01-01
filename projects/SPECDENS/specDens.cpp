// #define FAKE_HP
#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

struct fun_minuit_wraper_t :
  public minimizer_fun_t
{
  //! type of function to be passed
  using fun_t=
    function<double(const vector<double> &p)>;
  
  //! function
  fun_t fun;
  
public:
  //! constructor
  fun_minuit_wraper_t(fun_t fun) :
    fun(fun)
  {
  }
  
  //! compute the function
  double operator()(const vector<double> &p) const
  {
    return fun(p);
  }
  
  double Up() const
  {
    return 1;
  }
};

int main()
{
  set_njacks(15);
  PrecFloat::setDefaultPrecision(128);
  
  PrecFloat sigma=0.1;
  PrecFloat Estar=0.05;
  grace_file_t Dens("Dens.xmg");
  grace_file_t SigmaPlot("Sigma.xmg");
  
  const PrecFloat E0=0.0;
  
  const size_t T=192;
  const size_t tMin=1;
  const size_t nT=T/2;
  
  const PrecFloat lambda=0.2;
  const PrecFloat alpha=0.0;
  const int useBw=1;
  constexpr bool useTantalo=false;
  djvec_t corr(T/2+1);
  ifstream input("i");
  if(not input.good())
    CRASH("unable to open input");
  for(size_t t=0;t<=T/2;t++)
    {
      // const PrecFloat f1=exp(-(PrecFloat)t*0.4)+exp(-((PrecFloat)T-t)*0.4);
      // const PrecFloat f2=exp(-(PrecFloat)t*0.8)+exp(-((PrecFloat)T-t)*0.8);
      // corr[t]=f1+f2;
      double f,e;
      input>>f>>e;
      corr[t].fill_gauss({f,e},t+33);
    }
  
  gen_t g(31241);
  do
    {
      cout<<"Estar: "<<Estar.get()<<endl;
      
      TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0,lambda,sigma,useBw);
      TantaloBaccoRecoEngine recoEngine(pars,Estar);
      
      TantaloBaccoReco reco(recoEngine,Estar,corr);
      
      const string Etag=to_string(Estar.get());
      grace_file_t RecDelta("RecDelta"+Etag+".xmg");
      grace_file_t TargDelta("TargDelta"+Etag+".xmg");
      grace_file_t ErrDelta("ErrDelta"+Etag+".xmg");
      
      for(PrecFloat E=0;E<Estar*3;E+=Estar/40)
	{
	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
	  const PrecFloat sReco=reco.recoDelta(E);
	  RecDelta.write_xy(E.get(),sReco.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=sReco-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	}
      
      /// Verify the norm
      // PrecFloat norm1=0.0;
      // for(size_t iT=0;iT<nT;iT++)
      // 	norm1+=reco.R(iT)*reco.g[iT];
      // const PrecFloat norm1m1=norm1-1;
      // cout<<"norm-1: "<<norm1m1.get()<<endl;
      /// Find last coeff
      
      /// Compute the width
      PrecFloat width2=reco.g.transpose()*reco.B*reco.g;
      const PrecFloat width=sqrt(width2);
      cout<<"width: "<<width.get()<<endl;
      
      PrecFloat var=0.0,ave=0.0,norm2=0.0;
      const PrecFloat dE=Estar/160;
      for(PrecFloat E=0;E<Estar*8;E+=dE)
	{
	  const PrecFloat y=reco.recoDelta(E);
	  norm2+=y*y;
	  ave+=y*y*E;
	  var+=y*y*E*E;
	}
      var*=dE;
      ave*=dE;
      norm2*=dE;
      var/=norm2;
      ave/=norm2;
      var-=ave*ave;
      const PrecFloat sigmaEff=sqrt(var);
      cout<<"Ave: "<<ave.get()<<" "<<Estar.get()<<endl;
      cout<<"sigma: "<<sigmaEff.get()<<endl;
      cout<<"norm2: "<<norm2.get()<<endl;
      
      SigmaPlot.write_xy(Estar.get(),width.get());
      
      PrecFloat Err2=0.0;
      for(size_t iT=0;iT<nT;iT++)
	Err2+=sqr(reco.g(iT)*corr[iT+tMin].err());
      const PrecFloat F=(1-lambda)*width2+lambda*Err2;
      cout<<"Functionals: "<<F.get()<<endl;
      
      const djack_t rd=reco.recoDensity();
      Dens.write_ave_err(Estar.get(),rd.ave_err());
      
      grace_file_t gFile("g"+Etag+".xmg");
      for(size_t it=0;it<nT;it++)
	gFile.write_xy(it,reco.g[it].get());
      
      // const PrecFloat EstarPrime=sqr(Estar+sigma)/Estar;
      // const PrecFloat sigmaPrime=sigma*sqrt(EstarPrime/Estar);
      
      //cout<<(sigma+sigmaPrime).get()<<" "<<(EstarPrime-Estar).get()<<endl;
      
      minimizer_pars_t fitPars;
      vector<double> fitG(nT-1);
      for(size_t iT=0;iT<nT-1;iT++)
	fitPars.add("G"+to_string(iT),reco.g(iT).get(),fabs(reco.g(iT).get()/100));

      fun_minuit_wraper_t FCN([&reco,&corr,lambda](const vector<double>& p)
      {
	for(size_t iT=0;iT<nT-1;iT++)
	  reco.g[iT]=p[iT];
	
	PrecFloat lastCoeff=1.0;
	for(size_t iT=0;iT<nT-1;iT++)
	  lastCoeff-=reco.R(iT)*reco.g[iT];
	lastCoeff/=reco.R(nT-1);
	cout<<"lastCoeff: "<<lastCoeff.get()<<endl;
	
	/// Compute the width
	PrecFloat width2=reco.g.transpose()*reco.B*reco.g;
	const PrecFloat width=sqrt(width2);
	cout<<"width: "<<width.get()<<endl;
	
	const djack_t rd=reco.recoDensity();
	const ave_err_t ae=rd.ave_err();
	//(Estar.get(),rd.ave_err());
	
	// PrecFloat Err2=0.0;
	// for(size_t iT=0;iT<nT;iT++)
	//   Err2+=sqr(reco.g(iT)*corr[iT+tMin].err());
	const PrecFloat F=(1-lambda)*width*sqr(ae.ave())+lambda*ae.err();
	cout<<"Functionals: "<<F.get()<<endl;
	
	return F.get();
      });
      
      minimizer_t migrad(FCN,fitPars);
      migrad.minimize();
      
      Estar+=width;
      // Estar+=0.05;//sigmaEff;//EstarPrime;
      //sigma=sigmaPrime;
      
      // cout<<Estar.get()<<" "<<width.get()<<endl;
      cout<<"================="<<endl;
    }
  while(Estar<5);
  
  return 0;
}
