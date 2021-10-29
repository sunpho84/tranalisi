#include <tantaloBacco.hpp>

int main()
{
  PrecFloat::setDefaultPrecision(1024);
  
  set_njacks(10);
  
  const int T=96;
  const int tMin=1;
  const int nT=50;
  const double Efict=2.0;
  const int useBw=0;
  
  djvec_t corr(T);
  grace_file_t corrPlot("corr.xmg");
  for(int t=0;t<T;t++)
    {
      corr[t]=exp(-Efict*t)+exp(-Efict*(T-t))*useBw;
      corrPlot.write_xy(t,corr[t].ave());
    }
  
  const PrecFloat E0=0.1;
  grace_file_t densPlot("dens.xmg");
  grace_file_t sigmaPlot("sigma.xmg");
  PrecFloat sigma,sigma0=0.025;
  for(PrecFloat Estar=E0.get();Estar<4;Estar+=sigma/2)
    {
      const string Etag=to_string(Estar.get());
      grace_file_t RecDelta("RecDelta"+Etag+".xmg");
      grace_file_t TargDelta("TargDelta"+Etag+".xmg");
      grace_file_t ErrDelta("ErrDelta"+Etag+".xmg");
      
      sigma=sigma0*sqrt(Estar/E0);
      sigmaPlot.write_xy(Estar.get(),sigma.get());
      
      const bool useTantalo=true;
      const PrecFloat lambda=0.0;
      const PrecFloat alpha=0.0;
      TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0,lambda,sigma,useBw);
      
      TantaloBaccoRecoEngine recoEngine(pars,Estar);
      TantaloBaccoReco reco(recoEngine,Estar,corr);
      const djack_t d=reco.recoDensity();
      densPlot.write_ave_err(Estar.get(),d.ave_err());
      
      for(PrecFloat E=E0;E<5.0;E+=0.01)
	{
	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
	  const PrecFloat sReco=reco.recoDelta(E);
	  RecDelta.write_xy(E.get(),sReco.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=sReco-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	}
    }
  
  return 0;
}
