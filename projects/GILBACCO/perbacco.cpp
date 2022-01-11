#include <tantaloBacco.hpp>

int main()
{
  PrecFloat::setDefaultPrecision(1024);
  
  set_njacks(10);
  
  const int T=96;
  const int tMin=1;
  const int nT=40;
  const PrecFloat Efict=1.5;
  const int useBw=0;
  
  djvec_t corr(T);
  grace_file_t corrPlot("corr.xmg");
  for(int t=0;t<T;t++)
    {
      corr[t]=(exp(-Efict*t)+exp(-Efict*(T-t))*useBw).get();
      corrPlot.write_xy(t,corr[t].ave());
    }
  
  const PrecFloat E0=0.1;
  grace_file_t densPlot("dens.xmg");
  grace_file_t sigmaPlot("sigma.xmg");
  grace_file_t centPlot("cent.xmg");
  grace_file_t resPlot("res.xmg");
  grace_file_t convPlot("conv.xmg");
  PrecFloat sigma,sigma0=0.2;
  for(PrecFloat Estar=E0.get();Estar<4;Estar+=sigma/4)
    {
      cout<<Estar.get()<<endl;
      
      const string Etag=to_string(Estar.get());
      grace_file_t RecDelta("RecDelta"+Etag+".xmg");
      grace_file_t TargDelta("TargDelta"+Etag+".xmg");
      grace_file_t ErrDelta("ErrDelta"+Etag+".xmg");
      
      sigma=sigma0*sqrt(Estar/E0);
      sigmaPlot.write_xy(Estar.get(),sigma.get());
      
      const bool useTantalo=false;
      const PrecFloat lambda=0.0;
      const PrecFloat alpha=0.0;
      TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0*0,lambda,sigma,useBw);
      
      TantaloBaccoRecoEngine recoEngine(pars,Estar);
      TantaloBaccoReco reco(recoEngine,Estar,corr);
      djack_t d=reco.recoDensity();
      densPlot.write_ave_err(Estar.get(),d.ave_err());
      
      PrecFloat x=0,x2=0,in=0;
      for(PrecFloat E=0;E<5.0;E+=0.01)
	{
	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
	  const PrecFloat sReco=reco.recoDelta(E);
	  RecDelta.write_xy(E.get(),sReco.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=sReco-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	  
	  const PrecFloat w=abs(sReco);
	  x+=E*w;
	  x2+=E*E*w;
	  in+=w;
	}
      x/=in;
      x2/=in;
      x2-=x*x;
      
      const double cen=x.get();
      const double res=sqrt(x2).get();
      centPlot.write_xy(Estar.get(),cen);
      resPlot.write_xy(Estar.get(),res);
      
      const PrecFloat sConv=reco.recoDelta(Efict);
      convPlot.write_xy(Estar.get(),sConv.get());
      sigma=res;
    }
  
  return 0;
}
