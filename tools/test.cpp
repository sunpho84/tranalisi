#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

int main()
{
  set_njacks(15);
  PrecFloat::setDefaultPrecision(1024);
  
  PrecFloat sigma=0.025;
  PrecFloat Estar=0.1;
  grace_file_t Dens("/tmp/Dens.xmg");
  
  const PrecFloat E0=0.0;
  
  const size_t T=128;
  const size_t tMin=1;
  const size_t nT=T/2;
  
  const PrecFloat lambda=0.0;
  const PrecFloat alpha=0.0;
  const int useBw=1;
  constexpr bool useTantalo=true;
  vector<jack_t<PrecFloat>> corr(T/2+1);
  ifstream input("/tmp/i");
  for(size_t t=0;t<=T/2;t++)
    {
      //const PrecFloat f=exp(-(PrecFloat)t*0.4)+exp(-((PrecFloat)T-t)*0.4);
      double f;
      input>>f;
      corr[t]=f;
    }
  gen_t g(31241);
  do
    {
      TantaloBaccoPars pars(useTantalo,T,tMin,nT,E0,lambda,sigma,useBw);
      TantaloBaccoRecoEngine recoEngine(pars,Estar);
      
      TantaloBaccoReco reco(recoEngine,Estar,corr);
      
      const string Etag=to_string(Estar.get());
      grace_file_t RecDelta("/tmp/RecDelta"+Etag+".xmg");
      grace_file_t RecDelta2("/tmp/RecDelta2"+Etag+".xmg");
      grace_file_t TargDelta("/tmp/TargDelta"+Etag+".xmg");
      grace_file_t ErrDelta("/tmp/ErrDelta"+Etag+".xmg");
      
      PrecFloat y=0;
      for(size_t it=0;it<nT;it++)
	y+=reco.g[it]*exp(-0.4*(PrecFloat)it);
      
      for(PrecFloat E=0.0;E<Estar*3;E+=0.01)
	{
	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
	  const PrecFloat sReco=reco.recoDelta(E);
	  RecDelta.write_xy(E.get(),sReco.get());
	  TargDelta.write_xy(E.get(),sTarg.get());
	  const PrecFloat dS=sReco-sTarg;
	  ErrDelta.write_xy(E.get(),dS.get());
	  
	  PrecFloat sReco2=0.0;
	  for(size_t iT=0;iT<nT;iT++)
	    sReco2+=reco.g(iT)*reco.bT(iT+tMin,E);
	  sReco2+=(reco.g(2)*g.get_double(-1e-5,1e-5))*reco.bT(2+tMin,E);
	  RecDelta2.write_xy(E.get(),sReco2.get());
	}
      
      const PrecFloat rd=reco.recoDensity()[0];
      for(PrecFloat E=Estar-sigma;E<Estar+sigma;E+=0.01)
	Dens.write_xy(E.get(),rd.get());
      
      grace_file_t gFile("/tmp/g"+Etag+".xmg");
      for(size_t it=0;it<nT;it++)
	gFile.write_xy(it,reco.g[it].get());
      
      const PrecFloat EstarPrime=sqr(Estar+sigma)/Estar;
      const PrecFloat sigmaPrime=sigma*sqrt(EstarPrime/Estar);
      
      //cout<<(sigma+sigmaPrime).get()<<" "<<(EstarPrime-Estar).get()<<endl;
      
      Estar=EstarPrime;
      sigma=sigmaPrime;
      
      cout<<Estar.get()<<" "<<sigma.get()<<endl;
    }
  while(Estar<1);
  
  return 0;
}
