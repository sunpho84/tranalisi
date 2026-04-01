// #define FAKE_HP
#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

template <typename T>
T f(const T& aE,const T& t)
{
  return exp(-aE*t)*aE*aE*aE;
}

int main()
{
  PrecFloat::setDefaultPrecision(256);
  
  set_njacks(15);
  
  const int T=192;
  const bool hasBwSignal=false;
  const double a=1;
  const double aE0=0.5*a;
  const double aEPeak1=1.0*a;
  const double aEPeak2=1.5*a;
  djvec_t corr(T/2+1);
  for(int t=0;t<=T/2;t++)
    {
      corr[t]=0.0;
      corr[t]+=f<double>(aEPeak1,t);
      corr[t]+=f<double>(aEPeak2,t);
    }
  
  using namespace Bacco;
  
  CorrelatorPars correlatorPars(T,hasBwSignal,aE0,corr);
  grace_file_t targetDensPlot("/tmp/targetDensPlot.xmg");
  
  const int tMin=1,tMax=60;
  const double lambda=0.5;
  const double sigmaFrEStar=0.1;
  const double mFact=2;
  
  grace_file_t DeviationPlots;
  grace_file_t recoDensPlots;
  
  DeviationPlots.open("/tmp/deviation.xmg");
  recoDensPlots.open("/tmp/recoDensPlot.xmg");
  
  double aEStar=aE0;
  
  do
    {
      const double sigma=sigmaFrEStar*aEStar;
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
      cout<<"aEstar: "<<aEStar<<endl;
      cout<<"Sigma: "<<sigma<<endl;
      
      GenericDivE2Reconstructor gd2Reconstructor(mFact,correlatorPars,tMin,tMax,aEStar,lambda,aE0,
						 [aEStar,sigma](const PrecFloat& E)
						 {
						   return
						     exp(-sqr((E-aEStar)/sigma)/2)/(sqrt(2*precPi())*sigma);
						 });
      const Reconstruction gd2Reco=
	gd2Reconstructor.getReco();
      const double deviation=
	gd2Reconstructor.deviation(gd2Reco);
      cout<<"Deviation: "<<deviation<<endl;
      DeviationPlots.write_xy(aEStar,deviation);
      
      gd2Reco.plot("/tmp/gaussDivE2reco"+to_string(aEStar)+".xmg",aE0,4*aEStar);
      grace_file_t gPlot("/tmp/g"+to_string(aEStar)+".xmg");
      const auto& g=gd2Reco.g;
      const int nT=g.size();
      for(auto i=0;i<nT;i++)
	gPlot.write_xy(i,g(i).get());
      gd2Reconstructor.plot("/tmp/gaussDivE2target"+to_string(aEStar)+".xmg",aE0,4*aEStar);
      
      PrecFloat d=0;
      for(int iT=0;iT<nT;iT++)
	{
	  d+=g(iT)*f<PrecFloat>(aEPeak1,iT+tMin);
	  d+=g(iT)*f<PrecFloat>(aEPeak2,iT+tMin);
	}
      recoDensPlots.write_xy(aEStar,d.get());
      targetDensPlot.write_xy(aEStar,
			      gd2Reconstructor.targetFunction(aEPeak1)*pow(aEPeak1,3)+
			      gd2Reconstructor.targetFunction(aEPeak2)*pow(aEPeak2,3));
      
      aEStar+=sigma/4;
    }
  while(aEStar<4);
  
  return 0;
}
