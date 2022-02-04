 #define FAKE_HP
#include <tantaloBacco.hpp>

#include <tranalisi.hpp>
#include <minimizer_wrapper.hpp>

int main()
{
  input_file_t input("input.txt");
  const size_t precision=input.read<size_t>("Precision");
  const size_t T=input.read<size_t>("T");
  const size_t tMin=input.read<size_t>("tMin");
  const size_t tMax=input.read<size_t>("tMax");
  const size_t nT=tMax-tMin;
  cout<<"nT: "<<nT<<endl;
  const int hasBwSignal=input.read<int>("hasBwSignal");
  const string corrPath=input.read<string>("CorrPath");
  const double lambda=input.read<double>("Lambda");
  const double a=input.read<double>("a");
  const double za=input.read<double>("za");
  const double sigmaRatio=input.read<double>("SigmaFracEstar");
  const double E0InGeV=input.read<double>("E0InGev");
  const double EMaxInGeV=input.read<double>("EMaxInGeV");
  const int mFact=input.read<int>("MFact");
  
  set_njacks(30);
  const double aE0=E0InGeV*a;
  PrecFloat::setDefaultPrecision(precision);
  
  // grace_file_t Dens("plots/Dens.xmg");
  
  grace_file_t SigmaPlot("plots/Sigma.xmg");
  SigmaPlot.set_yaxis_label("\\xs\\0 [GeV]");
  SigmaPlot.set_xaxis_label("E [GeV]");
  
  grace_file_t DeviationPlot("plots/deviation.xmg");
  DeviationPlot.set_xaxis_label("E [GeV]");
  
  // grace_file_t RsFile("plots/Rs.xmg");
  // RsFile.set_title("R(s)");
  // RsFile.set_xaxis_label("E [GeV]");
  
  grace_file_t RsAltFile("plots/RsAlt.xmg");
  // grace_file_t RsFromPeakAltsFile("plots/RsFromPeaks.xmg");
  
  const PrecFloat alpha=0.0;
  const djvec_t corr=read_djvec(corrPath,T/2+1)*sqr(za)*12*sqr(M_PI);
  
  corr.ave_err().write("plots/Corr.xmg");
  effective_mass(corr,T/2,hasBwSignal?+1:0).ave_err().write("plots/EffMass.xmg");
  
  double aEStar=E0InGeV*a;
  do
    {
      const string Etag=
	to_string(aEStar);
      const double sigma=
	aEStar*sigmaRatio;
      
      cout<<"Estar: "<<aEStar<<" = "<<aEStar/a<<" GeV"<<endl;
      cout<<"Sigma: "<<sigma<<" = "<<sigma/a<<" GeV"<<endl;
      
      using namespace Bacco;
      
      CorrelatorPars correlatorPars(T,hasBwSignal,aE0,corr);
      // GaussReconstructor reconstructor(correlatorPars,tMin,tMax,aEStar,lambda,aE0,sigma);
      
      // const Reconstruction reco=
      // 	reconstructor.getReco();
      
      GenericDivE2Reconstructor gd2Reconstructor(mFact,correlatorPars,tMin,tMax,aEStar,lambda,aE0,
						 [aEStar,sigma](const PrecFloat& E)
						 {
						   return
						     exp(-sqr((E-aEStar)/sigma)/2)/(sqrt(2*precPi())*sigma);
						     // sigma/(sqr(E-aEStar)+sqr(sigma))/M_PI;
						     // ((E>aEStar-sigma/2) and (E<aEStar+sigma/2))/sigma;
						 });
      const Reconstruction gd2Reco=
	gd2Reconstructor.getReco();
      
      // cout<<"Normalization: "<<gd2Reconstructor.normalization()<<endl;
      gd2Reco.plot("plots/gaussDivE2reco"+to_string(aEStar)+".xmg",aE0,4*aEStar);
      gd2Reconstructor.plot("plots/gaussDivE2target"+to_string(aEStar)+".xmg",aE0,4*aEStar);
      grace_file_t gRecoPlot("plots/gaussreco"+to_string(aEStar)+".xmg");
      // gRecoPlot.write_line([&gd2Reconstructor](const double& E){return gd2Reconstructor.targetFunction(E)*E*E;},aE0,4*aEStar);
      // gRecoPlot.write_line([&gd2Reco](const double& E){return gd2Reco.smearingFunction(E).get()*E*E;},aE0,4*aEStar);
      gRecoPlot.write_line([&gd2Reconstructor,sigma,aEStar](const double& x)
      {
	const double E=x*sigma+aEStar;
	return sigma*gd2Reconstructor.targetFunction(E)*E*E;
      },(aE0-aEStar)/sigma,3,grace::GREEN4);
      gRecoPlot.write_line([&gd2Reco,sigma,aEStar](const double& x)
      {
	const double E=x*sigma+aEStar;
	return sigma*gd2Reco.smearingFunction(E).get()*E*E;
      },(aE0-aEStar)/sigma,3);
      gRecoPlot.write_polygon([&](const double& x)
      {
	djvec_t fCorr=corr;
	for(size_t iT=0;iT<fCorr.size();iT++)
	  fCorr[iT].fill_gauss(fCorr[iT].ave_err(),34634+iT);
	const double E=x*sigma+aEStar;
	djack_t f(0.0);
	for(size_t ijack=0;ijack<=njacks;ijack++)
	  {
	    PrecFloat pf=0.0;
	    for(size_t iT=0;iT<nT;iT++)
	      {
		const int t=iT+tMin;
		const djack_t &c=fCorr[t];
		const double r=c[ijack]/c.ave()-1;
		const double fluct=r*sqrt(lambda)+1;
		const PrecFloat& gt=gd2Reco.g[iT];
		pf+=gt*correlatorPars.bT(t,E)*fluct;
	      }
	    f[ijack]=pf.get()*E*E*sigma;
	  }
	return f;
      },(aE0-aEStar)/sigma,3);
      
      PrecFloat A=gslIntegrateUpToInfinity([&gd2Reconstructor,mFact](const PrecFloat& E)
      {
	return (sqr(gd2Reconstructor.preciseTargetFunction(E))*pow(E,mFact*2)).get();
      },aE0);
      PrecFloat B=0.0;
      
      correlatorPars.loopOverIrIt([&A,&B,&tMin,&g=gd2Reco.g,&corr,&gd2Reconstructor](const int& iR,
									       const int& iT,
									       const int& i)
      {
	const PrecFloat a=gd2Reconstructor.aFun(i);
	A+=g[iR]*a*g[iT];
	if(iR==iT)
	  B+=g[iR]*a*g[iT]*sqr(corr[iR+tMin].err()/corr[iR+tMin].ave());
	
      },tMin,nT);
      
      for(size_t iT=0;iT<nT;iT++)
	{
	  PrecFloat f=gd2Reconstructor.fFun(iT+tMin);
	  
	  if(hasBwSignal)
	    f+=
	      gd2Reconstructor.fFun(T-iT-tMin);
	  
	  A-=2*gd2Reco.g[iT]*f;
	}
      
      cout<<" "<<A<<" "<<B<<endl;
      
      // grace_file_t plotSmeFWE("plots/gaussDivE2reco"+to_string(aEStar)+"WithError.xmg");
      // plotSmeFWE.write_polygon([&](const double& E)
      // {
      // 	djack_t res;
      // 	for(size_t ijack=0;ijack<=njacks;ijack++)
      // 	  {
      // 	    PrecFloat s=0;
	    
      // 	    for(int iT=0;iT<gd2Reco.nT;iT++)
      // 	      s+=gd2Reco.g[iT]*correlatorPars.corr[iT][ijack]/correlatorPars.corr[iT].ave()*correlatorPars.bT(iT+tMin,E);
	    
      // 	    res[ijack]=s.get();
      // 	  }
      // 	return res;
      // },aE0,4*aEStar);
      
      // cout<<"Gauss norm: "<<reco.norm()<<endl;
      // cout<<"Gauss/E^2 norm: "<<gd2Reco.norm()<<endl;
      /// Compute the width
      // const double width=
      // 	reco.widthAssumingGaussianAround(aEStar);
      // cout<<"width: "<<width<<endl;
      
      // const double norm2Difference=
      // 	reconstructor.deviation(reco);
      // cout<<"norm2 of the difference between target and reconstructed: "<<norm2Difference<<endl;
      // const double gDivE2norm2Difference2=
      // 	gd2Reconstructor.deviation2(gd2Reco);
      // cout<<"norm2 of the difference2 between g/e^2target and reconstructed: "<<gDivE2norm2Difference2<<endl;
      const double gDivE2norm2Difference=
	gd2Reconstructor.deviation(gd2Reco);
      cout<<"norm2 of the difference between g/e^2target and reconstructed: "<<gDivE2norm2Difference<<endl;
      DeviationPlot.write_xy(aEStar,gDivE2norm2Difference);
      
      // const double myStatisticalError=
      // 	gd2Reco.myStatisticalError().get();
      // cout<<"my statistical error: "<<myStatisticalError<<endl;
      
      // const double deviation=
      // 	gDivE2norm2Difference2;
      
      // cout<<"deviation of reco from target: "<<widthDeviation<<endl;
      
      // const double minDev=1e-4;
      // const double maxDev=5e-4;
      // if(0 or ( deviation>minDev and
      // 		deviation<maxDev))
      // 	{
      // cout<<"norm of reconstructed function: "<<reco.norm()<<endl;
      // cout<<"square norm of reconstructed function: "<<reco.squareNorm()<<endl;
      // cout<<"projection with reconstructed function: "<<reconstructor.projectionWithReco(reco)<<endl;
      // cout<<"square norm of target function: "<<reconstructor.squareNorm()<<endl;
      
      // const double mean=
      //   reco.mean().get();
      // cout<<"mean: "<<mean<<endl;
      
      SigmaPlot.write_xy(aEStar/a,sigma/a);
      
      // Valid only on B64
      // PrecFloat s=0;
      // for(const auto& [ePeak,wPeak] :
      // 	{std::pair<double,double>{0.540665,0.271642},
      // 	 {0.687199, 1.73837},
      // 	 {0.801856, 2.54741},
      // 	 {0.893164, 1.05354},
      // 	 {0.998051, 0.201638},
      // 	 {1.12037 ,0.325096}})
      //   s+=wPeak*gd2Reco.smearingFunction(ePeak*a)*sqr(ePeak*a);
      //RsFromPeakAltsFile.write_xy(aEStar/a,s.get());
      
      // PrecFloat Err2=0.0;
      // for(size_t iT=0;iT<nT;iT++)
      // 	Err2+=sqr(reco.g[iT]*corr[iT+tMin].err());
      
      // const djack_t rd=
      // 	reco.recoDensity();
      // Dens.write_ave_err(aEStar,rd.ave_err());
      // const djack_t rs=rd/sqr(aEStar);
      // RsFile.write_ave_err(aEStar/a,rs.ave_err());
      RsAltFile.write_ave_err(aEStar/a,gd2Reco.recoDensity().ave_err());
      
      // reconstructor.plot("plots/TargDelta"+Etag+".xmg",aE0,aEStar*4);
      // reco.plot("plots/RecoDelta"+Etag+".xmg",aE0,aEStar*4);
      // grace_file_t TargDelta("plots/TargDelta"+Etag+".xmg");
      // grace_file_t ErrDelta("plots/ErrDelta"+Etag+".xmg");
      
      // for(double E=0;E<aEStar*3;E+=aEStar/40)
      // 	{
      // 	  const PrecFloat sTarg=
      // 	    gaussReconstructor.targetFunction(E);
      // 	  const PrecFloat sReco=
      // 	    gaussReco.smearingFunction(E);
      // 	  RecDelta.write_xy(E,sReco.get());
      // 	  TargDelta.write_xy(E,sTarg.get());
      // 	  const PrecFloat dS=sReco-sTarg;
      // 	  ErrDelta.write_xy(E,dS.get());
      // 	}
      
      grace_file_t gFile("plots/g"+Etag+".xmg");
      for(size_t it=0;it<nT;it++)
	gFile.write_xy(it,gd2Reco.g[it].get());
      
      // const PrecFloat aEStarPrime=sqr(aEStar+sigma)/aEStar;
      // const PrecFloat sigmaPrime=sigma*sqrt(aEStarPrime/aEStar);
      
      //cout<<(sigma+sigmaPrime)<<" "<<(aEStarPrime-aEStar)<<endl;
      
      // minimizer_pars_t fitPars;
      // vector<double> fitG(nT);
      // for(size_t iT=0;iT<nT;iT++)
      // 	fitPars.add("G"+to_string(iT),gd2Reco.g(iT).get(),fabs(gd2Reco.g(iT).get()/100));
      
      // djvec_t fCorr=corr;
      // for(size_t iT=0;iT<fCorr.size();iT++)
      // 	fCorr[iT].fill_gauss(fCorr[iT].ave_err(),34634+iT);
      
      // fun_minuit_wraper_t FCN([&](const vector<double>& p)
      // {
      // 	const double F=
      // 	  gslIntegrateUpToInfinity([&](const double& E)
      // 	  {
      // 	    const double t=gd2Reconstructor.targetFunction(E);
      // 	    double e=0;
      // 	    for(size_t ijack=0;ijack<njacks;ijack++)
      // 	      {
      // 		PrecFloat pf=0.0;
      // 		for(size_t iT=0;iT<nT;iT++)
      // 		  {
      // 		    const int t=iT+tMin;
      // 		    const djack_t &c=fCorr[t];
      // 		    const double r=c[ijack]/c.ave()-1;
      // 		    const double fluct=r*sqrt(lambda)+1;
      // 		    const PrecFloat& gt=p[iT];
      // 		    pf+=gt*correlatorPars.bT(t,E)*fluct;
      // 		  }
      // 		e+=sqr((pf.get()-t)*E*E);
      // 	      }
      // 	    return e;
      // 	  },aE0);
	
      // 	cout<<"Functionals: "<<F<<endl;
	
      // 	return F;
      // });
      
      // minimizer_t migrad(FCN,fitPars);
      // vector<double> pars=migrad.minimize();
      
      // /// Parameters of the reconstrucion
      // PrecVect g(nT);
      // for(size_t iT=0;iT<nT;iT++)
      // 	g[iT]=pars[iT];
      
      // Reconstruction fitReco(nT,tMin,correlatorPars,g);
      // gRecoPlot.write_line([&fitReco,sigma,aEStar](const double& x)
      // {
      // 	const double E=x*sigma+aEStar;
      // 	return sigma*fitReco.smearingFunction(E).get()*E*E;
      // },(aE0-aEStar)/sigma,3);
      // CRASH("");
      aEStar+=sigma/4;
      // const double newSigma=
      //   width*newaEStar/aEStar;
      
      
      cout<<"================="<<endl;
      // }
      // else
      //   {
      //     cout<<"Gaussian deviating too ";
      //     if(deviation<minDev)
      // 	{
      // 	  cout<<"little, shrinking it"<<endl;
      // 	  sigma/=1.04;
      // 	}
      //     else
      // 	{
      // 	      cout<<"much, enlarging it"<<endl;
      // 	      sigma*=1.05;
      // 	    }
      // 	}
    }
  while(aEStar/a<EMaxInGeV);
  
  return 0;
}
