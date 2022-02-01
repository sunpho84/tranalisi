 #define FAKE_HP
#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

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
  
  set_njacks(30);
  PrecFloat::setDefaultPrecision(precision);
  
  const double E0=2*0.140*a;
  
  double sigma=0.15*a;
  double Estar=E0;//0.6*a;
  grace_file_t Dens("plots/Dens.xmg");
  grace_file_t SigmaPlot("plots/Sigma.xmg");
  SigmaPlot.set_yaxis_label("\\xs\\0 [GeV]");
  SigmaPlot.set_xaxis_label("E [GeV]");
  grace_file_t RsFile("plots/Rs.xmg");
  RsFile.set_title("R(s)");
  RsFile.set_xaxis_label("E [GeV]");
  grace_file_t RsAltFile("plots/RsAlt.xmg");
  // grace_file_t RsFromPeakAltsFile("plots/RsFromPeaks.xmg");
  
  const PrecFloat alpha=0.0;
  const djvec_t corr=read_djvec(corrPath,T/2+1)*sqr(za)*12*sqr(M_PI);
  
  corr.ave_err().write("plots/Corr.xmg");
  effective_mass(corr,T/2,hasBwSignal?+1:0).ave_err().write("plots/EffMass.xmg");
  
  do
    {
      cout<<"Estar: "<<Estar<<" = "<<Estar/a<<" GeV"<<endl;
      cout<<"Sigma: "<<sigma<<" = "<<sigma/a<<" GeV"<<endl;
      
      using namespace Bacco;
      
      CorrelatorPars correlatorPars(T,hasBwSignal,E0,corr);
      GaussReconstructor reconstructor(correlatorPars,tMin,tMax,Estar,lambda,E0,sigma);
      
      const Reconstruction reco=
	reconstructor.getReco();
      
      GaussDivE2Reconstructor gd2Reconstructor(correlatorPars,tMin,tMax,Estar,lambda,E0,sigma);
      const Reconstruction gd2Reco=
	gd2Reconstructor.getReco();
      // cout<<"Normalization: "<<gd2Reconstructor.normalization()<<endl;
      gd2Reco.plot("plots/gaussDivE2reco"+to_string(Estar)+".xmg",E0,4*Estar);
      gd2Reconstructor.plot("plots/gaussDivE2target"+to_string(Estar)+".xmg",E0,4*Estar);
      grace_file_t plotSmeFWE("plots/gaussDivE2reco"+to_string(Estar)+"WithError.xmg");
      plotSmeFWE.write_polygon([&](const double& E)
      {
	djack_t res;
	for(size_t ijack=0;ijack<=njacks;ijack++)
	    {
	      PrecFloat s=0;
	      
	      for(int iT=0;iT<gd2Reco.nT;iT++)
		s+=gd2Reco.g[iT]*correlatorPars.corr[iT][ijack]/correlatorPars.corr[iT].ave()*correlatorPars.bT(iT+tMin,E);
	      
	      res[ijack]=s.get();
	    }
	return res;
      },E0,4*Estar);
      
      // cout<<"Gauss norm: "<<reco.norm()<<endl;
      // cout<<"Gauss/E^2 norm: "<<gd2Reco.norm()<<endl;
      /// Compute the width
      // const double width=
      // 	reco.widthAssumingGaussianAround(Estar);
      // cout<<"width: "<<width<<endl;
      
      const double norm2Difference=
	reconstructor.deviation(reco);
      cout<<"norm2 of the difference between target and reconstructed: "<<norm2Difference<<endl;
      const double gDivE2norm2Difference2=
	gd2Reconstructor.deviation2(gd2Reco);
      cout<<"norm2 of the difference2 between g/e^2target and reconstructed: "<<gDivE2norm2Difference2<<endl;
      const double gDivE2norm2Difference=
	gd2Reconstructor.deviation(gd2Reco);
      cout<<"norm2 of the difference between g/e^2target and reconstructed: "<<gDivE2norm2Difference<<endl;
      const double myStatisticalError=
	gd2Reco.myStatisticalError().get();
      cout<<"my statistical error: "<<myStatisticalError<<endl;
      
      const double deviation=
	gDivE2norm2Difference2;
      
      // cout<<"deviation of reco from target: "<<widthDeviation<<endl;
      
      const double minDev=1e-4;
      const double maxDev=5e-4;
      if(0 or ( deviation>minDev and
		deviation<maxDev))
	{
	  // cout<<"norm of reconstructed function: "<<reco.norm()<<endl;
	  // cout<<"square norm of reconstructed function: "<<reco.squareNorm()<<endl;
	  // cout<<"projection with reconstructed function: "<<reconstructor.projectionWithReco(reco)<<endl;
	  // cout<<"square norm of target function: "<<reconstructor.squareNorm()<<endl;
	  
	  // const double mean=
	  //   reco.mean().get();
	  // cout<<"mean: "<<mean<<endl;
	  
	  SigmaPlot.write_xy(Estar/a,sigma/a);
	  
	  /// Valid only on B64
	  // PrecFloat s=0;
	  // for(const auto& [ePeak,wPeak] :
	  // 	{std::pair<double,double>{0.540665,0.271642},
	  // 	 {0.687199, 1.73837},
	  // 	 {0.801856, 2.54741},
	  // 	 {0.893164, 1.05354},
	  // 	 {0.998051, 0.201638},
	  // 	 {1.12037 ,0.325096}})
	  //   s+=wPeak*gd2Reco.smearingFunction(ePeak*a)*sqr(ePeak*a);
	  // RsFromPeakAltsFile.write_xy(Estar/a,s.get());
	  
	  PrecFloat Err2=0.0;
	  for(size_t iT=0;iT<nT;iT++)
	    Err2+=sqr(reco.g[iT]*corr[iT+tMin].err());
	  
	  const djack_t rd=
	    reco.recoDensity();
	  Dens.write_ave_err(Estar,rd.ave_err());
	  const djack_t rs=rd/sqr(Estar);
	  RsFile.write_ave_err(Estar/a,rs.ave_err());
	  RsAltFile.write_ave_err(Estar/a,gd2Reco.recoDensity().ave_err());
	  
	  const string Etag=to_string(Estar);
	  reconstructor.plot("plots/TargDelta"+Etag+".xmg",E0,Estar*4);
	  reco.plot("plots/RecoDelta"+Etag+".xmg",E0,Estar*4);
	  // grace_file_t TargDelta("plots/TargDelta"+Etag+".xmg");
	  // grace_file_t ErrDelta("plots/ErrDelta"+Etag+".xmg");
	  
	  // for(double E=0;E<Estar*3;E+=Estar/40)
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
	    gFile.write_xy(it,reco.g[it].get());
	  
	  // const PrecFloat EstarPrime=sqr(Estar+sigma)/Estar;
	  // const PrecFloat sigmaPrime=sigma*sqrt(EstarPrime/Estar);
	  
	  //cout<<(sigma+sigmaPrime)<<" "<<(EstarPrime-Estar)<<endl;
	  
	  // minimizer_pars_t fitPars;
	  // vector<double> fitG(nT-1);
	  // for(size_t iT=0;iT<nT-1;iT++)
	  // 	fitPars.add("G"+to_string(iT),reco.g(iT).get(),fabs(reco.g(iT).get()/100));
	  
	  // fun_minuit_wraper_t FCN([&reco,&corr,nT,lambda](const vector<double>& p)
	  // {
	  // 	for(size_t iT=0;iT<nT-1;iT++)
	  // 	  reco.g[iT]=p[iT];
	  
	  // 	PrecFloat lastCoeff=1.0;
	  // 	for(size_t iT=0;iT<nT-1;iT++)
	  // 	  lastCoeff-=reco.R(iT)*reco.g[iT];
	  // 	lastCoeff/=reco.R(nT-1);
	  // 	cout<<"lastCoeff: "<<lastCoeff<<endl;
	  
	  // 	/// Compute the width
	  // 	PrecFloat width2=reco.g.transpose()*reco.B*reco.g;
	  // 	const PrecFloat width=sqrt(width2);
	  // 	cout<<"width: "<<width<<endl;
	  
	  // 	const djack_t rd=reco.recoDensity();
	  // 	const ave_err_t ae=rd.ave_err();
	  // 	//(Estar.get(),rd.ave_err());
	  
	  // 	// PrecFloat Err2=0.0;
	  // 	// for(size_t iT=0;iT<nT;iT++)
	  // 	//   Err2+=sqr(reco.g(iT)*corr[iT+tMin].err());
	  // 	const PrecFloat F=(1-lambda)*width*sqr(ae.ave())+lambda*ae.err();
	  // 	cout<<"Functionals: "<<F<<endl;
	  
	  // 	return F.get();
	  // });
	  
	  // minimizer_t migrad(FCN,fitPars);
	  // migrad.minimize();
	  
	  
	  const double newEstar=
	    Estar+sigma/4;
	  // const double newSigma=
	  //   width*newEstar/Estar;
	  
	  Estar=newEstar;
	  //sigma=newSigma;
	  
	  cout<<"================="<<endl;
	}
      else
	{
	  cout<<"Gaussian deviating too ";
	  if(deviation<minDev)
	    {
	      cout<<"little, shrinking it"<<endl;
	      sigma/=1.04;
	    }
	  else
	    {
	      cout<<"much, enlarging it"<<endl;
	      sigma*=1.05;
	    }
	}
    }
  while(Estar<1);
  
  return 0;
}
