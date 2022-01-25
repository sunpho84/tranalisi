// #define FAKE_HP
#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

int main()
{
  PrecFloat::setDefaultPrecision(1024);
  cout.precision(PrecFloat::getNDigits());
  
  const PrecFloat r=integrateUpToInfinite([](const PrecFloat& x){return exp(-sqr(x-0.5)/2);});
  cout<<r<<endl;
  const PrecFloat n=sqrt(precPi()/2)*erfc(-(PrecFloat)0.5/sqrt(PrecFloat(2)));
  cout<<n<<endl;
  
  // return 0;
  
  // testQ();

  // return 0;
  
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
  
  set_njacks(30);
  PrecFloat::setDefaultPrecision(precision);
  
  double sigma=0.300*a;
  double Estar=0.1*a;
  grace_file_t Dens("plots/Dens.xmg");
  grace_file_t SigmaPlot("plots/Sigma.xmg");
  grace_file_t RsFile("plots/Rs.xmg");
  RsFile.set_title("R(s)");
  RsFile.set_xaxis_label("E [GeV]");
  
  const double E0=0.0;
  
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
  effective_mass(corr,T/2,hasBwSignal?+1:0).ave_err().write("plots/EffMass.xmg");
  
  // gen_t g(31241);
  do
    {
      cout<<"Estar: "<<Estar<<endl;
      cout<<"Sigma: "<<sigma<<endl;
      
      using namespace Bacco;
      
      CorrelatorPars correlatorPars(T,hasBwSignal,E0,corr);
      GaussReconstructor reconstructor(correlatorPars,tMin,tMax,Estar,lambda,E0,sigma);
      GaussDivE2Reconstructor gd2Reconstructor(correlatorPars,tMin,tMax,Estar,lambda,E0,sigma);
  
      //BGReconstructor reconstructor(correlatorPars,tMin,tMax,Estar,lambda);
      
      const Reconstruction reco=
	reconstructor.getReco();
      
      /// Compute the width
      const double width=
	reco.widthAssumingGaussianAround(Estar);
      cout<<"width: "<<width<<endl;
      
      const double norm2Difference=
	reconstructor.deviation(reco);
      cout<<"norm2 of the difference between target and reconstructed: "<<norm2Difference<<endl;
      
      const double widthDeviation=
	width/sigma-1;
      
      cout<<"deviation of reco from target: "<<widthDeviation<<endl;

      const double lowDev=0.01;
      const double highDev=0.012;
      if(1 or (widthDeviation>lowDev and widthDeviation<highDev))
	{
	  // cout<<"norm of reconstructed function: "<<reco.norm()<<endl;
	  // cout<<"square norm of reconstructed function: "<<reco.squareNorm()<<endl;
	  // cout<<"projection with reconstructed function: "<<reconstructor.projectionWithReco(reco)<<endl;
	  // cout<<"square norm of target function: "<<reconstructor.squareNorm()<<endl;
	  
	  const double mean=
	    reco.mean().get();
	  cout<<"mean: "<<mean<<endl;
	  
	  SigmaPlot.write_xy(Estar/a,width/a);
	  
	  PrecFloat Err2=0.0;
	  for(size_t iT=0;iT<nT;iT++)
	    Err2+=sqr(reco.g[iT]*corr[iT+tMin].err());
	  
	  const djack_t rd=
	    reco.recoDensity();
	  Dens.write_ave_err(Estar,rd.ave_err());
	  const djack_t rs=rd/sqr(Estar);
	  RsFile.write_ave_err(Estar/a,rs.ave_err());
	  
	  const string Etag=to_string(Estar);
	  reconstructor.plot("plots/TargDelta"+Etag+".xmg",0,Estar*4);
	  reco.plot("plots/RecoDelta"+Etag+".xmg",0,Estar*4);
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
	    Estar+width/4;
	  // const double newSigma=
	  //   width*newEstar/Estar;
	  
	  Estar=newEstar;
	  //sigma=newSigma;
	  
	  cout<<"================="<<endl;
	}
      else
	{
	  if(widthDeviation>highDev)
	    sigma*=1.05;
	  if(widthDeviation<lowDev)
	    sigma/=1.048;
	  cout<<"Gaussian deviating too much, enlarging it"<<endl;
	}
    }
  while(Estar<1);
  
  return 0;
}
