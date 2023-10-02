#include <tranalisi.hpp>

#include <loadCorrCov.hpp>

#include <trambacco.hpp>

const size_t T=224;

/// Loads a correlator avoiding to compute covariance, dropping returned dummy vector

inline double thetaSigma(const double& sigma,
			 const double& x)
{
  return 1/(1+exp(-x/sigma));
}

inline double Kl(const double& sigma,
		 const double& x)
{
  return sqr(1-x*x)*thetaSigma(sigma,1-x)/x;
}

inline double Kt(const double& sigma,
		 const double& x)
{
  return (1+2*x*x)*Kl(sigma,x);
}

int main()
{
  set_njacks(50);
  
  const string ensName=
    "C.06.112";
  
  const string baseOut=
    "/home/francesco/QCD/LAVORI/TAU_VUS/";
  
  auto loadCorrCov=
    [&ensName,&baseOut](const string& chann)
    {
      const string baseIn=
	"/home/francesco/QCD/LAVORI/GM3/";
      
      const auto [readT,corr,cov]=
	::loadCorrCov(baseIn+"/"+ensName+"/data/"+chann,false);
      
      if(T!=readT)
	CRASH("T=%zu, read %zu",T,readT);
      
      corr.ave_err().write(baseOut+ensName+"/"+chann+".xmg");
      effective_mass(corr).ave_err().write(baseOut+ensName+"/eff_mass_"+chann+".xmg");
      
      return make_tuple(corr,cov);
    };
  
    auto loadCorr=
    [&loadCorrCov](const string& chann)
    {
      const auto [corr,cov]=
	loadCorrCov(chann);
      
      return corr;
    };
    
  const djvec_t P5P5_ll=
    loadCorr("mix_fixed_l_l_TM_P5P5");
  
  /// Estimates the lattice spacing
  const djack_t aInv=
    estimateA(aMlist.at('C'),P5P5_ll,T,baseOut+"/"+ensName+"/pion.xmg");
  cout<<"a^-1: "<<aInv.ave_err()<<" GeV^-1"<<endl;
  
#define LOAD(TAG) \
  const auto [c ## TAG ## _ls,cov ## TAG ## _ls]=\
    loadCorrCov("mix_fixed_l_s1_TM_" #TAG)

  LOAD(V0V0);
  LOAD(VKVK);
  LOAD(A0A0);
  LOAD(AKAK);
  LOAD(P5P5);
  
#undef LOAD
  
  const double EMinIntInGeV=0.27;
  const double EMinInt=EMinIntInGeV/aInv.ave();
  const double EMaxInt=6;
  
  const double sigma=0.03;
  
  const size_t tMin=1,tMax=36;
  const Basis<double> basis(T,tMin,tMax,+1);
  const size_t nT=basis.nT;
  
  djvec_t corr(nT);
  
  MatrixX<double> corrCov(nT,nT);
  for(size_t iT=0;iT<nT;iT++)
    {
      corr[iT]=cVKVK_ls[iT+tMin];
      for(size_t iS=0;iS<nT;iS++)
	corrCov(iT,iS)=covVKVK_ls[(iT+tMin)+(T/2+1)*(iS+tMin)];
    }
  
  /// Preconditioner of the problem
  vector<double> preco(nT);
  for(size_t iT=0;iT<basis.nT;iT++)
    preco[iT]=corr[iT].ave();
  
  const size_t RE_exp=2;
  
  /// Normalization of the correlation function, when assuming that R(E)=E^N
  djvec_t normOfT(nT);
  for(size_t iT=0;iT<nT;iT++)
    {
      /// Correlator assuming that R(E)=E^N
      const double cAss=
	gslIntegrateFromTo([iT,&basis](const double& E)
	{
	  return pow(E,RE_exp)*basis(iT,E);
	},0.5/aInv.ave(),EMaxInt);
      
      cout<<iT<<" "<<corr[iT].ave_err()<<" "<<cAss<<endl;
      
      normOfT[iT]=corr[iT]/cAss;
    }
  normOfT.ave_err().write(baseOut+"/"+ensName+"/norm.xmg");
  
  /// Time at which to normalize the correlator: we choose 2 GeV^-1
  const size_t tNorm=round(2*aInv.ave());
  cout<<"tNorm: "<<tNorm<<endl;
  
  /// Ansatz of the spectral density
  const auto specAns=
    [norm=normOfT[tNorm].ave()](const double& E)
    {
      return norm*pow(E,RE_exp);
    };
  
  const auto targetFunction=
    [sigma](const double& x)
    {
      return Kt(sigma,x);
    };
  
  const TrambaccoFunctional tf=
    getTrambaccoFunctional(basis,
			   specAns,
			   targetFunction,
			   EMinInt,EMaxInt,corrCov,preco);
  
  /// Track the results
  djack_t R;
  VectorXd g;
  
  const auto reconstruct=
    [&](const double& statOverSyst=1)
    {
      g=tf.getG(preco,statOverSyst);
      R=0;
      for(size_t iT=0;iT<nT;iT++)
	R+=g[iT]*corr[iT];
    };
  
  // Performs the stability check
  grace_file_t stabilityPlot(baseOut+"/"+ensName+"/stab.xmg");
  vector<pair<double,djack_t>> stab;
  for(double statOverSyst=512;statOverSyst>=1e-10;statOverSyst/=2)
    {
      reconstruct(statOverSyst);
      stab.emplace_back(log(statOverSyst),R);
    }
  
  // Final reconstruction
  const double lambdaExp=-2.5;
  reconstruct(exp(lambdaExp));
  
  for(const auto& [x,y] : stab)
    stabilityPlot.write_ave_err(x,y.ave_err());
  stabilityPlot.write_constant_band(lambdaExp-0.1,lambdaExp+0.1,R);
  
  /// Print coefficients
  grace_file_t gPlot(baseOut+"/"+ensName+"/g.xmg");
  for(size_t iT=0;iT<nT;iT++)
    gPlot.write_xy(iT,g[iT]);
  
  const auto interpReco=
    [g,&aInv,&nT,&basis](const double& EInGeV)
    {
      const double E=EInGeV/aInv.ave();
      double s=0;
      for(size_t iT=0;iT<nT;iT++)
	s+=g[iT]*basis(iT,E);
      
      return s;
    };
  
  // Plots the reconstruction of the target function
  grace_file_t recoPlot(baseOut+"/"+ensName+"/reco.xmg");
  recoPlot.write_line(interpReco,EMinInt*aInv.ave(),EMaxInt*aInv.ave());
  
  recoPlot.write_line([&](const double& EInGeV)
  {
    const double E=EInGeV/aInv.ave();
    return targetFunction(E);
  },EMinInt*aInv.ave(),EMaxInt*aInv.ave());
  
  return 0;
}
