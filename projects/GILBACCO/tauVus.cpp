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
      
      corr.ave_err().write(baseOut+ensName+"/"+chann+".xmg");
      
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
  
  const auto [cV0V0_ls,corrV0V0_ls]=
    loadCorrCov("mix_fixed_l_s1_TM_V0V0");
  const auto [cVKVK_ls,corrVKVK_ls]=
    loadCorrCov("mix_fixed_l_s1_TM_VKVK");
  
  const double EMinIntInGeV=0.09;
  const double EMinInt=EMinIntInGeV/aInv.ave();
  const double EMaxInt=4;
  
  const double sigma=0.3;
  
  /// Preconditioner of the problem
  vector<Real> preco(nT);
  for(size_t iT=0;iT<nT;iT++)
    preco[iT]=corr[iT].ave();
  
  const TrambaccoFunctional tf=
    getTrambaccoFunctional(basis,
			   specAns,
			   [sigma](const double& x)
			   {
			     return Kt(sigma,x);
			   }
			   ,EMinInt,EMaxInt,cVKVK_ls,preco);
  
  grace_file_t temp("/tmp/t.xmg");
  temp.write_line(bind(thetaSigma,sigma,std::placeholders::_1),-10,10);
  
  return 0;
}
