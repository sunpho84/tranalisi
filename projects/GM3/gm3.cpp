#define EXTERN
#define INIT_EXTERN_TO(ARGS...) ARGS

#include <aLaLuscherRepresentation.hpp>
#include <amu.hpp>
#include <cachedData.hpp>
#include <data.hpp>
#include <fit.hpp>
#include <kernels.hpp>
#include <params.hpp>
#include <rawData.hpp>
#include <renoConstants.hpp>

using namespace std;

void testMf()
{
  const djvec_t cP5P5=
    getAveForRego(0,nSources,idP5P5,REGO_TM);
  
  djack_t Z2,aMpi;
  two_pts_fit(Z2,aMpi,cP5P5,TH,tMinP5P5[0],tMaxP5P5[0],"plots/P5P5fitForMpi.xmg");
  
  const djack_t aFpi=
    2*amq*sqrt(Z2)/(aMpi*sinh(aMpi));
  
  console<<"aMpi: "<<aMpi<<endl;
  console<<"aFpi: "<<aFpi<<endl;
  
  const djack_t r=
    aMpi/aFpi;
  
  const djack_t r2=
    r*r;
  
  const double r2Phys=
    sqr(134.98/130.4);
  
  console<<" r2: "<<r2.ave_err()<<endl;
  console<<" r2Phys: "<<r2Phys<<endl;
  
  const djack_t r2frR2Phys=
    r2/r2Phys;
  
  console<<" r2/r2Phys: "<<r2frR2Phys.ave_err()<<endl;
  
  const djvec_t cV0P5=
    (getAve(0,nSources,idV0P5,0)-getAve(0,nSources,idV0P5,2))/2;
  cV0P5.ave_err().write("plots/corr_V0P5_usable.xmg");
  
  const djvec_t cV0P5der=
    symmetric_derivative(cV0P5);
  
  cV0P5der.ave_err().write("plots/corr_V0P5_der.xmg");
  
  const djvec_t mPCACEff=
    cV0P5der/(2*cP5P5);
  
  const djack_t mPCAC=
    constant_fit(mPCACEff,14,TH-1,"plots/eff_mPCAC.xmg");
  
  console<<"mPCAC: "<<mPCAC.ave_err()<<endl;
}

template <typename T>
T dualPartFun(const T& mRho,
	      const T& eDual,
	      const T& rDual,
	      const double& t)
{
  const T mRhoPlusEDual=
    mRho+eDual;
	      
  return
    5/(18*M_PI*M_PI)*
    rDual/(t*t*t)*
    exp(-mRhoPlusEDual*t)*
    (1+mRhoPlusEDual*t+sqr(mRhoPlusEDual)*t*t/2);
}

void fitVKVK()
{
  const double qMax=
    4.3;
  
  const int nIntervals=
    1000;
  
  HashedTanPhiAndDerivFunction hashedPhiAndDerivCalculator("lookupTables/phi.hash","lookupTables/phiDeriv.hash",qMax,nIntervals);
  hashedPhiAndDerivCalculator.plot("plots/phi.xmg","plots/phiDeriv.xmg");
  
  const djvec_t cP5P5=
    getAveForRego(0,nSources,idP5P5,REGO_TM);
  
  /// Charge factor of the correlator
  const double chargeFactor=
    5.0/9;
  
  djvec_t cVKVK=
    getAveForRego(0,nSources,idVKVK,REGO_TM)*sqr(Z[regoZId[REGO_TM]])*chargeFactor;
  cVKVK.ave_err().write("plots/corr_for_Rep.xmg");
  cVKVK[0]=cVKVK[1];
  
  jack_fit_t fitter;
  
  djvec_t pars(4);
  pars[0].fill_gauss(1.15,0.05,235235);
  pars[1].fill_gauss(0.17,0.05,7342);
  pars[2].fill_gauss(0.32,0.05,23423);
  pars[3].fill_gauss(28.1,0.1,32235);
  
  const size_t iRDual=
    fitter.add_fit_par_limits(pars[0],"RDual",pars[0].ave(),pars[0].err(), 0.4,1.5);
  const size_t iEThr=
    fitter.add_fit_par_limits(pars[1],"EThr",pars[1].ave(),pars[1].err(), 0.01,1.0);
  const size_t iMRho=
    fitter.add_fit_par_limits(pars[2],"MRho",pars[2].ave(),pars[2].err(), 0.01,1.0);
  const size_t iG2=
    fitter.add_fit_par_limits(pars[3],"g2",pars[3].ave(),pars[3].err(), 25.0,35.0);
  
  const djack_t aMPi=
    constant_fit(effective_mass(cP5P5),tMinP5P5[0],tMaxP5P5[0]);
  
  const int n=6;
  ALaLuscherRepresentationCached LuschRep(hashedPhiAndDerivCalculator,n);
  
  // {
  //   const double amRho=
  //     pars[iXMRho][0]*aMPi[0];
    
  //   const double& g2=
  //     pars[iG2][0];
    
  //   for(const auto& c : LuschRep(aMPi[0],L,amRho,g2).coeffs)
  //     cout<<c.weight<<" "<<c.energy<<endl;
  // }
  
  // CRASH("");
  
  auto fullRep=
    [&LuschRep,
     iRDual,
     iEThr,
     iMRho,
     iG2](const vector<double>& p,
	  const double& aMPi,
	  const double& t)
    {
      const double& rDual=
	p[iRDual];
      const double eThr=
	p[iEThr];
      const double amRho=
	p[iMRho];
      const double& g2=
	p[iG2];
      
      const double dualPart=
	dualPartFun(amRho,eThr,rDual,t);
      
      const auto LuschRepPars=
	LuschRep(aMPi,L,amRho,g2);
      
      static bool b=true;
      
      if(b)
	for(auto& p : LuschRepPars.coeffs)
	  cout<<p.weight<<" "<<p.energy<<endl;
      
      b=false;
      
      const double& LuscherPart=
	LuschRepPars(t);
      
      return
	std::pair<double,double>{dualPart,LuscherPart};
    };
  
  const size_t tMinFit=
    5;
  
  size_t tFit=
    tMinFit;
  
  while(tFit<THp1 and cVKVK[tFit].err()<0.1*fabs(cVKVK[tFit].ave()))
    {
      fitter.add_point(cVKVK[tFit],
		       [&aMPiAll=aMPi,
			&fullRep,
			tFit](const vector<double>& p,
			   const size_t ijack)->double
		       {
			 const double aMPi=
			   aMPiAll[ijack];
			 
			 const auto [dualPart,LuscherPart]=
			   fullRep(p,aMPi,tFit);
			 
			 return
			   dualPart+LuscherPart;
		       });
      tFit++;
    };
  const size_t tMaxFit=
    tFit;
  
  const djvec_t oldPars=
    pars;
  
  fitter.fit();
  
  vector<vector<double>> parsPerJack(njacks+1,vector<double>(4));
  for(int i=0;i<4;i++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      parsPerJack[ijack][i]=pars[i][ijack];
  
  auto jackUnbindRep=
    [&parsPerJack,
     &aMPi,
     &fullRep](const double& t) -> tuple<djack_t,djack_t>
    {
      djack_t dual,Luscher;
      
      for(size_t ijack=0;ijack<=njacks;ijack++)
	tie(dual[ijack],Luscher[ijack])=
	  fullRep(parsPerJack[ijack],aMPi[ijack],t);
      
      return
	{dual,Luscher};
    };
  
  djvec_t cVKVKDual(THp1);
  djvec_t cVKVKLuscher(THp1);
  djvec_t cVKVKFull(THp1);
  for(size_t t=0;t<=TH;t++)
    tie(cVKVKDual[t],cVKVKLuscher[t])=
      jackUnbindRep(t);
  
  cVKVKFull=
    cVKVKLuscher+cVKVKDual;
  
  auto jackFullRep=
    [&jackUnbindRep](const double& t) -> djack_t
    {
      const auto [dual,Luscher]=
	jackUnbindRep(t);
      
      return
	dual+Luscher;
    };
  
  cVKVKDual.ave_err().write("plots/dualRep.xmg");
  cVKVKLuscher.ave_err().write("plots/LuscherRep.xmg");
  cVKVKFull.ave_err().write("plots/fullRep.xmg");
  
  cout<<oldPars.ave_err()<<endl;
  cout<<"---"<<endl;
  cout<<pars.ave_err()<<endl;
  
  auto effMassFun=
    [&jackFullRep](const double& tPlot)
    {
      const double dTPlot=
	1;
      
      const djack_t y=
	jackFullRep(tPlot);
      
      const djack_t yPlusDt=
	jackFullRep(tPlot+dTPlot);
      
      const djack_t e=
	effective_mass(y,yPlusDt,tPlot,TH,1,0,dTPlot);
      
      return
	e;
    };
  
  const djack_t& rDual=
	pars[iRDual];
  const djack_t eThr=
	pars[iEThr];
  const djack_t amRho=
	pars[iMRho];
  const djack_t g=
    sqrt(pars[iG2]);
  
  cout<<"rDual: "<<rDual.ave_err()<<endl;
  cout<<"eThr: "<<eThr.ave_err()<<endl;
  cout<<"amRho: "<<amRho.ave_err()<<endl;
  cout<<"g: "<<g.ave_err()<<endl;
  
  grace_file_t plotFit("plots/aLaLuscherFit.xmg");
  plotFit.write_vec_ave_err(effective_mass(cVKVK).ave_err());
  plotFit.set_no_line();
  plotFit.write_polygon(effMassFun,tMinFit,tMaxFit);
  plotFit.write_polygon(effMassFun,tMaxFit,TH-1,grace::GREEN4);
}

int main(int narg,char **arg)
{
  init(narg,arg);
  
  loadData(narg,arg);
  
  Z=determineRenoConst();
  
  if(canUseRawData)
    {
      analyzeRawData();
      convertForSilvano();
    }
  
  computeAmu();
  
  testMf();
  
  fitVKVK();
  
  close();
  
  return 0;
}
