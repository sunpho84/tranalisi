#define EXTERN
#define INIT_EXTERN_TO(ARGS...) ARGS

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

int main(int narg,char **arg)
{
  init(narg,arg);
  
  loadData(narg,arg);
  
  const djvec_t Z=
    determineRenoConst();
  
  if(canUseRawData)
    {
      analyzeRawData();
      convertForSilvano();
    }
  
  computeAmu(Z);
  
  testMf();
  
  close();
  
  return 0;
}
