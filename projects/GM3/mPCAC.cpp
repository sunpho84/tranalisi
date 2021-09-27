#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <fit.hpp>
#include <meas_vec.hpp>

#include <data.hpp>
#include <physLine.hpp>

/// Compute the PCAC mass
void computeMPCAC()
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
  
  const double fPiPhys=
    0.1304;
  
  const double mPiPhys=
    0.13498;
  
  const double r2Phys=
    sqr(mPiPhys/0.1304);
  
  const djack_t aFromFPi=
    aFpi/fPiPhys;
  
  console<<" r2: "<<r2.ave_err()<<endl;
  console<<" r2Phys: "<<r2Phys<<endl;
  console<<" a from aFpi: "<<aFromFPi.ave_err()<<endl;
  
  *aFromPhysLine=
    jackCall(aFromMPiFPi,aMpi,aFpi);
  console<<" a from phys line analysis: "<<aFromPhysLine->ave_err()<<endl;
  
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

