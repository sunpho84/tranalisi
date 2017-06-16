#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_CORRECTIONS
 #include <corrections.hpp>

void set_pr_bil_a2(gaz_t iaz)
{
  double *cS=pr_bil_a2[iZS];
  double *cV=pr_bil_a2[iZV];
  double *cA=pr_bil_a2[iZA];
  double *cP=pr_bil_a2[iZP];
  double *cT=pr_bil_a2[iZT];
  
  cS[0]=-EpsS2[iaz][1];
  cS[1]=-1.0/4.0;
  cS[2]=13.0/24.0+C2[iaz]/2.0;
  
  cP[0]=-EpsP2[iaz][1];
  cP[1]=cS[1];
  cP[2]=cP[1];
  
  cV[0]=EpsV2[iaz][4]+(EpsV2[iaz][1]+EpsV2[iaz][7])/4.0;
  cV[1]=(-3.0-6.0*c1[iaz]+2*C2[iaz])/8.0;
  
}
