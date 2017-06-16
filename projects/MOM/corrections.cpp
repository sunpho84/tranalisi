#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_CORRECTIONS
 #include <corrections.hpp>

void set_pr_bil_a2(gaz_t iaz)
{
  double *cS=pr_bil_a2_c[iZS];
  double *cV=pr_bil_a2_c[iZV];
  double *cA=pr_bil_a2_c[iZA];
  double *cP=pr_bil_a2_c[iZP];
  double *cT=pr_bil_a2_c[iZT];
  
  cS[0]=-EpsS2[iaz][1];
  cS[1]=-1.0/4.0;
  cS[2]=13.0/24.0+C2[iaz]/2.0;
  
  cP[0]=-EpsP2[iaz][1];
  cP[1]=cS[1];
  cP[2]=cP[1];
  
  cV[0]=EpsV2[iaz][4]+(EpsV2[iaz][1]+EpsV2[iaz][7])/4.0;
  cV[1]=(-3.0-6.0*c1[iaz]+2*C2[iaz])/8.0;
  cV[2]=13.0/32.0+C2[iaz]/3;
  
  cA[0]=EpsA2[iaz][4]+(EpsA2[iaz][1]+EpsA2[iaz][7])/4.0;
  cA[1]=cV[1];
  cA[2]=cV[2];
  
  cT[0]=(EpsT2prime[iaz][1]-EpsT2prime[iaz][4])/4.0+EpsT2prime[iaz][7]+(-7.0/3.0-4.0*C2[iaz]/3.0)/12.0;
  cT[1]=-5.0/12.0-c1[iaz]+C2[iaz]/3.0;
  cT[2]=(-13.0+5.0*C2[iaz])/18.0;
}
