#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_CORRECTIONS
 #include <MOM2/corrections.hpp>

void set_pr_bil_a2(gaz::type_t iaz,gf::type_t gf_g)
{
  // static bool init=false;
  // static gaz::type_t old_iaz;
  // static gf::type_t old_gf_g;
  
  // if(not init or old_iaz!=iaz or old_gf_g!=gf_g)
  //   {
  //     init=true;
  //     old_iaz=iaz;
  //     old_gf_g=gf_g;
      
      double *cS=pr_bil_a2_c[iS];
      double *cV=pr_bil_a2_c[iV];
      double *cA=pr_bil_a2_c[iA];
      double *cP=pr_bil_a2_c[iP];
      double *cT=pr_bil_a2_c[iT];
      
      double lambda=gf::lambda[gf_g];
      
      cS[0]=EpsS2[iaz][1]-lambda*2.27358943;
      cS[1]=-1.0/4.0+lambda*3.0/4.0;
      cS[2]=13.0/24.0+C2[iaz]/2.0-lambda/8.0;
      
      cP[0]=EpsP2[iaz][1]-lambda*0.83810121;
      cP[1]=-1.0/4.0+lambda/4.0;
      cP[2]=13.0/24.0+C2[iaz]/2.0-lambda/8.0;
      
      //including Lorentz index normalization (4,4,6)
      
      cV[0]=(EpsV2[iaz][4]-lambda*0.8110353+(EpsV2[iaz][1]+lambda/8.0+EpsV2[iaz][7]+lambda*0.2436436)/4.0)/4.0;
      cV[1]=(11.0/240.0-c1[iaz]/2.0-C2[iaz]/60.0+lambda/8.0+(-53.0/120.0+11.0*C2[iaz]/10.0-149.0/120.0-c1[iaz]-C2[iaz]/30.0+lambda/4.0)/4.0)/4.0;
      cV[2]=(3.0/80.0+C2[iaz]/10.0+lambda*5.0/48.0+(-101.0/60.0+11.0*C2[iaz]/15.0+lambda/3.0-1.0/60.0+2*C2[iaz]/5.0+lambda/12.0-3.0/40.0-C2[iaz]/5.0-lambda*5.0/24.0)/4.0)/4.0;
      
      cA[0]=(EpsA2[iaz][4]-lambda*1.7465235+(EpsA2[iaz][1]+lambda/8.0+EpsA2[iaz][7]+lambda*1.1146200)/4.0)/4.0;
      cA[1]=(-109.0/240-c1[iaz]/2.0-C2[iaz]/60.0+lambda*5.0/8.0+(-53.0/120.0+11.0*C2[iaz]/10.0+91.0/120.0-c1[iaz]-C2[iaz]/30.0-lambda*3.0/4.0)/4.0)/4.0;
      cA[2]=(3.0/80.0+C2[iaz]/10.0+lambda*5.0/48.0+(-101.0/60.0+11.0*C2[iaz]/15.0+lambda/3.0-1.0/60.0+2.0*C2[iaz]/5.0+lambda/12.0-3.0/40.0-C2[iaz]/5.0-lambda*5.0/24.0)/4.0)/4.0;
      
      cT[0]=((EpsT2prime[iaz][1]+lambda/4.0-EpsT2prime[iaz][4]+lambda*1.12097643)/4.0+EpsT2prime[iaz][7]-lambda*1.2194576+(-7.0/3.0-4.0*C2[iaz]/3.0-lambda/2.0)/12.0)/6.0;
      cT[1]=(-5.0/12.0-c1[iaz]+C2[iaz]/3.0+lambda/4.0)/6.0;
      cT[2]=((-26.0+10.0*C2[iaz]+9.0*lambda)/36.0)/6.0;
    // }
}
