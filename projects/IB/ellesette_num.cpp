#include <string>
#include <tranalisi.hpp>

const size_t L=24;
const size_t T=48,TH=T/2;

djvec_t read(const string& suff)
{
  return read_djvec("jacks/M0_R0_M0_R0"+suff,T,0).symmetrized();
}

int main()
{
  set_njacks(15);
  
  const djvec_t isoQcd=read("");
  
  const int ndelta=4,ncontr=3;
  
  djvec_t contr[ndelta][ncontr];
  const string deltaTag[4]={"","0","1","2"};
  const string contrTag[3]={"uu","ud","dd"};
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      contr[0][icontr]=isoQcd;
      for(int idelta=1;idelta<ndelta;idelta++)
	contr[idelta][icontr]=read("_"+contrTag[icontr]+"_"+deltaTag[idelta]);
    }
  
  const int npions=2;
  djvec_t corr[ndelta][npions];
  for(int idelta=0;idelta<ndelta;idelta++)
    {
      corr[idelta][0]=(contr[idelta][0]+contr[idelta][2])/2.0;
      corr[idelta][1]=contr[idelta][1];
    }
  
  const int nder=3;
  const double delta[ndelta]={0,0.0004,0.0008,0.002};
  djvec_t corrSecond[nder][npions];
  for(int ider=0;ider<nder;ider++)
    for(int ipion=0;ipion<npions;ipion++)
      corrSecond[ider][ipion]=(corr[ider+1][ipion]-isoQcd)/(2*sqr(delta[ider+1]));
  
  for(int ider=0;ider<nder;ider++)
    for(int ipion=0;ipion<npions;ipion++)
      corrSecond[ider][ipion].ave_err().write("plots/corrSecondPi"+to_string(ipion)+"Der"+to_string(ider)+".xmg");
  
  const djvec_t eff_isoQcd=effective_mass(isoQcd);
  const djack_t m=constant_fit(eff_isoQcd,10,15,"plots/eff_mass_isoQcd.xmg");
  cout<<"M: "<<m.ave_err()<<endl;
  
  for(int ider=0;ider<nder;ider++)
    {
      const djvec_t corrSecondDiff=corrSecond[ider][1]-corrSecond[ider][0];
      corrSecondDiff.ave_err().write("plots/corrSecondDiff"+to_string(ider)+".xmg");
      
      const djvec_t corrSecondDiffRat=corrSecondDiff/isoQcd;
      corrSecondDiffRat.ave_err().write("plots/corrSecondDiffRat"+to_string(ider)+".xmg");
      const djvec_t eff_slope_corrSecondDiff=effective_slope(corrSecondDiffRat,eff_isoQcd,TH);
      eff_slope_corrSecondDiff.ave_err().write("plots/corrSecondDiffSlope"+to_string(ider)+".xmg");
    }
  
  const djvec_t U_U=read("_uu_0");
  const djvec_t D_D=read("_dd_0");
  const djvec_t U_D=read("_ud_0");
  
  const double dm=0.0004;
  const djvec_t test=(U_D-(U_U+D_D)/2)/(2*dm*dm*isoQcd);
  const djvec_t eff_test=effective_slope(test,eff_isoQcd,TH);
  eff_test.ave_err().write("plots/test.xmg");
  
  return 0;
}
