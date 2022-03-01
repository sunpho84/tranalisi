#include <tranalisi.hpp>

int T=128;

djvec_t loadOPPO(const string& tag,const int parity=+1)
{
  return ((read_djvec("data/mes_contr_"+tag,T,16)+read_djvec("data/mes_contr_"+tag,T,19))/2).symmetrize(parity);
}

djvec_t loadSAME(const string& tag,const int parity=+1)
{
  return read_djvec("data/mes_contr_"+tag,T,17).symmetrize(parity);
}

int main()
{
  const double amq=0.21000;
  const int tMinP5P5=40;
  const int tMaxP5P5=50;
  
  set_njacks(50);
  
  const djvec_t corrP5P5[2]={loadOPPO("P5P5"),loadSAME("P5P5")};
  
  // P5P5_oppo.ave_err().write("plots/P5P5_OPPO.xmg");
  // P5P5_same.ave_err().write("plots/P5P5_SAME.xmg");
  
  const djvec_t corrA0P5[2]={-loadOPPO("A0P5",-1),-loadSAME("A0P5",-1)};
  // const djvec_t P5A0_oppo=loadOPPO("P5A0",-1);
  
  // A0P5_oppo.ave_err().write("plots/A0P5_OPPO.xmg");
  // A0P5_same.ave_err().write("plots/A0P5_SAME.xmg");
  // P5A0_oppo.ave_err().write("plots/P5A0_OPPO.xmg");
  
  djack_t mP[2],ZA0[2],ZP5[2];
  for(int iMes=0;iMes<2;iMes++)
    two_pts_SL_fit(ZP5[iMes],ZA0[iMes],mP[iMes],corrA0P5[iMes],corrP5P5[iMes],T/2,41,57,combine("plots/A0P5FitMes%zu.xmg",iMes),-1,+1);
  const djvec_t derP5A0TM=-symmetric_derivative(corrA0P5[0]);
  derP5A0TM.ave_err().write("plots/derP5A0_regoTM.xmg");
  const djvec_t ZvSilvCorr=2*amq*corrP5P5[0]/derP5A0TM;
  const djack_t ZvSilv=constant_fit(ZvSilvCorr,tMinP5P5,tMaxP5P5,"plots/ZvSilv.xmg");
  cout<<"ZvSilv: "<<ZvSilv.ave_err()<<endl;
  
  const djack_t ZaSilvCorrectingFactor=ZP5[0]/ZP5[1]*mP[1]*sinh(mP[1])/(mP[0]*sinh(mP[0]));
  cout<<"ZaSilvCorrectingFactor: "<<ZaSilvCorrectingFactor.ave_err()<<endl;
  
  const djvec_t derP5A0OS=-symmetric_derivative(corrA0P5[1]);
  // derP5A0OS.ave_err().write("plots/derP5A0_regoOS.xmg");
  const djvec_t ZaSilvUnCorr=2*amq*corrP5P5[1]/derP5A0OS;
  const djvec_t ZaSilvCorr=ZaSilvUnCorr*ZaSilvCorrectingFactor;
  // const djack_t ZaSilvUn=constant_fit(ZaSilvUnCorr,tMinP5P5[0],tMaxP5P5[0],"plots/ZaSilvUncorr.xmg");
  const djack_t ZaSilv=constant_fit(ZaSilvCorr,18,30,"plots/ZaSilv.xmg");
  cout<<"ZaSilv: "<<ZaSilv.ave_err()<<endl;
  // Zsilv[1]=ZaSilv;
  
  return 0;
}
