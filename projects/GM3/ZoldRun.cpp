#include <tranalisi.hpp>

int T;

djvec_t loadOPPO(const string& tag,const int parity=+1)
{
  return ((read_djvec("data/mes_contr_"+tag,T,20)+read_djvec("data/mes_contr_"+tag,T,23))/2).symmetrize(parity);
}

djvec_t loadSAME(const string& tag,const int parity=+1)
{
  return read_djvec("data/mes_contr_"+tag,T,21).symmetrize(parity);
}

int main()
{
  input_file_t input("input");
  T=input.read<int>("T");
  const double amq=input.read<double>("amq");
  const int tMinP5P5=input.read<double>("tMin");
  const int tMaxP5P5=input.read<double>("tMax");
  
  set_njacks(input.read<int>("nJacks"));
  
  const djvec_t corrP5P5[2]={loadOPPO("P5P5"),loadSAME("P5P5")};
  
  const djvec_t corrA0P5[2]={-loadOPPO("A0P5",-1),-loadSAME("A0P5",-1)};
  
  djack_t mP[2],ZA0[2],ZP5[2];
  for(int iMes=0;iMes<2;iMes++)
    two_pts_SL_fit(ZP5[iMes],ZA0[iMes],mP[iMes],corrA0P5[iMes],corrP5P5[iMes],T/2,tMinP5P5,tMaxP5P5,combine("plots/A0P5FitMes%d.xmg",iMes),-1,+1);
  const djvec_t derP5A0TM=-symmetric_derivative(corrA0P5[0]);
  derP5A0TM.ave_err().write("plots/derP5A0_regoTM.xmg");
  const djvec_t ZvSilvCorr=2*amq*corrP5P5[0]/derP5A0TM;
  const djack_t ZvSilv=constant_fit(ZvSilvCorr,tMinP5P5,tMaxP5P5,"plots/ZvSilv.xmg");
  cout<<"ZvSilv: "<<ZvSilv.ave_err()<<endl;
  
  const djack_t ZaSilvCorrectingFactor=ZP5[0]/ZP5[1]*mP[1]*sinh(mP[1])/(mP[0]*sinh(mP[0]));
  // cout<<"ZaSilvCorrectingFactor: "<<ZaSilvCorrectingFactor.ave_err()<<endl;
  
  const djvec_t derP5A0OS=-symmetric_derivative(corrA0P5[1]);
  const djvec_t ZaSilvUnCorr=2*amq*corrP5P5[1]/derP5A0OS;
  const djvec_t ZaSilvCorr=ZaSilvUnCorr*ZaSilvCorrectingFactor;
  const djack_t ZaSilv=constant_fit(ZaSilvCorr,tMinP5P5,tMaxP5P5,"plots/ZaSilv.xmg");
  cout<<"ZaSilv: "<<ZaSilv.ave_err()<<endl;
  
  return 0;
}
