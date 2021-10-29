#include <tranalisi.hpp>

const size_t T=48,TH=T/2;

int main()
{
  set_njacks(15);
  
  const djvec_t P5P5_TM=
    (read_djvec("data/mes_contr_M0_R0_0_M0_R0_0_P5P5",T,0).symmetrized()+
     read_djvec("data/mes_contr_M0_R1_0_M0_R1_0_P5P5",T,0).symmetrized())/2;
  
  const djvec_t A0P5_TM=
    (read_djvec("data/mes_contr_M0_R0_0_M0_R0_0_A0P5",T,0).symmetrized(-1)+
     read_djvec("data/mes_contr_M0_R1_0_M0_R1_0_A0P5",T,0).symmetrized(-1))/2;
  
  // const djvec_t P5A0_TM=
  //   (read_djvec("data/mes_contr_M0_R0_0_M0_R0_0_P5A0",T,0).symmetrized(-1)+
  //    read_djvec("data/mes_contr_M0_R1_0_M0_R1_0_P5A0",T,0).symmetrized(-1))/2;
  
  P5P5_TM.ave_err().write("plots/P5P5_TM.xmg");
  // P5A0_TM.ave_err().write("plots/P5A0_TM.xmg");
  A0P5_TM.ave_err().write("plots/A0P5_TM.xmg");
  
  const djvec_t P5P5_OS=
    read_djvec("data/mes_contr_M0_R0_0_M0_R1_0_P5P5",T,0).symmetrized();
  
  const djvec_t S0P5_OS=
    read_djvec("data/mes_contr_M0_R0_0_M0_R1_0_S0P5",T,0).symmetrized();
  
  const djvec_t A0P5_OS=
    read_djvec("data/mes_contr_M0_R0_0_M0_R1_0_A0P5",T,0).symmetrized(-1);
  
  // const djvec_t P5A0_OS=
  //   read_djvec("data/mes_contr_M0_R1_0_M0_R1_0_P5A0",T,0).symmetrized(-1);
  
  effective_mass(P5P5_TM,TH,-1).ave_err().write("plots/P5P5_TM.xmg");
  // effective_mass(P5A0_TM,TH,-1).ave_err().write("plots/P5A0_TM.xmg");
  effective_mass(A0P5_TM,TH,-1).ave_err().write("plots/A0P5_TM.xmg");
  
  effective_mass(P5P5_OS,TH,+1).ave_err().write("plots/P5P5_OS.xmg");
  // effective_mass(P5A0_OS,TH,-1).ave_err().write("plots/P5A0_OS.xmg");
  effective_mass(A0P5_OS,TH,-1).ave_err().write("plots/A0P5_OS.xmg");
  // effective_mass(
		 S0P5_OS// ,TH,+1)
    .ave_err().write("plots/S0P5_OS.xmg");
  
  const double amq=0.0100;
  
  const djvec_t derA0P5TM=symmetric_derivative(A0P5_TM);
  derA0P5TM.ave_err().write("plots/derP5A0_regoTM.xmg");
  const djvec_t ZvSilvCorr=2*amq*P5P5_TM/derA0P5TM;
  const djack_t ZvSilv=constant_fit(ZvSilvCorr,12,23,"plots/ZvSilv.xmg");
  cout<<"ZvSilv: "<<ZvSilv.ave_err()<<endl;
  
  djack_t mP_OS,mP_TM;
  djack_t ZA0_OS,ZA0_TM;
  djack_t ZP5_OS,ZP5_TM;
  
  two_pts_SL_fit(ZP5_TM,ZA0_TM,mP_TM,A0P5_TM,P5P5_TM,TH,12,23,"plots/A0P5FitMesTM.xmg",-1,+1);
  two_pts_SL_fit(ZP5_OS,ZA0_OS,mP_OS,A0P5_OS,P5P5_OS,TH,12,23,"plots/A0P5FitMesOS.xmg",-1,+1);
  const djack_t ZaSilvCorrectingFactor=ZP5_TM/ZP5_OS*mP_OS*sinh(mP_OS)/(mP_TM*sinh(mP_TM));
  cout<<"ZaSilvCorrectingFactor: "<<ZaSilvCorrectingFactor.ave_err()<<endl;
  
  const djvec_t derA0P5OS=symmetric_derivative(A0P5_OS);
  derA0P5OS.ave_err().write("plots/derP5A0_regoOS.xmg");
  const djvec_t ZaSilvUnCorr=2*amq*P5P5_OS/derA0P5OS;
  const djvec_t ZaSilvCorr=ZaSilvUnCorr*ZaSilvCorrectingFactor;
  const djack_t ZaSilvUn=constant_fit(ZaSilvUnCorr,12,23,"plots/ZaSilvUncorr.xmg");
  const djack_t ZaSilv=constant_fit(ZaSilvCorr,12,23,"plots/ZaSilv.xmg");
  cout<<"ZaSilv: "<<ZaSilv.ave_err()<<endl;
  
  const djack_t fPfromP_TM=2.0*ZP5_TM*amq/sqr(mP_TM);
  const djack_t fPfromA_TM=ZA0_TM/mP_TM;
  const djack_t Zv=fPfromP_TM/fPfromA_TM;
  cout<<"Zv: "<<Zv.ave_err()<<endl;

  // const djack_t fPfromP_OS=2.0*ZP5_OS*amq/sqr(mP_OS);
  const djack_t fPfromA_OS=ZA0_OS/mP_OS;
  const djack_t Za=fPfromP_TM/fPfromA_OS;
  cout<<"Za: "<<Za.ave_err()<<endl;
  
  return 0;
}
