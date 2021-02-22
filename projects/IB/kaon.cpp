#include <tranalisi.hpp>

const size_t T=48,TH=T/2;
const double ainv=2.2;

const double eu=+2.0/3;
const double ed=-1.0/3;
const double es=-1.0/3;
const double adml=0.00032378;

const double Z_QED=(6.0*log(mu_MS/ainv)-22.5954)/(16.0*sqr(M_PI));

int main()
{
  set_njacks(22);
  
  const djvec_t M0_M0=read_djvec("mes_contr_P5P5_M0_R0_0_M0_R0_0",T).symmetrized();
  const djvec_t M0_MM=read_djvec("mes_contr_P5P5_M0_R0_0_M0_R0_MASS",T).symmetrized();
  const djvec_t M0_MQ=read_djvec("mes_contr_P5P5_M0_R0_0_M0_R0_QED",T).symmetrized();
  const djvec_t M0_N0=read_djvec("mes_contr_P5P5_M0_R0_0_N0_R0_0",T).symmetrized();
  const djvec_t M0_NM=read_djvec("mes_contr_P5P5_M0_R0_0_N0_R0_MASS",T).symmetrized();
  const djvec_t M0_NQ=read_djvec("mes_contr_P5P5_M0_R0_0_N0_R0_QED",T).symmetrized();
  const djvec_t MF_MF=read_djvec("mes_contr_P5P5_M0_R0_F_M0_R0_F",T).symmetrized();
  const djvec_t MF_NF=read_djvec("mes_contr_P5P5_M0_R0_F_N0_R0_F",T).symmetrized();
  const djvec_t N0_M0=read_djvec("mes_contr_P5P5_N0_R0_0_M0_R0_0",T).symmetrized();
  const djvec_t N0_MM=read_djvec("mes_contr_P5P5_N0_R0_0_M0_R0_MASS",T).symmetrized();
  const djvec_t N0_MQ=read_djvec("mes_contr_P5P5_N0_R0_0_M0_R0_QED",T).symmetrized();
  const djvec_t NF_MF=read_djvec("mes_contr_P5P5_N0_R0_F_M0_R0_F",T).symmetrized();
  
  const djvec_t eff_pion=effective_mass(M0_M0);
  const djvec_t eff_kaon=effective_mass(N0_M0);
  
  const djack_t MP=constant_fit(ainv*eff_pion,10,18,"plots/pion.xmg");
  const djack_t MK=constant_fit(ainv*eff_kaon,10,18,"plots/kaon.xmg");
  
  const djack_t Mss=sqrt(2*MK*MK-MP*MP);
  cout<<" Mss: "<<smart_print(Mss)<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t K_IB_corr=-N0_MM;
  const djvec_t K_IB_ratio=K_IB_corr/M0_N0;
  const djvec_t K_IB_slope=-effective_slope(K_IB_ratio,eff_kaon,TH);
  const djack_t dMK_IB=constant_fit(K_IB_slope,10,18,"plots/kaon_ib_corr.xmg");
  
  const double ZP=0.529;
  const djack_t dM2K_IB=dMK_IB*2*MK*ZP;
  cout<<"dM2k_ib="<<dM2K_IB.ave_err()<<" fig.12 https://arxiv.org/pdf/1704.06561.pdf ok "<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  djack_t ZP_fact;
  ZP_fact.fill_gauss({1.5,0.2,23492});
  
  const double aml=0.0100;
  const djvec_t K_QED_corr=N0_MQ*(eu*eu-ed*ed)+NF_MF*es*(eu-ed)-N0_MM*aml*(sqr(ed)-sqr(eu))*Z_QED*ZP_fact;
  const djvec_t a=N0_MM*aml*(sqr(ed)-sqr(eu))*Z_QED*ZP_fact;
  const djvec_t b=(N0_MQ*(eu*eu-ed*ed)+NF_MF*es*(eu-ed));
  
  effective_slope(djvec_t(a/M0_N0),eff_kaon,TH).ave_err().write("plots/a.xmg");
  effective_slope(djvec_t(b/M0_N0),eff_kaon,TH).ave_err().write("plots/b.xmg");
  
  const djvec_t K_QED_ratio=K_QED_corr/M0_N0;
  const djvec_t K_QED_slope=-ainv*effective_slope(K_QED_ratio,eff_kaon,TH)*e2;
  
  const djack_t dMK_QED=constant_fit(K_QED_slope,10,18,"plots/kaon_qed_corr.xmg");
  
  const djack_t dM2K_QED=dMK_QED*2*MK;
  cout<<"dM2k="<<dM2K_QED.ave_err()<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t P_QED_corr=MF_MF*sqr(eu-ed);
  const djvec_t P_QED_ratio=-P_QED_corr/M0_M0;

  // djack_t Z2,M,DZ2_fr_Z2,SL;
  // two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL,M0_M0,P_QED_corr,TH,10,18,"plots/a.xmg","plots/b.xmg");
  
  const djvec_t P_QED_slope=-ainv*effective_slope(P_QED_ratio,eff_pion,TH)*e2;
  
  const djack_t dMP_QED=constant_fit(P_QED_slope,10,18,"plots/pion_qed_corr.xmg");
  
  const djack_t dM2P_QED=dMP_QED*MP;
  cout<<"dM2p="<<dM2P_QED.ave_err()<<" fig.6 https://arxiv.org/pdf/1704.06561.pdf ok "<<endl;
  
  return 0;
  
}
