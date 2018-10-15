#include <tranalisi.hpp>

size_t ib;
size_t NT;
size_t TH;
size_t tmin;
size_t tmax;

double a_GeVMOne[3]={0.4497,0.4137,0.3142};

void set_NT(int ext_NT)
{
  NT=ext_NT;
  TH=NT/2;
  tmin=TH/4;
  tmax=TH;
}

enum ins_t{_,P,S,T,F,FF};
const string ins_tag[6]={"","_P","_S","_T","_F","_FF"};
enum bil_t{P5P5,V0P5};
const string bil_tag[2]={"P5P5","V0P5"};

const double tau3[2]={-1.0,1.0};

dcompl_t coeff_to_read_ins(ins_t ins,const size_t r)
{
  switch(ins)
    {
    case P:
      return dcompl_t(0,tau3[r]);
      break;
    case S:
      return -1.0;
      break;
    default:
      return 1.0;
    }
}

dcompl_t coeff_to_read_bil(const bil_t bil,const size_t r)
{
  const dcompl_t ri[2]={{1.0,0.0},{0.0,tau3[r]}};
  
  return ri[bil];
}

djvec_t load(const size_t rbw,const ins_t ins_bw,const size_t rfw,const ins_t ins_fw,const bil_t bil)
{
  const dcompl_t coeff=conj(coeff_to_read_ins(ins_bw,rbw))*coeff_to_read_ins(ins_fw,rfw)*coeff_to_read_bil(bil,rfw);
  int par=1;
  
  if(bil==V0P5) par=-par;
  
  string name="mes_contr";
  name+="_S";
  name+=to_string(rbw);
  name+=ins_tag[ins_bw];
  name+="_S";
  name+=to_string(rfw);
  name+=ins_tag[ins_fw];
  name+="_"+bil_tag[bil];
  
  int i;
  double d;
  if(coeff.real())
    {
      i=0;
      d=coeff.real();
    }
  else
    {
      i=1;
      d=coeff.imag();
    }
  
  const djvec_t out=d*read_djvec("jacks/"+name,NT,i);//.symmetrized(par*0);
  
  out.ave_err().write("plots/"+name+".xmg");
  
  return out;
}

djvec_t load(const size_t rdiff,const ins_t ins_bw,const ins_t ins_fw,const bil_t bil)
{
  djvec_t out(TH+1);
  
  const size_t nr=2;
  for(size_t rbw=0;rbw<nr;rbw++)
    {
      const int rfw=(rbw+rdiff)%2;
      out+=load(rbw,ins_bw,rfw,ins_fw,bil);
    }
  out/=nr;
  
  string name="mes_contr_rdiff"+to_string(rdiff);
  name+=ins_tag[ins_bw];
  name+=ins_tag[ins_fw];
  name+="_"+bil_tag[bil];
  
  out.ave_err().write("plots/"+name+".xmg");
  
  return out;
}

int main()
{
  input_file_t input("analysis.txt");
  const size_t read_NT=input.read<size_t>("T");
  ib=input.read<size_t>("ib");
  set_NT(read_NT);
  input.close();
  
  set_njacks(15);
  
  const int SAME=0;
  const int OPPO=1;
  
  //load correlators with the same r for P5P5
  const djvec_t P5P5_SAME_LO=load(SAME,_,_,P5P5);
  const djvec_t P5P5_SAME_P=2.0*load(SAME,_,P,P5P5);
  const djvec_t P5P5_SAME_S=2.0*load(SAME,_,S,P5P5);
  const djvec_t P5P5_SAME_QED=2.0*(load(SAME,_,T,P5P5)+load(SAME,_,FF,P5P5))+load(SAME,F,F,P5P5);
  const djvec_t eff_P5P5_SAME_LO=effective_mass(P5P5_SAME_LO);
  //load correlators with opposite r for P5P5
  const djvec_t P5P5_OPPO_LO=load(OPPO,_,_,P5P5);
  const djvec_t P5P5_OPPO_P=2.0*load(OPPO,_,P,P5P5);
  const djvec_t P5P5_OPPO_S=2.0*load(OPPO,_,S,P5P5);
  const djvec_t P5P5_OPPO_QED=2.0*(load(OPPO,_,T,P5P5)+load(OPPO,_,FF,P5P5))+load(OPPO,F,F,P5P5);
  //load correlators with same r for V0P5
  const djvec_t V0P5_SAME_LO=load(SAME,_,_,V0P5);
  const djvec_t V0P5_SAME_P=2.0*load(SAME,_,P,V0P5);
  const djvec_t V0P5_SAME_S=2.0*load(SAME,_,S,V0P5);
  const djvec_t V0P5_SAME_QED=2.0*(load(SAME,_,T,V0P5)+load(SAME,_,FF,V0P5))+load(SAME,F,F,V0P5);
  
  //build the system of equations
  const djvec_t a=djvec_t(symmetric_derivative(V0P5_SAME_QED)/P5P5_SAME_LO-symmetric_derivative(V0P5_SAME_LO)/sqr(P5P5_SAME_LO)*P5P5_SAME_QED).subset(0,TH-1);
  const djvec_t b=djvec_t(symmetric_derivative(V0P5_SAME_S)/P5P5_SAME_LO-symmetric_derivative(V0P5_SAME_LO)/sqr(P5P5_SAME_LO)*P5P5_SAME_S).subset(0,TH-1);
  const djvec_t c=djvec_t(symmetric_derivative(V0P5_SAME_P)/P5P5_SAME_LO-symmetric_derivative(V0P5_SAME_LO)/sqr(P5P5_SAME_LO)*P5P5_SAME_P).subset(0,TH-1);
  const djvec_t d=effective_slope(djvec_t(P5P5_SAME_QED/P5P5_SAME_LO),eff_P5P5_SAME_LO,TH);
  const djvec_t e=effective_slope(djvec_t(P5P5_SAME_S/P5P5_SAME_LO),eff_P5P5_SAME_LO,TH);
  const djvec_t f=effective_slope(djvec_t(P5P5_SAME_P/P5P5_SAME_LO),eff_P5P5_SAME_LO,TH);

  //determine deltam ignoring mass insertion
  const djvec_t deltam_cr_corr_simple=-a/c;
  const djack_t deltam_cr_simple=constant_fit(deltam_cr_corr_simple,tmin,tmax,"plots/deltam_cr_corr_SAME_simple.xmg");
  cout<<"deltam_cr_simple: "<<smart_print(deltam_cr_simple.ave_err())<<endl;
  
  //solve the system
  const djvec_t den=b*f-c*e;
  const djvec_t deltam_tm_corr_full=(-a*f+c*d)/den;
  const djvec_t deltam_cr_corr_full=(-b*d+a*e)/den;
  
  //fit the coefficients
  const djack_t deltam_tm_full=constant_fit(deltam_tm_corr_full,tmin,tmax,"plots/deltam_tm_corr_SAME_full.xmg");
  const djack_t deltam_cr_full=constant_fit(deltam_cr_corr_full,tmin,tmax,"plots/deltam_cr_corr_SAME_full.xmg");
  cout<<"deltam_tm_full: "<<smart_print(deltam_tm_full.ave_err())<<endl;
  cout<<"deltam_cr_full: "<<smart_print(deltam_cr_full.ave_err())<<endl;
  
  //select simple or full calculation
  enum{SIMPLE,FULL};
  const bool simple_full=FULL;
  const djack_t deltam_tm=(simple_full==SIMPLE)?0.0:deltam_tm_full;
  const djack_t deltam_cr=(simple_full==SIMPLE)?deltam_cr_simple:deltam_cr_full;
  
  //build corrected correlators with same r
  const djvec_t P5P5_SAME_CORR=P5P5_SAME_QED+deltam_tm*P5P5_SAME_S+deltam_cr*P5P5_SAME_P;
  const djvec_t P5P5_SAME_REL_CORR=P5P5_SAME_CORR/P5P5_SAME_LO;
  P5P5_SAME_REL_CORR.ave_err().write("plots/P5P5_SAME_CORR.xmg");
  djack_t G2P_SAME,M_SAME,DG2_FR_G2_SAME,M_SL_SAME;
  two_pts_with_ins_ratio_fit(G2P_SAME,M_SAME,DG2_FR_G2_SAME,M_SL_SAME,P5P5_SAME_LO,P5P5_SAME_CORR,TH,tmin,tmax,"plots/P5P5_SAME_effmass.xmg","plots/P5P5_SAME_CORR_eff_slope.xmg");
  const djack_t DG2_FR_G2_SAME_constfit=constant_fit(P5P5_SAME_REL_CORR,tmin,tmax,"plots/P5P5_SAME_CORR_constant_fit.xmg");
  
  //build corrected correlators with opposite r
  const djvec_t P5P5_OPPO_CORR=P5P5_OPPO_QED+deltam_tm*P5P5_OPPO_S+deltam_cr*P5P5_OPPO_P;
  const djvec_t P5P5_OPPO_REL_CORR=P5P5_OPPO_CORR/P5P5_OPPO_LO;
  P5P5_OPPO_REL_CORR.ave_err().write("plots/P5P5_OPPO_CORR.xmg");
  djack_t G2P_OPPO,M_OPPO,DG2_FR_G2_OPPO,M_SL_OPPO;
  two_pts_with_ins_ratio_fit(G2P_OPPO,M_OPPO,DG2_FR_G2_OPPO,M_SL_OPPO,P5P5_OPPO_LO,P5P5_OPPO_CORR,TH,tmin,tmax,"plots/P5P5_OPPO_effmass.xmg","plots/P5P5_OPPO_CORR_eff_slope.xmg");
  
  //build the correlator representing the correction of Zs/Zp
  // const djvec_t dZs_fr_Zp_corr=(P5P5_SAME_REL_CORR-P5P5_OPPO_REL_CORR)/2.0;
  // const djack_t dZs_fr_Zp=constant_fit(dZs_fr_Zp_corr,tmin,tmax,"plots/dZs_fr_Zp.xmg");
  // cout<<"dZs_fr_Zp rough: "<<smart_print(dZs_fr_Zp)<<endl;
  
  const djack_t dZs_fr_Zp_fitted=(DG2_FR_G2_SAME-DG2_FR_G2_OPPO)/2.0;
  cout<<"dZs_fr_Zp fitted: "<<smart_print(dZs_fr_Zp_fitted)<<endl;
  cout<<M_SAME.ave()/a_GeVMOne[ib]<<" "<<dZs_fr_Zp_fitted<<endl;
  
  // const djack_t dZs_fr_Zp_fitted_with_same_const=(DG2_FR_G2_SAME_constfit-DG2_FR_G2_OPPO)/2.0;
  // cout<<"dZs_fr_Zp fitted with \"same\" constant: "<<smart_print(dZs_fr_Zp_fitted_with_same_const)<<endl;
  
  return 0;
}
