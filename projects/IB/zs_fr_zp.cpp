#include <tranalisi.hpp>

const int NT=48;
const int TH=NT/2;

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
  
  set_njacks(15);
  
  const int SAME=0;
  const int OPPO=1;
  
  //load correlators with the same r for P5P5
  const djvec_t P5P5_SAME_LO=load(SAME,_,_,P5P5);
  const djvec_t P5P5_SAME_P=2.0*load(SAME,_,P,P5P5);
  const djvec_t P5P5_SAME_S=2.0*load(SAME,_,S,P5P5);
  const djvec_t P5P5_SAME_QED=2.0*(load(SAME,_,T,P5P5)+load(SAME,_,FF,P5P5))+load(SAME,F,F,P5P5);
  const djvec_t eff_P5P5_LO=effective_mass(P5P5_SAME_LO);
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
  const djvec_t d=effective_slope(djvec_t(P5P5_SAME_QED/P5P5_SAME_LO),eff_P5P5_LO,TH);
  const djvec_t e=effective_slope(djvec_t(P5P5_SAME_S/P5P5_SAME_LO),eff_P5P5_LO,TH);
  const djvec_t f=effective_slope(djvec_t(P5P5_SAME_P/P5P5_SAME_LO),eff_P5P5_LO,TH);
  
  //solve the system
  const djvec_t den=b*f-c*e;
  const djvec_t deltam_tm_corr=(-a*f+c*d)/den;
  const djvec_t deltam_cr_corr=(-b*d+a*e)/den;
  
  //fit the coefficients
  const djack_t deltam_cr=constant_fit(deltam_cr_corr,TH/4,TH*3/4,"plots/deltam_cr_corr_SAME.xmg");
  const djack_t deltam_tm=constant_fit(deltam_tm_corr,TH/4,TH*3/4,"plots/deltam_tm_corr_SAME.xmg");
  cout<<"deltam_cr: "<<smart_print(deltam_cr.ave_err())<<endl;
  cout<<"deltam_tm: "<<smart_print(deltam_tm.ave_err())<<endl;
  
  //build corrected correlators with same r
  const djvec_t P5P5_SAME_CORR=P5P5_SAME_QED+deltam_tm*P5P5_SAME_S+deltam_cr*P5P5_SAME_P;
  const djvec_t P5P5_SAME_REL_CORR=P5P5_SAME_CORR/P5P5_SAME_LO;
  P5P5_SAME_REL_CORR.ave_err().write("plots/P5P5_SAME_CORR.xmg");
  
  //build corrected correlators with opposite r
  const djvec_t P5P5_OPPO_CORR=P5P5_OPPO_QED+deltam_tm*P5P5_OPPO_S+deltam_cr*P5P5_OPPO_P;
  const djvec_t P5P5_OPPO_REL_CORR=P5P5_OPPO_CORR/P5P5_OPPO_LO;
  P5P5_OPPO_REL_CORR.ave_err().write("plots/P5P5_OPPO_CORR.xmg");
  
  //build the correlator representing the correction of Zs/Zp
  const djvec_t dZs_fr_Zp_corr=(P5P5_SAME_REL_CORR-P5P5_OPPO_REL_CORR)/2.0;
  const djack_t dZs_fr_Zp=constant_fit(dZs_fr_Zp_corr,TH/4,TH*3/4,"plots/dZs_fr_Zp.xmg");
  cout<<"dZs_fr_Zp: "<<smart_print(dZs_fr_Zp)<<endl;
  
  return 0;
}
