#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

size_t T,TH,tmin;
double kappa;
constexpr size_t _njacks=15;

djvec_t read(const string &name,const bool REIM,int par)
{
  return read_djvec("data/mes_contr_"+name,T,REIM).symmetrized(par);
}

int main()
{
  set_njacks(_njacks);
  
  input_file_t input("input.txt");
  T=input.read<size_t>("T");
  TH=T/2;
  tmin=input.read<size_t>("Tmin");
  kappa=input.read<double>("Kappa");
  
  djvec_t V0P5=read("S0__S0_V0P5",IM,-1);
  djvec_t P5P5=read("S0__S0_P5P5",RE,1);
  djvec_t V0_P_P5=read("S0_P_S0__S0_V0P5",RE,-1);
  
  V0P5.ave_err().write("plots/V0P5.xmg");
  V0_P_P5.ave_err().write("plots/V0_P_P5.xmg");
  
  //pion mass
  djack_t MPi=constant_fit(effective_mass(P5P5),tmin,TH,"plots/P5P5.xmg");
  cout<<"MPi: "<<smart_print(MPi)<<endl;
  
  //m_cr
  djvec_t m_cr_corr_fw=forward_derivative(V0P5)/(2.0*P5P5);
  djvec_t m_cr_corr_symm=(forward_derivative(V0P5)+backward_derivative(V0P5))/(4.0*P5P5);
  //fit
  djack_t m_cr_fw=constant_fit(m_cr_corr_fw,tmin,TH,"plots/m_cr_fw_t.xmg");
  djack_t m_cr_symm=constant_fit(m_cr_corr_symm,tmin,TH,"plots/m_cr_symm_t.xmg");
  djack_t m_cr_symm_p2=constant_fit(m_cr_corr_symm,tmin+2,TH,"plots/m_cr_symm_t_p2.xmg");
  //out
  cout<<"m_cr_fw: "<<smart_print(m_cr_fw.ave_err())<<endl;
  cout<<"m_cr_symm: "<<smart_print(m_cr_symm.ave_err())<<endl;
  cout<<"m_cr_symm_p2: "<<smart_print(m_cr_symm_p2.ave_err())<<endl;
  
  djvec_t num_deltam_cr=V0P5;
  djvec_t den_deltam_cr=V0_P_P5;
  djvec_t deltam_cr_corr=num_deltam_cr/(2.0*den_deltam_cr);
  djack_t deltam_cr=constant_fit(deltam_cr_corr,tmin,TH,"plots/deltam_cr_t.xmg");
  
  cout<<"deltam_cr: "<<smart_print(deltam_cr.ave_err())<<endl;
  
  djack_t kappa_true=kappa/(1+2*deltam_cr*kappa);
  cout<<"true kappa: "<<smart_print(kappa_true.ave_err())<<" old: "<<kappa<<endl;
  
  return 0;
}
