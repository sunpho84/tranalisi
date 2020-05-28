#include <tranalisi.hpp>

double adml,ainv;
size_t T,TH;

const double qu=+2.0/3,qd=-1.0/3;

string jacks_path;
string plots_path;

size_t tfit_beg,tfit_end;

djvec_t read_bar_contr(const string& what)
{
  return read_djvec(combine("%s/bar_contr_%s",jacks_path.c_str(),what.c_str()),T).symmetrized(-1);
}

djvec_t read_bar(const string& what)
{
  return read_bar_contr("dir_"+what)-read_bar_contr("exc_"+what);
}

void write_plot(const string& name,const djvec_t& a)
{
  a.ave_err().write(combine("%s/bar_contr_%s.xmg",plots_path.c_str(),name.c_str()));
}

int main()
{
  //udu
  
  input_file_t input("input");
  const size_t ext_njacks=input.read<size_t>("NJacks");
  
  jacks_path=input.read<string>("JacksPath");
  plots_path=input.read<string>("PlotsPath");
  T=input.read<size_t>("T");
  adml=input.read<double>("adml");
  ainv=input.read<double>("ainv");
  tfit_beg=input.read<double>("tfitBeg");
  tfit_end=input.read<double>("tfitEnd");
  
  TH=T/2;
  set_njacks(ext_njacks);
  
  const djvec_t bar_000=read_bar("000");
  write_plot("000",bar_000);
  
  const djvec_t eff=effective_mass(bar_000);
  write_plot("eff",eff);
  
  const djvec_t P5P5=read_djvec(jacks_path+"/mes_contr_P5P5_M0_R0_0_M0_R0_0",T).symmetrized(1);
  const djvec_t V0P5=read_djvec(jacks_path+"/mes_contr_V0P5_M0_R0_0_M0_R0_0",T).symmetrized(-1);
  const djvec_t mPCAC_corr=forward_derivative(V0P5)/P5P5;
  mPCAC_corr.ave_err().write(plots_path+"/mPCAC.xmg");
  
  const djvec_t P5P5_QED=read_djvec(jacks_path+"/mes_contr_P5P5_M0_R0_0_M0_R0_QED",T).symmetrized(1);
  const djvec_t V0P5_QED=read_djvec(jacks_path+"/mes_contr_V0P5_M0_R0_0_M0_R0_QED",T).symmetrized(-1);
  V0P5_QED.ave_err().write(plots_path+"/V0P5_QED.xmg");
  
  const djvec_t deltam_cr=forward_derivative(V0P5_QED)/P5P5-forward_derivative(V0P5)/sqr(P5P5)*P5P5_QED;
  deltam_cr.ave_err().write(plots_path+"/deltam_cr.xmg");
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t bar_Q00=read_bar("Q00")/bar_000;
  const djvec_t bar_0Q0=read_bar("0Q0")/bar_000;
  const djvec_t bar_00Q=read_bar("00Q")/bar_000;
  
  const djvec_t bar_Q=(bar_Q00+bar_00Q-bar_0Q0)*(qd*qd-qu*qu);
  
  write_plot("Q00",bar_Q00);
  write_plot("0Q0",bar_0Q0);
  write_plot("00Q",bar_00Q);
  
  write_plot("Q",bar_Q);
  
  const djvec_t eff_slope_Q=-effective_slope(bar_Q,eff,TH);
  write_plot("slope_Q",eff_slope_Q);
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t bar_0FF=read_bar("0FF")/bar_000;
  const djvec_t bar_F0F=read_bar("F0F")/bar_000;
  const djvec_t bar_FF0=read_bar("FF0")/bar_000;
  
  const djvec_t bar_FF=bar_F0F*(qd*qd-qu*qu);
  
  write_plot("0FF",bar_0FF);
  write_plot("F0F",bar_F0F);
  write_plot("FF0",bar_FF0);
  
  write_plot("FF",bar_FF);
  
  const djvec_t eff_slope_FF=-effective_slope(bar_FF,eff,TH);
  write_plot("slope_FF",eff_slope_FF);
  
  const djvec_t romiti=-(eff_slope_FF+eff_slope_Q);
  write_plot("romiti",romiti);
  
  const djvec_t corr_QED=(eff_slope_FF+eff_slope_Q)*e2*ainv*1e3;
  const djack_t NP_QED=constant_fit(corr_QED,
				    tfit_beg,tfit_end,plots_path+"/corr_QED.xmg");
  cout<<"QED: "<<smart_print(NP_QED)<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t bar_S00=read_bar("S00")/bar_000;
  const djvec_t bar_0S0=read_bar("0S0")/bar_000;
  const djvec_t bar_00S=read_bar("00S")/bar_000;
  
  const djvec_t bar_S=(-bar_0S0+bar_S00+bar_00S);
  
  write_plot("S00",bar_S00);
  write_plot("0S0",bar_0S0);
  write_plot("00S",bar_00S);
  
  write_plot("S",bar_S);
  
  cout<<" dml: "<<adml*ainv*1e3/0.423<<endl;
  
  const djvec_t eff_slope_S=-effective_slope(bar_S,eff,TH)*adml*ainv*1e3;
  write_plot("slope_S",eff_slope_S);
  const djack_t test=constant_fit((djvec_t)(-effective_slope(bar_S,eff,TH)*0.423),tfit_beg,tfit_end,plots_path+"/test.xmg");
  cout<<"Test: "<<smart_print(test)<<endl;
  
  
  const djack_t corr_IB=constant_fit(eff_slope_S,tfit_beg,tfit_end,plots_path+"/corr_IB.xmg");
  cout<<"IB: "<<corr_IB.ave_err()<<" MeV"<<endl;
  
  const djvec_t corr_TOT=eff_slope_S+corr_QED;
  
  const djack_t NP_diff=constant_fit(corr_TOT,tfit_beg,tfit_end);
  grace_file_t NP_diff_plot(plots_path+"/corr_TOT.xmg");
  NP_diff_plot.write_constant_band(tfit_beg,tfit_end,NP_diff);
  NP_diff_plot.write_vec_ave_err(corr_TOT.ave_err());
  const double NP_diff_phys=1.293;
  NP_diff_plot.write_line([NP_diff_phys](double x){return NP_diff_phys;},0,TH,grace::BLACK);
  cout<<"Total: "<<smart_print(NP_diff)<<" MeV    (exp: "<<NP_diff_phys<<" MeV)"<<endl;
  
  return 0;
}
