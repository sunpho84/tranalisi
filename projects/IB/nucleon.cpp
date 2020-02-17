#include <tranalisi.hpp>

const size_t T=48,TH=T/2;

const double qu=+2.0/3,qd=-1.0/3;

const string jacks_path="jacks_16";
const string plots_path="plots_16";

djvec_t read_dir(const string& what)
{
  return read_djvec(combine("%s/bar_contr_dir_%s",jacks_path.c_str(),what.c_str()),T).symmetrized(-1);
}

void write_plot(const string& name,const djvec_t& a)
{
  a.ave_err().write(combine("%s/bar_contr_dir_%s.xmg",plots_path.c_str(),name.c_str()));
}

int main()
{
  //udu
  
  set_njacks(15);
  
  const djvec_t dir_000=read_dir("000");
  
  const djvec_t eff=effective_mass(dir_000);
  write_plot("eff",eff);
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t dir_Q00=read_dir("Q00")/dir_000;
  const djvec_t dir_0Q0=read_dir("0Q0")/dir_000;
  const djvec_t dir_00Q=read_dir("00Q")/dir_000;
  
  const djvec_t dir_Q=(dir_Q00+dir_00Q-dir_0Q0)*(qd*qd-qu*qu);
  
  write_plot("Q00",dir_Q00);
  write_plot("0Q0",dir_0Q0);
  write_plot("00Q",dir_00Q);
  
  write_plot("Q",dir_Q);
  
  const djvec_t eff_slope_Q=effective_slope(dir_Q,eff,TH);
  write_plot("slope_Q",eff_slope_Q);
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t dir_0FF=read_dir("0FF")/dir_000;
  const djvec_t dir_F0F=read_dir("F0F")/dir_000;
  const djvec_t dir_FF0=read_dir("FF0")/dir_000;
  
  const djvec_t dir_FF=dir_F0F*(qd*qd-qu*qu);
  
  write_plot("0FF",dir_0FF);
  write_plot("F0F",dir_F0F);
  write_plot("FF0",dir_FF0);
  
  write_plot("FF",dir_FF);
  
  const djvec_t eff_slope_FF=effective_slope(dir_FF,eff,TH);
  write_plot("slope_FF",eff_slope_FF);
  
  const djvec_t corr_QED=(eff_slope_FF+eff_slope_Q)*e2*2.2;
  write_plot("corr_QED",corr_QED);
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t dir_S00=read_dir("S00")/dir_000;
  const djvec_t dir_0S0=read_dir("0S0")/dir_000;
  const djvec_t dir_00S=read_dir("00S")/dir_000;
  
  const djvec_t dir_S=(-dir_0S0+dir_S00+dir_00S)*2;
  
  write_plot("S00",dir_S00);
  write_plot("0S0",dir_0S0);
  write_plot("00S",dir_00S);
  
  write_plot("S",dir_S);
  
  const djvec_t eff_slope_S=effective_slope(dir_S,eff,TH);
  write_plot("slope_S",eff_slope_S);
  
  
  return 0;
}
