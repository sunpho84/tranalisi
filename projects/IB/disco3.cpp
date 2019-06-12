#include <tranalisi.hpp>

int T,TH,L;
int vol;

int main(int narg,char **arg)
{
  string name="input.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  T=input.read<int>("T");
  TH=T/2;
  L=input.read<int>("L");
  vol=T*L*L*L;
  
  range_t conf_range;
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  set_njacks(15);
  
  const char temp[]="jacks/mes_contr_%s_R0_0_%s_R0_0_P5P5";
  const djvec_t l_l=read_djvec(combine(temp,"S_M0","S_M0"),T).symmetrized();
  const djvec_t s_l=read_djvec(combine(temp,"S_M1","S_M0"),T).symmetrized();
  const djvec_t l_s=read_djvec(combine(temp,"S_M0","S_M1"),T).symmetrized();
  const djvec_t s_s=read_djvec(combine(temp,"S_M1","S_M1"),T).symmetrized();
  const djvec_t sFs_lFl=read_djvec(combine(temp,"SFS_M1","SFS_M0"),T).symmetrized();
  const djvec_t lFl_sFs=read_djvec(combine(temp,"SFS_M0","SFS_M1"),T).symmetrized();
  const djvec_t l3l_l=read_djvec(combine(temp,"SXS_M0","S_M0"),T).symmetrized();
  const djvec_t s3s_l=read_djvec(combine(temp,"SXS_M1","S_M0"),T).symmetrized();
  const djvec_t l3l_s=read_djvec(combine(temp,"SXS_M0","S_M1"),T).symmetrized();
  const djvec_t s3s_s=read_djvec(combine(temp,"SXS_M1","S_M1"),T).symmetrized();
  const djvec_t s_l3l=read_djvec(combine(temp,"S_M1","SXS_M0"),T).symmetrized();
  const djvec_t l_s3s=read_djvec(combine(temp,"S_M0","SXS_M1"),T).symmetrized();
  
  const djvec_t r_s3s_l=s3s_l/s_l;
  const djvec_t r_s_l3l=s_l3l/s_l;
  // const djvec_t rexc=sFs_sFs/s_s;
  djack_t Z2,M,DZ2_fr_Z2_3,SL_3,DZ2_fr_Z2_exc,SL_exc;
  s_l3l.ave_err().write("plots/s_l3l.xmg");
  
  l3l_l.ave_err().write("plots/l3l_l.xmg");
  s3s_s.ave_err().write("plots/s3s_s.xmg");
  
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2_3,SL_3,s_l,s_l3l,TH,10,TH,"plots/s_l.xmg","plots/r_s_l3l.xmg");
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2_exc,SL_exc,s_l,sFs_lFl,TH,10,TH,"plots/s_l.xmg","plots/r_sFs_lFl.xmg");
  // two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2_3,SL_3,l_l,l3l_l,TH,10,TH,"plots/l_l.xmg","plots/r_l3l_l.xmg");
  // two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2_3,SL_3,s_s,s3s_s,TH,10,TH,"plots/s_s.xmg","plots/r_s3s_s.xmg");
  
  djack_t SL_tot=SL_exc/3-SL_3;
  
  // cout<<SL_3.ave_err()<<" "<<SL_3.err()/SL_3.ave()<<endl;
  // cout<<SL_exc.ave_err()<<" "<<SL_exc.err()/SL_exc.ave()<<endl;
  // cout<<SL_tot.ave_err()<<" "<<SL_tot.err()/SL_tot.ave()<<endl;

  const djack_t r=SL_3*3/SL_exc;
  cout<<r.ave_err()<<endl;
  
  // r_s3s_l.ave_err().write("r_s3s_l.xmg");
  // r_s_l3l.ave_err().write("r_s_l3l.xmg");
  // r3.ave_err().write("/tmp/r3.xmg");
  // rexc.ave_err().write("/tmp/rexc.xmg");
  
  return 0;
}
