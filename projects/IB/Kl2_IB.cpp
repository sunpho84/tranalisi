#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int T;
index_t ind;
index_t ind_hl;

djvec_t load(string path,size_t im1,size_t im2,size_t reim,int rpar,int spat_par)
{
  djvec_t corr_r0=read_djvec(path,T,ind({im1,im2,0,reim}));
  djvec_t corr_r1=read_djvec(path,T,ind({im1,im2,1,reim}));
  
  return djvec_t(corr_r0+rpar*corr_r1).symmetrized(spat_par)/(1.0+abs(rpar));
}

djvec_t load_P5P5(string ins,size_t im1,size_t im2,size_t reim=RE,int rpar=1,int spat_par=1)
{return load(combine("corr%s_P5P5",ins.c_str()),im1,im2,reim,rpar,spat_par);}

djvec_t load_V0P5(string ins,size_t im1,size_t im2,size_t reim=IM,int rpar=-1,int spat_par=-1)
{return load(combine("corr%s_V0P5",ins.c_str()),im1,im2,reim,rpar,spat_par);}

int main(int narg,char **arg)
{
  if(narg<3) CRASH("Use %s T ibeta",arg[0]);

  int ibeta;
  T=atoi(arg[1]);
  ibeta=atoi(arg[2]);
  int TH=T/2;
  
  ind.set_ranges({3,3,2,2});
  const size_t nrl=2,norie=2,nrq=2,nqrev=2,nqins=3,nleptons=6;
  ind_hl.set_ranges({nrl,norie,nrq,nqrev,nqins,nleptons,2/*reim*/});
  set_njacks(15);
  
  raw_file_t input;
  
  switch(ibeta)
    {
    case 0: input.open("../../input_t_A_32.txt","r"); break;
    case 1: input.open("../../input_t_A_24.txt","r"); break;
    case 2: input.open("../../input_t_A_20.txt","r"); break;
    case 3: input.open("../../input_t_B_32.txt","r"); break;
    case 4: input.open("../../input_t_B_24.txt","r"); break;
    case 5: input.open("../../input_t_D.txt","r"); break;
    default: CRASH("Undefined ibeta"); break;
    }
 
  size_t pi_tmin=input.read<int>("pi_tmin");
  size_t pi_tmax=input.read<int>("pi_tmax");
  size_t k_tmin=input.read<int>("k_tmin");
  size_t k_tmax=input.read<int>("k_tmax");
  size_t D_tmin=input.read<int>("D_tmin");
  size_t D_tmax=input.read<int>("D_tmax");
  size_t Ds_tmin=input.read<int>("Ds_tmin");
  size_t Ds_tmax=input.read<int>("Ds_tmax");
  size_t ud_tmin=input.read<int>("ud_tmin");
  size_t ud_tmax=input.read<int>("ud_tmax");
  size_t s_tmin=input.read<int>("s_tmin");
  size_t s_tmax=input.read<int>("s_tmax");
  size_t c_tmin=input.read<int>("c_tmin");
  size_t c_tmax=input.read<int>("c_tmax");
  
  //pion

  djvec_t pi_plain=load_P5P5("00",0,0);
  djvec_t pi_exch=load_P5P5("LL",0,0);
  djvec_t pi_self=load_P5P5("0M",0,0);
  djvec_t pi_tad=load_P5P5("0T",0,0);
  djvec_t pi_s=load_P5P5("0S",0,0);
  djvec_t pi_p=-load_P5P5("0P",0,0,IM,-1);
  
  djack_t pi_Z,pi_M,pi_A_exch,pi_SL_exch,pi_A_selftad,pi_SL_selftad,pi_A_s,pi_SL_s,pi_A_p,pi_SL_p;
  two_pts_with_ins_ratio_fit(pi_Z,pi_M,pi_A_s,pi_SL_s,pi_plain,pi_s,TH,pi_tmin,pi_tmax,"pi_mass.xmg","pi_s.xmg");
  two_pts_with_ins_ratio_fit(pi_Z,pi_M,pi_A_p,pi_SL_p,pi_plain,pi_p,TH,pi_tmin,pi_tmax,"pi_mass.xmg","pi_p.xmg");
  two_pts_with_ins_ratio_fit(pi_Z,pi_M,pi_A_selftad,pi_SL_selftad,pi_plain,djvec_t(pi_self+pi_tad),TH,pi_tmin,pi_tmax,"pi_mass.xmg","pi_selftad.xmg");
  two_pts_with_ins_ratio_fit(pi_Z,pi_M,pi_A_exch,pi_SL_exch,pi_plain,pi_exch,TH,pi_tmin,pi_tmax,"pi_mass.xmg","pi_exch.xmg");

  cout<<"pi_M: "<<pi_M.ave_err()<<endl;
  cout<<"pi_A_exch: "<<pi_A_exch.ave_err()<<endl;
  cout<<"pi_SL_exch: "<<pi_SL_exch.ave_err()<<endl;
  cout<<"pi_A_selftad: "<<pi_A_selftad.ave_err()<<endl;
  cout<<"pi_SL_selftad: "<<pi_SL_selftad.ave_err()<<endl;
  cout<<"pi_A_s: "<<pi_A_s.ave_err()<<endl;
  cout<<"pi_SL_s: "<<pi_SL_s.ave_err()<<endl;
  cout<<"pi_A_p: "<<pi_A_p.ave_err()<<endl;
  cout<<"pi_SL_p: "<<pi_SL_p.ave_err()<<endl;

  raw_file_t pi_obs_file("pi_obs","w");
  pi_obs_file.bin_write(pi_M);
  pi_obs_file.bin_write(pi_SL_exch);
  pi_obs_file.bin_write(pi_SL_selftad);
  pi_obs_file.bin_write(pi_SL_s);
  pi_obs_file.bin_write(pi_SL_p);
  
  //kaon

  djvec_t k_plain=load_P5P5("00",0,1);
  djvec_t k_exch=load_P5P5("LL",0,1);
  djvec_t k_self=load_P5P5("0M",1,0);
  djvec_t k_tad=load_P5P5("0T",1,0);
  djvec_t k_self_revins=load_P5P5("0M",0,1);
  djvec_t k_tad_revins=load_P5P5("0T",0,1);
  djvec_t k_s=load_P5P5("0S",1,0);
  djvec_t k_p=-load_P5P5("0P",1,0,IM,-1);
  djvec_t k_p_revins=-load_P5P5("0P",0,1,IM,-1);
  
  djack_t k_Z,k_M,k_A_exch,k_SL_exch,k_A_selftad,k_SL_selftad,k_A_selftad_revins,k_SL_selftad_revins,k_A_s,k_SL_s,k_A_p,k_SL_p,k_A_p_revins,k_SL_p_revins;
  
  // djack_t M_k,Z_k;

  // djvec_t k_corr_00=load_P5P5("00",1,0);
  // two_pts_true_fit(Z_k,M_k,k_corr_00,TH,k_tmin,k_tmax,"Mk.xmg");
  // cout<<"k_Z: "<<Z_k.ave_err()<<endl;
  // cout<<"k_M: "<<M_k.ave_err()<<endl;
  
  two_pts_with_ins_ratio_fit(k_Z,k_M,k_A_exch,k_SL_exch,k_plain,k_exch,TH,k_tmin,k_tmax,"kaon_mass.xmg","k_exch.xmg");
  two_pts_with_ins_ratio_fit(k_Z,k_M,k_A_s,k_SL_s,k_plain,k_s,TH,k_tmin,k_tmax,"kaon_mass.xmg","k_s.xmg");
  two_pts_with_ins_ratio_fit(k_Z,k_M,k_A_p,k_SL_p,k_plain,k_p,TH,k_tmin,k_tmax,"kaon_mass.xmg","k_p.xmg");
  two_pts_with_ins_ratio_fit(k_Z,k_M,k_A_selftad_revins,k_SL_selftad_revins,k_plain,djvec_t(k_self_revins+k_tad_revins),TH,k_tmin,k_tmax,"kaon_mass.xmg","k_selftad_revins.xmg");
  two_pts_with_ins_ratio_fit(k_Z,k_M,k_A_p_revins,k_SL_p_revins,k_plain,k_p_revins,TH,k_tmin,k_tmax,"kaon_mass.xmg","k_p_revins.xmg");
  two_pts_with_ins_ratio_fit(k_Z,k_M,k_A_selftad,k_SL_selftad,k_plain,djvec_t(k_self+k_tad),TH,k_tmin,k_tmax,"kaon_mass.xmg","k_selftad.xmg");
  
  cout<<"k_M: "<<k_M.ave_err()<<endl;
  cout<<"k_A_exch: "<<k_A_exch.ave_err()<<endl;
  cout<<"k_SL_exch: "<<k_SL_exch.ave_err()<<endl;
  cout<<"k_A_selftad: "<<k_A_selftad.ave_err()<<endl;
  cout<<"k_SL_selftad: "<<k_SL_selftad.ave_err()<<endl;
  cout<<"k_A_selftad_revins: "<<k_A_selftad_revins.ave_err()<<endl;
  cout<<"k_SL_selftad_revins: "<<k_SL_selftad_revins.ave_err()<<endl;
  cout<<"k_A_s: "<<k_A_s.ave_err()<<endl;
  cout<<"k_SL_s: "<<k_SL_s.ave_err()<<endl;
  cout<<"k_A_p: "<<k_A_p.ave_err()<<endl;
  cout<<"k_SL_p: "<<k_SL_p.ave_err()<<endl;
  cout<<"k_A_p_revins: "<<k_A_p_revins.ave_err()<<endl;
  cout<<"k_SL_p_revins: "<<k_SL_p_revins.ave_err()<<endl;

  raw_file_t k_obs_file("k_obs","w");
  k_obs_file.bin_write(k_M);
  k_obs_file.bin_write(k_SL_exch);
  k_obs_file.bin_write(k_SL_selftad);
  k_obs_file.bin_write(k_SL_s);
  k_obs_file.bin_write(k_SL_p);
  k_obs_file.bin_write(k_SL_selftad_revins);
  k_obs_file.bin_write(k_SL_p_revins);

  //D meson

  djvec_t D_plain=load_P5P5("00",2,0);
  djvec_t D_exch=load_P5P5("LL",0,2);
  djvec_t D_self=load_P5P5("0M",2,0);
  djvec_t D_tad=load_P5P5("0T",2,0);
  djvec_t D_self_revins=load_P5P5("0M",0,2);
  djvec_t D_tad_revins=load_P5P5("0T",0,2);
  djvec_t D_s=load_P5P5("0S",2,0);
  djvec_t D_p=-load_P5P5("0P",2,0,IM,-1);
  djvec_t D_p_revins=-load_P5P5("0P",0,2,IM,-1);

  djack_t D_Z,D_M,D_A_exch,D_SL_exch,D_A_selftad,D_SL_selftad,D_A_selftad_revins,D_SL_selftad_revins,D_A_s,D_SL_s,D_A_p,D_SL_p,D_A_p_revins,D_SL_p_revins;
  
  // djack_t M_D,Z_D;
  
  // djvec_t D_corr_00=load_P5P5("00",2,0);
  // two_pts_true_fit(Z_D,M_D,D_corr_00,TH,D_tmin,D_tmax,"MD.xmg");
  // cout<<"D_Z: "<<Z_D.ave_err()<<endl;
  // cout<<"D_M: "<<M_D.ave_err()<<endl;
  
  two_pts_with_ins_ratio_fit(D_Z,D_M,D_A_exch,D_SL_exch,D_plain,D_exch,TH,D_tmin,D_tmax,"D_mass.xmg","D_exch.xmg");
  two_pts_with_ins_ratio_fit(D_Z,D_M,D_A_s,D_SL_s,D_plain,D_s,TH,D_tmin,D_tmax,"D_mass.xmg","D_s.xmg");
  two_pts_with_ins_ratio_fit(D_Z,D_M,D_A_p,D_SL_p,D_plain,D_p,TH,D_tmin,D_tmax,"D_mass.xmg","D_p.xmg");
  two_pts_with_ins_ratio_fit(D_Z,D_M,D_A_selftad_revins,D_SL_selftad_revins,D_plain,djvec_t(D_self_revins+D_tad_revins),TH,D_tmin,D_tmax,"D_mass.xmg","D_selftad_revins.xmg");
  two_pts_with_ins_ratio_fit(D_Z,D_M,D_A_p_revins,D_SL_p_revins,D_plain,D_p_revins,TH,D_tmin,D_tmax,"D_mass.xmg","D_p_revins.xmg");
  two_pts_with_ins_ratio_fit(D_Z,D_M,D_A_selftad,D_SL_selftad,D_plain,djvec_t(D_self+D_tad),TH,D_tmin,D_tmax,"D_mass.xmg","D_selftad.xmg");
  
  cout<<"D_M: "<<D_M.ave_err()<<endl;
  cout<<"D_A_exch: "<<D_A_exch.ave_err()<<endl;
  cout<<"D_SL_exch: "<<D_SL_exch.ave_err()<<endl;
  cout<<"D_A_selftad: "<<D_A_selftad.ave_err()<<endl;
  cout<<"D_SL_selftad: "<<D_SL_selftad.ave_err()<<endl;
  cout<<"D_A_selftad_revins: "<<D_A_selftad_revins.ave_err()<<endl;
  cout<<"D_SL_selftad_revins: "<<D_SL_selftad_revins.ave_err()<<endl;
  cout<<"D_A_s: "<<D_A_s.ave_err()<<endl;
  cout<<"D_SL_s: "<<D_SL_s.ave_err()<<endl;
  cout<<"D_A_p: "<<D_A_p.ave_err()<<endl;
  cout<<"D_SL_p: "<<D_SL_p.ave_err()<<endl;
  cout<<"D_A_p_revins: "<<D_A_p_revins.ave_err()<<endl;
  cout<<"D_SL_p_revins: "<<D_SL_p_revins.ave_err()<<endl;

  raw_file_t D_obs_file("D_obs","w");
  D_obs_file.bin_write(D_M);
  D_obs_file.bin_write(D_SL_exch);
  D_obs_file.bin_write(D_SL_selftad);
  D_obs_file.bin_write(D_SL_s);
  D_obs_file.bin_write(D_SL_p);
  D_obs_file.bin_write(D_SL_selftad_revins);
  D_obs_file.bin_write(D_SL_p_revins);

  //Ds meson

  djack_t Ds_M,Ds_Z;

  djvec_t Ds_corr_00=load_P5P5("00",2,1);
  two_pts_true_fit(Ds_Z,Ds_M,Ds_corr_00,TH,Ds_tmin,Ds_tmax,"Ds_mass.xmg");
  cout<<"Ds_Z: "<<Ds_Z.ave_err()<<endl;
  cout<<"Ds_M: "<<Ds_M.ave_err()<<endl;

  raw_file_t Ds_obs_file("Ds_obs","w");
  Ds_obs_file.bin_write(Ds_M);

  ///////////////////////Delta m critico/////////////////////

  //ud
  
  djvec_t ud_contr_ward_cr=djvec_t(load_V0P5("LL",0,0)+2.0*djvec_t(load_V0P5("0M",0,0)+load_V0P5("0T",0,0)));
  djvec_t ud_contr_ward_cr_wo=djvec_t(load_V0P5("LL",0,0)+2.0*djvec_t(load_V0P5("0M",0,0)));
  djvec_t ud_num_deltam_cr=forward_derivative(ud_contr_ward_cr);
  djvec_t ud_num_deltam_cr_wo=forward_derivative(ud_contr_ward_cr_wo);
  ud_num_deltam_cr.ave_err().write("ud_num_deltam_cr.xmg");
  djvec_t ud_denom_deltam_cr=forward_derivative(load_V0P5("0P",0,0,RE,1));
  ud_denom_deltam_cr.ave_err().write("ud_denom_deltam_cr.xmg");
  
  djvec_t ud_deltam_cr=djvec_t(-ud_num_deltam_cr/ud_denom_deltam_cr)/2.0;
  ud_deltam_cr.ave_err().write("ud_deltam_cr_t.xmg");
  djvec_t ud_deltam_cr_wo=djvec_t(-ud_num_deltam_cr_wo/ud_denom_deltam_cr)/2.0;
  ud_deltam_cr.ave_err().write("ud_deltam_cr_t_wo.xmg");

  djack_t ud_fit_deltam_cr=constant_fit(ud_deltam_cr,ud_tmin,ud_tmax,"ud_fit_deltam_cr.xmg");
  djack_t ud_fit_deltam_cr_wo=constant_fit(ud_deltam_cr_wo,ud_tmin,ud_tmax,"ud_fit_deltam_cr_wo.xmg");
  
  raw_file_t ud_out_fit_deltam_cr("ud_fit_deltam_cr","w");
  ud_out_fit_deltam_cr.bin_write(ud_fit_deltam_cr);
  raw_file_t ud_out_fit_deltam_cr_wo("ud_fit_deltam_cr_wo","w");
  ud_out_fit_deltam_cr_wo.bin_write(ud_fit_deltam_cr_wo);

  cout<<"ud_deltam_cr:"<<" "<<ud_fit_deltam_cr.ave_err()<<endl;
  cout<<"ud_deltam_cr_wo:"<<" "<<ud_fit_deltam_cr_wo.ave_err()<<endl;

  //s

  djvec_t s_contr_ward_cr=djvec_t(load_V0P5("LL",1,1)+2.0*djvec_t(load_V0P5("0M",1,1)+load_V0P5("0T",1,1)));
  djvec_t s_num_deltam_cr=forward_derivative(s_contr_ward_cr);
  s_num_deltam_cr.ave_err().write("s_num_deltam_cr.xmg");
  djvec_t s_denom_deltam_cr=forward_derivative(load_V0P5("0P",1,1,RE,1));
  s_denom_deltam_cr.ave_err().write("s_denom_deltam_cr.xmg");
  
  djvec_t s_deltam_cr=djvec_t(-s_num_deltam_cr/s_denom_deltam_cr)/2.0;
  s_deltam_cr.ave_err().write("s_deltam_cr_t.xmg");

  djack_t s_fit_deltam_cr=constant_fit(s_deltam_cr,s_tmin,s_tmax,"s_fit_deltam_cr.xmg");
  
  raw_file_t s_out_fit_deltam_cr("s_fit_deltam_cr","w");
  s_out_fit_deltam_cr.bin_write(s_fit_deltam_cr);

  cout<<"s_deltam_cr:"<<" "<<s_fit_deltam_cr.ave_err()<<endl;

  //c

  djvec_t c_contr_ward_cr=djvec_t(load_V0P5("LL",2,2)+2.0*djvec_t(load_V0P5("0M",2,2)+load_V0P5("0T",2,2)));
  djvec_t c_num_deltam_cr=forward_derivative(c_contr_ward_cr);
  c_num_deltam_cr.ave_err().write("c_num_deltam_cr.xmg");
  djvec_t c_denom_deltam_cr=forward_derivative(load_V0P5("0P",2,2,RE,1));
  c_denom_deltam_cr.ave_err().write("c_denom_deltam_cr.xmg");
  
  djvec_t c_deltam_cr=djvec_t(-c_num_deltam_cr/c_denom_deltam_cr)/2.0;
  c_deltam_cr.ave_err().write("c_deltam_cr_t.xmg");

  djack_t c_fit_deltam_cr=constant_fit(c_deltam_cr,c_tmin,c_tmax,"c_fit_deltam_cr.xmg");
  
  raw_file_t c_out_fit_deltam_cr("c_fit_deltam_cr","w");
  c_out_fit_deltam_cr.bin_write(c_fit_deltam_cr);

  cout<<"c_deltam_cr:"<<" "<<c_fit_deltam_cr.ave_err()<<endl;  

  return 0;

}
