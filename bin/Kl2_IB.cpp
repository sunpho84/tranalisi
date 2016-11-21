#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int T;
index_t<4> ind;

djvec_t load(const char *path,size_t im1,size_t im2,size_t reim=RE,int rpar=1,int spat_par=1)
{
  djvec_t corr_r0=read_djvec(path,T,ind({im1,im2,0,reim}));
  djvec_t corr_r1=read_djvec(path,T,ind({im1,im2,1,reim}));

  return djvec_t(corr_r0+rpar*corr_r1).symmetrized(spat_par)*0.5;
}
  
int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use %s T",arg[0]);
  
  T=atoi(arg[1]);
  int TH=T/2;
  
  ind.set_ranges({3,3,2,2});  
  set_njacks(15);
  
  //pion
  
  djvec_t pi_ratio=load("corrLL_P5P5",0,0)/load("corr00_P5P5",0,0);

  grace_file_t out_pi("pi_ratio.xmg");
  out_pi<<pi_ratio.ave_err();

  
  //kaon


  djvec_t k_ratio_exch=load("corrLL_P5P5",0,1)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_self=load("corr0M_P5P5",1,0)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_tad=load("corr0T_P5P5",1,0)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_s=load("corr0S_P5P5",1,0)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_p=load("corr0P_P5P5",1,0,IM,-1)/load("corr00_P5P5",1,0);


  grace_file_t out_k_exch("k_ratio_exch.xmg");
  grace_file_t out_k_self("k_ratio_self.xmg");
  grace_file_t out_k_tad("k_ratio_tad.xmg");
  grace_file_t out_k_s("k_ratio_s.xmg");
  grace_file_t out_k_p("k_ratio_p.xmg");
  out_k_exch<<k_ratio_exch.ave_err();
  out_k_self<<k_ratio_self.ave_err();
  out_k_tad<<k_ratio_tad.ave_err();
  out_k_s<<k_ratio_s.ave_err();
  out_k_p<<k_ratio_p.ave_err();

  //D meson

 
  djvec_t D_ratio_exch=load("corrLL_P5P5",0,2)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_self=load("corr0M_P5P5",2,0)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_tad=load("corr0T_P5P5",2,0)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_s=load("corr0S_P5P5",2,0)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_p=load("corr0P_P5P5",2,0,IM,-1)/load("corr00_P5P5",2,0);


  grace_file_t out_D_exch("D_ratio_exch.xmg");
  grace_file_t out_D_self("D_ratio_self.xmg");
  grace_file_t out_D_tad("D_ratio_tad.xmg");
  grace_file_t out_D_s("D_ratio_s.xmg");
  grace_file_t out_D_p("D_ratio_p.xmg");
  out_D_exch<<D_ratio_exch.ave_err();
  out_D_self<<D_ratio_self.ave_err();
  out_D_tad<<D_ratio_tad.ave_err();
  out_D_s<<D_ratio_s.ave_err();
  out_D_p<<D_ratio_p.ave_err();

  djack_t M,A,SL;
  two_pts_with_ins_ratio_fit(M,A,SL,load("corrLL_P5P5",0,0),load("corr00_P5P5",0,0),TH,10,24,"test.xmg","test_ins.xmg");
  cout<<"M: "<<M.ave_err()<<endl;
  cout<<"A: "<<A.ave_err()<<endl;
  cout<<"SL: "<<SL.ave_err()<<endl;
  


  
  
  return 0;
}
