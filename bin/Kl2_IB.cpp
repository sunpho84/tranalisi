#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use %s T",arg[0]);
  
  int T=atoi(arg[1]);
  

  
  set_njacks(15);
  
 

  index_t<4> ind;
  ind.set_ranges({3,3,2,2});
  
  //test
  djvec_t dpi000=read_djvec("corr00_P5P5",T,ind({0,0,0,0}));
  djvec_t dpi001=read_djvec("corr00_P5P5",T,ind({0,0,1,0}));
  djvec_t dpi0S0=read_djvec("corr0S_P5P5",T,ind({0,0,0,0}));

  djvec_t jratio=dpi0S0/dpi000;
  
  grace_file_t out_test("test_ratio.xmg");
  out_test<<jratio.ave_err();

  //pion

  djvec_t dpiLL0=read_djvec("corrLL_P5P5",T,ind({0,0,0,0}));
  djvec_t dpiLL1=read_djvec("corrLL_P5P5",T,ind({0,0,1,0}));

  djvec_t pi_ratio=(dpiLL0+dpiLL1)/(dpi000+dpi001);

  grace_file_t out_pi("pi_ratio.xmg");
  out_pi<<pi_ratio.ave_err();

  //kaon

  //0,1,0,0 è identico a 1,0,0,0
  djvec_t dk000=read_djvec("corr00_P5P5",T,ind({1,0,0,0}));
  djvec_t dk001=read_djvec("corr00_P5P5",T,ind({1,0,1,0}));

  //1,0,0,0 è identico a 0,1,0,0
  djvec_t dkLL0=read_djvec("corrLL_P5P5",T,ind({0,1,0,0}));
  djvec_t dkLL1=read_djvec("corrLL_P5P5",T,ind({0,1,1,0}));

  djvec_t dk0M0=read_djvec("corr0M_P5P5",T,ind({1,0,0,0}));
  djvec_t dk0M1=read_djvec("corr0M_P5P5",T,ind({1,0,1,0}));

  djvec_t dk0T0=read_djvec("corr0T_P5P5",T,ind({1,0,0,0}));
  djvec_t dk0T1=read_djvec("corr0T_P5P5",T,ind({1,0,1,0}));

  djvec_t dk0S0=read_djvec("corr0S_P5P5",T,ind({1,0,0,0}));
  djvec_t dk0S1=read_djvec("corr0S_P5P5",T,ind({1,0,1,0}));

  djvec_t dk0P0=read_djvec("corr0P_P5P5",T,ind({1,0,0,1}));
  djvec_t dk0P1=read_djvec("corr0P_P5P5",T,ind({1,0,1,1}));

  djvec_t k_ratio_exch=(dkLL0+dkLL1)/(dk000+dk001);
  djvec_t k_ratio_self=(dk0M0+dk0M1)/(dk000+dk001);
  djvec_t k_ratio_tad=(dk0T0+dk0T1)/(dk000+dk001);
  djvec_t k_ratio_s=(dk0S0+dk0S1)/(dk000+dk001);
  djvec_t k_ratio_p=(dk0P0+dk0P1)/(dk000+dk001);


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
  
  

  
  
  return 0;
}
