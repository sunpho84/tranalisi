#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=48;

djvec_t load(const string &name)
{
  djvec_t out=read_djvec("data/corr"+name+"_ss",T,0).symmetrized();
  out.ave_err().write("plots/"+name+".xmg");
  return out;
}

double err_rat(const djvec_t &full,const djvec_t &diag)
{
  double out=0;
  for(size_t iel=0;iel<full.size();iel++) out+=sqr(full[iel].err()/diag[iel].err());
  return sqrt(out/full.size());
}

int main()
{
  set_njacks(15);
  
  djvec_t LL_diag_P5P5=load("LL_diag_P5P5");
  djvec_t LL_P5P5=load("LL_P5P5");
  cout<<"LL_P5P5_err_ratio: "<<err_rat(LL_P5P5,LL_diag_P5P5)<<endl;
  
  djvec_t LL_diag_VKVK=load("LL_diag_VKVK");
  djvec_t LL_VKVK=load("LL_VKVK");
  cout<<"LL_VKVK_err_ratio: "<<err_rat(LL_VKVK,LL_diag_VKVK)<<endl;
  
  djvec_t OM_diag_P5P5=load("0M_diag_P5P5");
  djvec_t OM_P5P5=load("0M_P5P5");
  cout<<"0M_P5P5_err_ratio: "<<err_rat(OM_P5P5,OM_diag_P5P5)<<endl;
  
  djvec_t OM_diag_VKVK=load("0M_diag_VKVK");
  djvec_t OM_VKVK=load("0M_VKVK");
  cout<<"0M_VKVK_err_ratio: "<<err_rat(OM_VKVK,OM_diag_VKVK)<<endl;
  
  //djvec_t LL_offdiag_VKVK=LL_VKVK-LL_diag_VKVK;
  //LL_offdiag_VKVK.ave_err().write("plots/LL_offdiag_VKVK.xmg");
  
  return 0;
}
