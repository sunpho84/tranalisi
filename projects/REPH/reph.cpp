#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  set_njacks(15);
  
  const int T=48;
  
  djvec_t corr(T);
  
  corr.bin_read("jacks/oPPo-ss");
  
  cout<<corr.ave_err()<<endl;
  
  return 0;
}
