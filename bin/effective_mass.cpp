#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  if(narg<4) CRASH("Use %s file T njacks iel=0");
  
  //setup the number of jackknives from command line
  set_njacks(atoi(arg[3]));
  
  int T=atoi(arg[2]);
  
  djvec_t data(T);
  
  data.bin_read(arg[1]);
  
  return 0;
}
