#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  if(narg<4) CRASH("Use %s file T njacks iel=0 par=1",arg[0]);
  
  //parse all pars
  const char *path=arg[1];
  int T=atoi(arg[2]);
  int ext_njacks=atoi(arg[3]);
  int iel=0;
  if(narg>=4) iel=atoi(arg[4]);
  int par=1;
  if(narg>=5) par=atoi(arg[5]);
  
  //setup the number of jackknives from command line
  set_njacks(ext_njacks);
  //read
  auto data=read_djvec(path,T,iel);
  //write average and error
  cout<<effective_mass(data.symmetrized(par)).ave_err()<<endl;
  
  return 0;
}
