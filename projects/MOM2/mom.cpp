#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <scan.hpp>
#include <timings.hpp>

int main(int narg,char **arg)
{
  string path="analysis.txt";
  if(narg>1)
    path=arg[1];
  
  scan::scan(path);
  
  //print time statistics
  cout<<ts<<endl;
  
  return 0;
}
