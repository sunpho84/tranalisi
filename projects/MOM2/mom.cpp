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
  scan::scan("analysis.txt");
  
  //print time statistics
  cout<<ts<<endl;
  
  return 0;
}
