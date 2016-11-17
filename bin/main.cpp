#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <valarray>

using namespace std;

int main(int narg,char **arg)
{
  set_njacks(15);
  djvec_t data=effective_mass(read_djvec("corr00_P5P5",48,0).symmetrized(0),1);
  
  constant_fit(data,10,23,"test.xmg");

  valarray<djack_t> f;
  
  return 0;
}
