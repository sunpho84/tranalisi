#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <valarray>

using namespace std;

int main(int narg,char **arg)
{
  set_njacks(15);
  djvec_t corr=read_djvec("corr00_P5P5",48,0).symmetrized(1);
  djvec_t effm=effective_mass(corr);
  
  constant_fit(effm,10,23,"test.xmg");

  djack_t Z,M;
  two_pts_fit(Z,M,corr, 10,24, 24,"mass.xmg","Z.xmg");
  cout<<"Z: "<<Z.ave_err()<<endl;
  cout<<"M: "<<M.ave_err()<<endl;

  two_pts_migrad_fit(Z,M,corr,10,24,24,"mass.xmg");
  cout<<"Z: "<<Z.ave_err()<<endl;
  cout<<"M: "<<M.ave_err()<<endl;
  
  return 0;
}
