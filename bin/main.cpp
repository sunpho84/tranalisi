#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <valarray>

using namespace std;

const int T=48,TH=24;

int main(int narg,char **arg)
{
  set_njacks(15);
  djvec_t corr_00=read_djvec("corr00_P5P5",T,0).symmetrized(1);
  djvec_t effm_00=effective_mass(corr_00,TH);
  
  constant_fit(effm_00,10,23,"test_00.xmg");

  djack_t Z,M;
  two_pts_fit(Z,M,corr_00,TH,10,24,"mass.xmg","Z.xmg");
  cout<<"Z: "<<Z.ave_err()<<endl;
  cout<<"M: "<<M.ave_err()<<endl;

  two_pts_migrad_fit(Z,M,corr_00,TH,10,24,"mass.xmg");
  cout<<"Z: "<<Z.ave_err()<<endl;
  cout<<"M: "<<M.ave_err()<<endl;
  
  ///////////////////// test slope //////////////////////////////
  
  djvec_t corr_0S=read_djvec("corr0S_P5P5",T,0).symmetrized(1);
  
  djack_t A,SL;
  two_pts_with_ins_ratio_fit(M,A,SL,corr_0S,corr_00,TH,10,24,"test_0S.xmg");
  cout<<"M: "<<M.ave_err()<<endl;
  cout<<"A: "<<A.ave_err()<<endl;
  cout<<"SL: "<<SL.ave_err()<<endl;
  
  return 0;
}
