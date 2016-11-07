#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <boot.hpp>
#include <oper.hpp>

#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  set_njacks(2);
  
  dbvec_t a(3),b(3);
  dboot_t c;
  a=3;
  b=1;
  c=4;

  dbvec_t d;

  uminus(1);
  
  d=-pow(log(b/(a+c)),2)+2;
  //has_iterator<boot_t<double>>::value;

  cout<<d[0][1]<<endl;

  ave_err_t av=d[0].ave_err();
  cout<<av<<endl;
  
  //cout<<has_method_size<boot_t<double>>::value<<endl;
  
  //decltype(boot_t<double>::iterator(),0);
  
  dboot_t test1(gauss_filler_t(3,2,124142));
  dboot_t test2(gauss_filler_t(3,2,124142));
  cout<<(test1/test2).ave_err()<<endl;
  cout<<test2.ave_err()<<endl;

  boot_init_t binit;binit.fill(23252523);
  vector<vector<double>> ciccio={{3,4},{4,5}};
  
  dbvec_t bciccio(binit,ciccio);
  
  cout<<bciccio.ave_err()<<endl;
  
  return 0;
}
