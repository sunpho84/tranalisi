#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <boot.hpp>
#include <oper.hpp>

#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  dbvec_t a(3),b(3),c(3);
  a=3;
  b=1;
  c=4;

  dbvec_t d;

  d=a;
  d=-pow(log(b/(a+c)),2);
  //has_iterator<boot_t<double>>::value;

  cout<<d[0][1]<<endl;
  
  
  //cout<<has_method_size<boot_t<double>>::value<<endl;
  
  //decltype(boot_t<double>::iterator(),0);
  
  return 0;
}
