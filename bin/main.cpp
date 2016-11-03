#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <boot.hpp>
#include <oper.hpp>

#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  boot_t<double> a,b,c;
  cout<<a.size()<<endl;
  
  auto d=b+(a+c);
  
  return 0;
}
