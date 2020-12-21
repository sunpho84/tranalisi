#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <tranalisi.hpp>

#include "node.hpp"

#include "fitter/fitter_parser.hpp"
#include "operations.hpp"

using namespace std;

vector<double> pars;

void init_scanner();
void destroy_scanner();
// int main()
// {
  
//   //   real val(10);
//   // uexp e(exp,&val);
//   // bexp s(sum,&val,&e);
//   // uexp ss(sin,&e);
//   // cout<<ss.eval()<<endl;
  
//   return 0;
// }

int i=0;
string str("8*(8--7)");

void fitter_lex_internal(char* buf,int &result,int max_size)
{
  result=(str[i]!='\0');
  *buf=str[i++];
}

int main(int narg,char **arg)
{
  if(narg<2)
    CRASH("Use %s formula",arg[0]);
  
  for(int iarg=1;iarg<narg;iarg++)
    {
      str=arg[iarg];
      i=0;
      
      init_scanner();
      
      fitter_parser_parse(nullptr);
      
      destroy_scanner();
    }
  
  return 0;
}

