#include <tranalisi.hpp>

#include <phiFunction.hpp>

#include <iostream>

using namespace std;

int main()
{
  const int nIntervals=200;
  
  PhiCalculator phiCalculator(0.,15,nIntervals);
  
  grace_file_t out("/tmp/test.xmg");
  for(double x=0;x<15;x+=0.1)
    out.write_xy(x,phiCalculator(x));
  
  return 0;
}
