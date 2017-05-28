#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t NCOL=3;
const size_t NSPIN=4;
const int NSPINCOL=NSPIN*NCOL;

using jprop_t=Matrix<djack_t,NSPINCOL,NSPINCOL>;
using dirac_t=Matrix<dcomplex,NSPIN,NSPIN>;

namespace Eigen
{
  template<> struct NumTraits<djack_t>
    : NumTraits<double> 
  {
    typedef djack_t Real;
    typedef djack_t NonInteger;
    typedef djack_t Nested;
    enum {
      IsComplex=0,
      IsInteger=0,
      IsSigned=1,
      RequireInitialization=1,
      ReadCost=1,
      AddCost=3,
      MulCost=3
    };
  };
}

//! specify hot to print a djack_t
ostream& operator<<(ostream &out,const djack_t &v)
{return out<<v.ave_err();}




int main()
{
  set_njacks(15);
  
  jprop_t a,b;
  
  a.Zero();
  b=b.Identity();
  a+=b;
  
  a*=b;
  
  cout<<b<<endl;
  
  return 0;
}
