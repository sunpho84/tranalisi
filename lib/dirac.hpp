#ifndef _DIRAC_HPP
#define _DIRAC_HPP

#include <jack.hpp>
#include <vector>
#include <eigen3/Eigen/Dense>

#ifndef EXTERN_DIRAC
 #define EXTERN_DIRAC extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using namespace Eigen;
using std::array;

const size_t NCOL=3; //< Number of colors
const size_t NSPIN=4; //< Number of spin components
const int NSPINCOL=NSPIN*NCOL; //< Total number of spin and color components

using jprop_t=Matrix<djack_t,NSPINCOL,NSPINCOL>; //< Matrix of djack_t

//! Dirac matrix
class Dirac_t : public Matrix<dcomplex,NSPIN,NSPIN>
{
public:
  Dirac_t(const vector<int> &pars) //: Matrix<dcomplex,NSPIN,NSPIN>()
  {for(size_t i=0;i<NSPIN;i++) (*this)(i,pars[i])=dcomplex(pars[NSPIN+i],pars[2*NSPIN+i]);}
};

namespace Eigen
{
  //! Traits for using djack_t as a scalar type for Eigen
  template<> struct NumTraits<djack_t>
    : NumTraits<double>
  {
    typedef djack_t Real;
    typedef djack_t NonInteger;
    typedef djack_t Nested;
    enum
      {
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

//! number of matrices in Clifford basis
const size_t nGamma=16;

//! sixteen Clifford basis
EXTERN_DIRAC vector<Dirac_t> Gamma INIT_TO
({
  vector<int>({{0,1,2,3 , 1,1,1,1 , 0,0,0,0}}),
  vector<int>({{3,2,1,0 , 0,0,0,0 , -1,-1,1,1}}),
  vector<int>({{3,2,1,0 , -1,1,1,-1 , 0,0,0,0}}),
  vector<int>({{2,3,0,1 , 0,0,0,0 , -1,1,1,-1}}),
  vector<int>({{2,3,0,1 , -1,-1,-1,-1 , 0,0,0,0}}),
  vector<int>({{0,1,2,3 , 1,1,-1,-1 , 0,0,0,0}}),
  vector<int>({{3,2,1,0 , 0,0,0,0 , 1,1,1,1}}),
  vector<int>({{3,2,1,0 , 1,-1,1,-1 , 0,0,0,0}}),
  vector<int>({{2,3,0,1 , 0,0,0,0 , 1,-1,1,-1}}),
  vector<int>({{2,3,0,1 , 1,1,-1,-1 , 0,0,0,0}}),
  vector<int>({{1,0,3,2 , 0,0,0,0 , -1,-1,1,1}}),
  vector<int>({{1,0,3,2 , -1,1,1,-1 , 0,0,0,0}}),
  vector<int>({{0,1,2,3 , 0,0,0,0 , -1,1,1,-1}}),
  vector<int>({{1,0,3,2 , 0,0,0,0 , 1,1,1,1}}),
  vector<int>({{1,0,3,2 , 1,-1,1,-1 , 0,0,0,0}}),
  vector<int>({{0,1,2,3 , 0,0,0,0 , 1,-1,1,-1}})});

#undef EXTERN_DIRAC
#undef INIT_TO

#endif
