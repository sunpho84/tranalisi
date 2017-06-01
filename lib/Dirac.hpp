#ifndef _DIRAC_HPP
#define _DIRAC_HPP

#include <jack.hpp>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

//! sparse matrix for Dirac
using Dirac_t=SparseMatrix<dcomplex>;

//! sixteen Clifford basis
extern vector<Dirac_t> Gamma;

//! get a single Gamma
Dirac_t init_Gamma(const int *irow,const int *re,const int *im);

#undef EXTERN_DIRAC
#undef INIT_TO

#endif
