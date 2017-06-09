#ifndef _DIRAC_HPP
#define _DIRAC_HPP

#include <index.hpp>
#include <jack.hpp>

#include <compl.hpp>
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

using prop_t=Matrix<dcompl_t,NSPINCOL,NSPINCOL>; //< Matrix of complex
using jprop_t=Matrix<cdjack_t,NSPINCOL,NSPINCOL>; //< Matrix of djack_t

namespace Eigen
{
  //! Traits for using cdjack_t as a scalar type for Eigen
  template<> struct NumTraits<cdjack_t>
    : NumTraits<double>
  {
    typedef djack_t Real;
    typedef cdjack_t NonInteger;
    typedef cdjack_t Nested;
    enum
      {
	IsComplex=1,
	IsInteger=0,
	IsSigned=1,
	RequireInitialization=1,
	ReadCost=30,
	AddCost=120,
	MulCost=120
      };
  };
}

//! clusterize a jprop_t
void clusterize(jprop_t &j,size_t clust_size=1);

//! put into a cluster
void put_into_cluster(jprop_t &jprop,const prop_t &prop,size_t iclust);

//! put into a jackknife
void put_into_jackknife(jprop_t &jprop,const prop_t &prop,size_t ijack);

//! get from a jackkinfe
prop_t get_from_jackknife(const jprop_t &jprop,size_t ijack);

//! number of matrices in Clifford basis
const size_t nGamma=16;

//! sparse matrix for Dirac
using Dirac_t=SparseMatrix<dcompl_t>;

//! sixteen Clifford basis
extern vector<Dirac_t> Gamma,vit_Gamma;

EXTERN_DIRAC const size_t igmu[4] INIT_TO({4,1,2,3});

//! get a single Gamma
Dirac_t init_Gamma(const int *irow,const int *re,const int *im);

//! return the colorspin index
inline size_t isc(size_t is,size_t ic)
{return ic+NCOL*is;}

//! convert to Vittorio's basis
prop_t convert_to_Vit_basis(const prop_t &p);

//! convert from Vittorio's basis
prop_t convert_from_Vit_basis(const prop_t &p);

#undef EXTERN_DIRAC
#undef INIT_TO

#endif
