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

using lprop_t=Matrix<dcompl_t,NSPIN,NSPIN>; //< Matrix of complex with spin
using qprop_t=Matrix<dcompl_t,NSPINCOL,NSPINCOL>; //< Matrix of complex, with spin and color
using jqprop_t=jack_t<qprop_t>;

//! number of matrices in Clifford basis
const size_t nGamma=16;

//! sparse matrix for Dirac
using Dirac_t=SparseMatrix<dcompl_t>;

//! sixteen Clifford basis
extern vector<Dirac_t> Gamma,vit_Gamma;

const size_t igmu[4]={4,1,2,3};

//! get a single Gamma
Dirac_t init_Gamma(const int *irow,const int *re,const int *im);

//! return the spincolor index
inline size_t isc(size_t is,size_t ic)
{return ic+NCOL*is;}

//! convert to Vittorio's basis
qprop_t convert_to_Vit_basis(const qprop_t &p);

//! convert from Vittorio's basis
qprop_t convert_from_Vit_basis(const qprop_t &p);

//! invert the jprop
jqprop_t invert(const jqprop_t &in);

#undef EXTERN_DIRAC
#undef INIT_TO

#endif
