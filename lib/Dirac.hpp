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
extern vector<Dirac_t> quaGamma,lepGamma;

const size_t igmu[4]={4,1,2,3};

//! get a single Gamma
template <int Ncol>
Dirac_t init_Gamma(const int *irow,const int *re,const int *im)
{
  //! matrix being built
  Dirac_t m(NSPIN*Ncol,NSPIN*Ncol);
  
  //! list of entries
  vector<Triplet<dcompl_t>> tr;
  
  for(size_t ispin=0;ispin<NSPIN;ispin++)
    for(size_t icol=0;icol<Ncol;icol++)
      tr.push_back(Triplet<dcompl_t>(ispin*Ncol+icol,irow[ispin]*Ncol+icol,dcompl_t(re[ispin],im[ispin])));
  
  m.setFromTriplets(tr.begin(),tr.end());
  
  return m;
}

//! type of data to initialize Dirac matrix
typedef int Gamma_data_t[3][NSPIN];

//! initialize the gamma basis
template <int Ncol>
vector<Dirac_t> init_all_Gamma(const Gamma_data_t *data)
{
  //! returned list of Gamma
  vector<Dirac_t> out;
  
  for(size_t iGamma=0;iGamma<nGamma;iGamma++)
    out.push_back(init_Gamma<Ncol>(data[iGamma][0],data[iGamma][1],data[iGamma][2]));
  
  return out;
}

//! return the spincolor index
inline size_t isc(size_t is,size_t ic)
{return ic+NCOL*is;}

//! invert the jprop
jqprop_t invert(const jqprop_t &in);

#undef EXTERN_DIRAC
#undef INIT_TO

#endif
