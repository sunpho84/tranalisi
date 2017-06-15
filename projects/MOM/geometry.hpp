#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <array>
#include <map>
#include <vector>

#include <Dirac.hpp>

#ifndef EXTERN_GEOMETRY
 #define EXTERN_GEOMETRY extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using namespace std;

//! number of dimensions
const int NDIM=4;

//! holds a coordinate
using coords_t=array<int,NDIM>;

//! holds a momentum
class p_t : public array<double,NDIM>
{
public:
  //! return a single component of ptilde
  double tilde_mu(size_t mu) const
  {return sin((*this)[mu]);}
  
  //! return the ptilde
  p_t tilde()
  {
    p_t out;
    for(size_t mu=0;mu<NDIM;mu++) out[mu]=tilde_mu(mu);
    return out;
  }
  
  //! return \sum_mu(p_mu^x)
  template <class T>
  double normx(T x)
  {
    double out=0;
    for(size_t mu=0;mu<NDIM;mu++) out+=pow((*this)[mu],x);
    return out;
  }
  
  //! return p2
  double norm2() {return normx(2);}
  
  //! return p4
  double norm4() {return normx(4);}
  
  //! return /\sum_nu(p_nu^2)
  double p4_fr_p2()
  {return norm4()/norm2();}
};

//! return pslash
inline prop_t slash(const p_t &p)
{return Gamma[4]*p[0]+Gamma[1]*p[1]+Gamma[2]*p[2]+Gamma[3]*p[3];}

//! phase of momenta
extern double ph_mom[NDIM];

//! holds the components of a momentum
class imom_t : public coords_t
{
public:
  //! return a single component of p
  double p_mu(size_t mu,const coords_t &L) const
  {return 2*M_PI*((*this)[mu]+ph_mom[mu])/L[mu];}
  
  //! return the p
  p_t p(const coords_t &L) const
  {
    p_t out;
    for(size_t mu=0;mu<NDIM;mu++) out[mu]=p_mu(mu,L);
    return out;
  }
  
  //! return the quadrant
  size_t imir()
  {
    size_t out=0;
    for(size_t mu=0;mu<NDIM;mu++) out+=((*this)[mu])<<mu;
    return out;
  }
};

EXTERN_GEOMETRY coords_t L; //!< lattice sizes
EXTERN_GEOMETRY size_t V; //!< lattice volume
EXTERN_GEOMETRY vector<imom_t> imoms; //!< list of momenta
EXTERN_GEOMETRY vector<size_t> Np_class; //< sum of Np over the class

//! reads the list of momenta from file
void get_list_of_moms(const string &path);

//! list of class of equivalence of momenta
EXTERN_GEOMETRY vector<pair<size_t,vector<size_t>>> equiv_imoms;

//! fills the list of class of equivalence
void get_class_of_equiv_moms();

//! get a vector of pt2 for indep moms
vector<double> get_indep_pt2();

//! list all smom pairs
void list_all_smom_pairs();

#undef EXTERN_GEOMETRY
#undef INIT_TO

#endif
