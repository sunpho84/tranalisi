#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <array>
#include <map>
#include <vector>

#include <Dirac.hpp>
#include <meas_vec.hpp>

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
  
  //! return a single component of phat
  double hat_mu(size_t mu) const
  {return 2.0*sin((*this)[mu]/2.0);}
  
  //! return the ptilde
  p_t tilde() const
  {
    p_t out;
    for(size_t mu=0;mu<NDIM;mu++) out[mu]=tilde_mu(mu);
    return out;
  }
  
  //! return the phat
  p_t hat() const
  {
    p_t out;
    for(size_t mu=0;mu<NDIM;mu++) out[mu]=hat_mu(mu);
    return out;
  }
  
  //! return \sum_mu(p_mu^x)
  template <class T>
  double normx(T x) const
  {
    double out=0;
    for(size_t mu=0;mu<NDIM;mu++) out+=pow((*this)[mu],x);
    return out;
  }
  
  //! return p2
  double norm2() const
  {return normx(2);}
  
  //! return p4
  double norm4() const
  {return normx(4);}
  
  //! return \sum_nu(p_nu^4)/\sum_nu(p_nu^2)^2
  double p4_fr_p22()
  {return norm4()/sqr(norm2());}
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
  
  //! return the number of components different from 0
  size_t Np() const
  {
    size_t out=0;
    for(size_t mu=0;mu<NDIM;mu++) out+=((*this)[mu]!=0 or fabs(ph_mom[mu])>1.0e-10);
    return out;
  }
  
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
EXTERN_GEOMETRY vector<bool> filt_moms; //!< store if a momentum passed the filter

//! reads the list of momenta from file
void set_list_of_moms(const string &path,double thresh=0.28);

//! list of class of equivalence of momenta
EXTERN_GEOMETRY vector<pair<size_t,vector<size_t>>> equiv_imoms;

//! list of class of equivalence of momenta
EXTERN_GEOMETRY vector<size_t> iequiv_mom_of_ifilt;

//! fills the list of class of equivalence
void set_class_of_equiv_moms();

//! set a list of imom from filter moms
void set_filtered_moms(const double thresh=0.28);

//! get a vector of pt2 for filtered moms
vector<double> get_filtered_pt2();

//! get a vector of pt2 for indep moms
vector<double> get_indep_pt2();

//! get a vector of pt2 for all moms
vector<double> get_pt2();

//! list all smom pairs
void list_all_smom_pairs();

//! average all equivalent momenta, using an index to access moments (assumed to be last component)
djvec_t average_equiv_moms(const djvec_t &in,const index_t &oth_ind_mom_ind,const index_t &oth_mom_ind);

//! take a vector of indep moms and returns the list of filtered moms
djvec_t get_filtered_moms(const djvec_t &in);

//! return the mirrorized site
size_t get_mir_mom(size_t imom,size_t imir);

#undef EXTERN_GEOMETRY
#undef INIT_TO

#endif
