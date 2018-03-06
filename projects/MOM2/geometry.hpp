#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <array>
#include <map>
#include <vector>

#include <Dirac.hpp>
#include <meas_vec.hpp>
#include <types.hpp>

#ifndef EXTERN_GEOMETRY
 #define EXTERN_GEOMETRY extern
 #define INIT_GEOMETRY_TO(...)
#else
 #define INIT_GEOMETRY_TO(...) = __VA_ARGS__
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
  {
    return sin((*this)[mu]);
  }
  
  //! return a single component of phat
  double hat_mu(size_t mu) const
  {
    return 2.0*sin((*this)[mu]/2.0);
  }
  
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
  {
    return normx(2);
  }
  
  //! return p4
  double norm4() const
  {
    return normx(4);
  }
  
  //! return \sum_nu(p_nu^4)/\sum_nu(p_nu^2)^2
  double p4_fr_p22()
  {
    return norm4()/sqr(norm2());
  }
  
  //! operator to print a imom_t
  friend ostream& operator<<(ostream &os,const p_t &m)
  {
    os<<"{"<<m[0];
    for(size_t mu=1;mu<NDIM;mu++) os<<","<<m[mu];
    os<<"}";
    
    return os;
  }
};

//! return pslash with a given set of gamma
inline auto slash(const vector<Dirac_t> &Gamma,const p_t &p)
{
  return Gamma[4]*p[0]+Gamma[1]*p[1]+Gamma[2]*p[2]+Gamma[3]*p[3];
}

//! return pslash for quark
inline qprop_t qua_slash(const p_t &p)
{
  return slash(quaGamma,p);
}

//! return pslash
inline lprop_t lep_slash(const p_t &p)
{
  return slash(lepGamma,p);
}

//! phase of momenta
extern double ph_mom[NDIM];

//! kind of boundary condition supported
namespace temporal_bc
{
  enum type_t{ANTIPERIODIC,PERIODIC};
  const map<string,tuple<type_t,double>> decr{{"Antiperiodic",{ANTIPERIODIC,0.5}},{"Periodic",{PERIODIC,0.0}}};
  EXTERN_GEOMETRY type_t bc;
  const size_t n=decr.size();
  PROVIDE_DECRYPTER;
}

//! holds the components of a momentum
class imom_t : public coords_t
{
public:
  //! return a single component of p
  double p_mu(size_t mu,const coords_t &L) const
  {
    return 2*M_PI*((*this)[mu]+ph_mom[mu])/L[mu];
  }
  
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
  
  //! operator to print a imom_t
  friend ostream& operator<<(ostream &os,const imom_t &m)
  {
    os<<"{"<<m[0];
    for(size_t mu=1;mu<NDIM;mu++) os<<","<<m[mu];
    os<<"}";
    
    return os;
  }
};

#endif
