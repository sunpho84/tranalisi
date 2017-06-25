#ifndef _PROP_HPP
#define _PROP_HPP

#include <geometry.hpp>

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using vprop_t=vector<prop_t>; //! vector of propagators
using vjprop_t=vector<jprop_t>; //! vector of jackkniffed props

EXTERN_PROP size_t nm; //!< number of masses
EXTERN_PROP size_t nr; //!< number of r
EXTERN_PROP size_t nmr; //total number of m and r

//! add the prop of a given conf on the jackknife
void build_jackknifed_prop(vjprop_t &jprop,const vprop_t &prop,size_t ijack);

//! index of the (im,ir) prop
EXTERN_PROP index_t mr_ind;
void set_mr_ind(size_t nm,size_t nr);

//! read a propagator
void read_prop(prop_t &prop,raw_file_t &file,const dcompl_t &fact);

//! holds a given l r and momentum, for a fized conf
class m_r_mom_conf_props_t
{
public:
  prop_t prop_0; //!< propagator with no insertion
  
  prop_t prop_FF; //!< propagator with 2 photon insertions
  prop_t prop_F; //!< propagator with 1 photon insertion
  prop_t prop_T; //!< propagator with Tadpole insertion
  prop_t prop_S; //!< propagator with Scalar insertion
  prop_t prop_P; //!< propagator with Pseudoscalar insertion
};

//! return the tag of a prop
string get_prop_tag(size_t im,size_t ir,const string &ins);

//! holds a given m and r, all jackks and moms
class jm_r_mom_props_t
{
public:
  vjprop_t jprop_0; //!< propagator with no photon, given momentum, all m and r
  
  vjprop_t jprop_2; //!< propagator with 2 photons, Tadpole, Pseudoscalar and possibly Scalar insertions
  vjprop_t jprop_P; //!< propagator with Pseudoscalar inserion
  vjprop_t jprop_S; //!< propagator with Scalar insertion
};
  
EXTERN_PROP index_t m_r_ind; //!| index of m,r
EXTERN_PROP index_t conf_hit_ind; //!< index of conf and hit

//! set all jackknifed props
void set_jprops(bool set_QED);

//! clusterize the propagators
void clusterize_all_mr_props(bool use_QED,size_t clust_size);

//! compute the inverse of all mr props
vjprop_t get_all_mr_props_inv(const vjprop_t &jprop);

#undef EXTERN_PROP
#undef INIT_TO

#endif
