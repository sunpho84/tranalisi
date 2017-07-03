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

//! read a propagator
void read_prop(prop_t *prop,raw_file_t &file,const dcompl_t &fact);

//! holds a given l r and momentum, for a fized conf
class m_r_mom_conf_props_t
{
public:
  prop_t LO; //!< propagator with no insertion
  
  prop_t FF; //!< propagator with 2 photon insertions
  prop_t F; //!< propagator with 1 photon insertion
  prop_t T; //!< propagator with Tadpole insertion
  prop_t S; //!< propagator with Scalar insertion
  prop_t P; //!< propagator with Pseudoscalar insertion
};

//! return the tag of a prop
string get_prop_tag(size_t im,size_t ir,const string &ins);

//! holds a given m and r, all jackks
class jm_r_mom_props_t
{
  //! return a list of pointers to all internal data
  vector<jprop_t*> get_all_ptrs(bool use_QED)
  {
    vector<jprop_t*> out={&LO};
    if(use_QED)
      for(auto &p : {&EM,&P,&S})
	out.push_back(p);
    return out;
  }
  
public:
  jprop_t LO; //!< propagator with no photon, given momentum, all m and r
  
  jprop_t EM; //!< propagator with 2 photons, Tadpole, Pseudoscalar and possibly Scalar insertions
  jprop_t P; //!< propagator with Pseudoscalar inserion
  jprop_t S; //!< propagator with Scalar insertion
  
  //! clusterize the propagators
  void clusterize_all_mr_props(bool use_QED,size_t clust_size)
  {
    //cout<<"Clusterizing all props, clust_size="<<clust_size<<endl;
    for(auto &p : get_all_ptrs(use_QED))
      clusterize(*p,clust_size);
  }
};

//! add the prop of a given conf on the jackknife
void build_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,const vector<m_r_mom_conf_props_t> &props,bool set_QED,const index_t &im_r_ijack_ind,const index_t &im_r_ind);

//! compute the inverse of all mr props
vjprop_t get_all_mr_props_inv(const vjprop_t &jprop);

//! clusterize all props
void clusterize_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,bool use_QED,size_t clust_size);

//! return the propagator at maximal twist
prop_t free_prop(const imom_t &pi,double mu,double kappa,size_t r);

//! finish computing em prop
void finish_jprops_EM(vector<jm_r_mom_props_t> &jprops,const djack_t &deltam_cr);

#undef EXTERN_PROP
#undef INIT_TO

#endif
