#ifndef _PROP_HPP
#define _PROP_HPP

#include <geometry.hpp>

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using vprop_t=vector<prop_t>; //! vector of propagators (spanning moms)
using vjprop_t=vector<jprop_t>; //! vector of jackkniffed props (spanning moms)

EXTERN_PROP size_t nm; //!< number of masses
EXTERN_PROP size_t nr; //!< number of r
EXTERN_PROP size_t nmr; //total number of m and r

//! read a propagator from a given file
vprop_t read_prop(const string &template_path,size_t iconf,size_t ihit);

//! add the prop of a given conf on the jackknife
void build_jackknifed_prop(vjprop_t &jprop,const vprop_t &prop,size_t ijack);

//! index of the (im,ir) prop
EXTERN_PROP index_t mr_ind;
void set_mr_ind(size_t nm,size_t nr);

EXTERN_PROP vector<vprop_t> conf_prop_0; //!< propagator with no insertion for a given conf

EXTERN_PROP vector<vprop_t> conf_prop_FF; //!< propagator with 2 photon insertions
EXTERN_PROP vector<vprop_t> conf_prop_F; //!< propagator with 1 photon insertion
EXTERN_PROP vector<vprop_t> conf_prop_T; //!< propagator with Tadpole insertion
EXTERN_PROP vector<vprop_t> conf_prop_S; //!< propagator with Scalar insertion
EXTERN_PROP vector<vprop_t> conf_prop_P; //!< propagator with Pseudoscalar insertion

//! return the tag of a prop
string get_prop_tag(size_t im,size_t ir,const string &ins);

//! set the propagators for each conf
void set_conf_props(bool set_QED);

//! read propagators with a certain insertion, for all mass and r
void read_all_mr_INS_props(vector<vprop_t> &conf_prop,const string &template_path,const string &ins);

//! read all propagators for LO
EXTERN_PROP function<void(const string &)> read_all_mr_LO_props INIT_TO(bind(read_all_mr_INS_props,conf_prop_0,_1,"0"));

//! read all propagators for QED
void read_all_mr_QED_props(const string &template_path);

//! read all propagators, LO and QED
void read_all_mr_props(bool read_QED,const string &template_path);

EXTERN_PROP vector<vjprop_t> jprop_0; //!< propagator with no photon
EXTERN_PROP vector<vjprop_t> jprop_2; //!< propagator with 2 photons, Tadpole, Pseudoscalar and possibly Scalar insertions
EXTERN_PROP vector<vjprop_t> jprop_P; //!< propagator with Pseudoscalar inserion
EXTERN_PROP vector<vjprop_t> jprop_S; //!< propagator with Scalar insertion

//! set all jackknifed props
void set_jprops(bool set_QED);

//! takes all conf prop and put them into appropriate jprop
void build_all_mr_jackknifed_props(bool use_QED,size_t ijack);

//! clusterize the propagators
void clusterize_all_mr_props(bool use_QED,size_t clust_size);

//! compute the inverse of all mr props
vector<vjprop_t> get_all_mr_props_inv(const vector<vjprop_t> &jprop);

#undef EXTERN_PROP
#undef INIT_TO

#endif
