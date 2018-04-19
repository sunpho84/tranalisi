#ifndef _CONTRACTIONS_HPP
#define _CONTRACTIONS_HPP

#include <meas_vec.hpp>

//! read a given contraction
djvec_t get_contraction(const string &combo,const string &ID,vector<size_t> &conf_list,const size_t reim,const int tpar,const string &prop_hadr_path);

//! compute the dshift to m cr
djack_t compute_deltam_cr(vector<size_t> &conf_list,const size_t tmin,const size_t tmax,const size_t im,const size_t nr,const bool use_QED,const string &prop_hadr_path);

//! compute a meson mass
djack_t compute_meson_mass(vector<size_t> &conf_list,const size_t tmin,const size_t tmax,const size_t im1,const size_t im2,const size_t nr,const string &prop_hadr_path);


#endif
