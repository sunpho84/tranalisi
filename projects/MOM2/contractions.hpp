#ifndef _CONTRACTIONS_HPP
#define _CONTRACTIONS_HPP

#include <tranalisi.hpp>

const int ODD=-1,EVN=1;

//! compute the dshift to m cr
djack_t compute_deltam_cr(vector<size_t> &conf_list,const size_t tmin,const size_t tmax,const size_t im,const size_t nr,const string &prop_hadr_path);

//! compute a meson mass
djack_t compute_meson_mass(vector<size_t> &conf_list,const size_t tmin,const size_t tmax,const size_t im1,const size_t im2,const size_t nr,const string &prop_hadr_path);

#endif
