#ifndef _CONTRACTIONS_HPP
#define _CONTRACTIONS_HPP

#include <meas_vec.hpp>

//! read a given contraction
djvec_t get_contraction(const string &combo,const string &ID,vector<size_t> &conf_list,const size_t reim,const int tpar);

//! compute the dshift to m cr
djack_t compute_deltam_cr(vector<size_t> &conf_list,const size_t tmin,const size_t tmax,const size_t im,const size_t nr);


#endif
