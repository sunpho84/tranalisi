#ifndef _GEVP_HPP
#define _GEVP_HPP

#include <meas_vec.hpp>

//! Solves the generalized eigenvalue problem for a matrix d X d of correlators
vector<djvec_t> gevp(const vector<djvec_t> &d,const size_t t0);

#endif
