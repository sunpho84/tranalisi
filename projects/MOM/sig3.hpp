#ifndef _SIG3_HPP
#define _SIG3_HPP

#include <meas_vec.hpp>

#include <geometry.hpp>
#include <types.hpp>

//! compute the form factor with gamma5
djvec_t compute_sig3(const vjprop_t &jprop_inv);

#endif
