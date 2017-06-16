#ifndef _ZQ_HPP
#define _ZQ_HPP

#include <meas_vec.hpp>

#include <geometry.hpp>
#include <types.hpp>

//! compute Zq
djvec_t compute_Zq(const vjprop_t &jprop_inv);

#endif
