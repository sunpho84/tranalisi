#ifndef _ZQ_SIG1_HPP
#define _ZQ_SIG1_HPP

#include <meas_vec.hpp>

#include <geometry.hpp>
#include <types.hpp>

//! compute Zq through sigma1
djvec_t compute_Zq_sig1(const vjprop_t &jprop_inv);

#endif
