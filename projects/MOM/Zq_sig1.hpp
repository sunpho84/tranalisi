#ifndef _ZQ_SIG1_HPP
#define _ZQ_SIG1_HPP

#include <meas_vec.hpp>

#include <geometry.hpp>
#include <types.hpp>

//! compute Zq through sigma1 for a given jack
double compute_Zq_sig1(const prop_t &prop_inv,size_t imom);

//! compute Zq through sigma1 for a given momentum
djack_t compute_Zq_sig1(const jprop_t &jprop_inv,size_t imom);

#endif
