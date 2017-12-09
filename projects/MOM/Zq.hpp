#ifndef _ZQ_HPP
#define _ZQ_HPP

#include <meas_vec.hpp>

#include <geometry.hpp>
#include <prop.hpp>

//! compute for a single jack
double compute_Zq(const prop_t &prop_inv,const size_t imom);

//! compute for a djack_t
djack_t compute_Zq(const jprop_t &jprop_inv,const size_t imom);

#endif
