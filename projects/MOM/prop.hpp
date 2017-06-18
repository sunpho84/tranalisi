#ifndef _PROP_HPP
#define _PROP_HPP

#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

//! compute the inverse of the passed propagator
vjprop_t get_prop_inv(const vjprop_t &jprop);

//! read a propagator from a given file
vprop_t read_prop(const string &template_path,size_t iconf);

//! add the prop of a given conf on the jackknife
void build_jackknifed_prop(vjprop_t &jprop,const vprop_t &prop,size_t ijack);

//! clusterize the propagator
void clusterize_prop(vjprop_t &jprop,size_t clust_size=1);

#endif
