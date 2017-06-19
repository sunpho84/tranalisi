#ifndef _ZBIL_HPP
#define _ZBIL_HPP

#include <meas_vec.hpp>

#include <geometry.hpp>
#include <types.hpp>

#ifndef EXTERN_ZBIL
 #define EXTERN_ZBIL extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

//! bilinears Z
enum iZbil{iZS,iZA,iZP,iZV,iZT};
const string Zbil_tag="SAPVT";

//! number of Z
const size_t nZbil=5;

//! compute the projected bilinears
vector<djvec_t> compute_proj_bil(const vjprop_t &jprop_inv1,vector<jverts_t> &jverts,const vjprop_t &jprop_inv2);

//! compute all 16 vertices, with a given pair of propagators
void build_jackknifed_verts(vector<jverts_t> &jverts,const vprop_t &prop1,const vprop_t &prop2,size_t ijack);

//! clusterize the vertexes
void clusterize_verts(vector<jverts_t> &jverts,size_t clust_size=1);

#undef EXTERN_ZBIL
#undef INIT_TO

#endif
