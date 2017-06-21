#ifndef _ZBIL_HPP
#define _ZBIL_HPP

#include <Dirac.hpp>
#include <meas_vec.hpp>

#include <geometry.hpp>
#include <prop.hpp>

#ifndef EXTERN_ZBIL
 #define EXTERN_ZBIL extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using vjbil_vert_t=vector<jprop_t>; //! vector of jackkniffed bilinear vertex (spanning moms)

//! bilinears Z
enum iZbil{iZS,iZA,iZP,iZV,iZT};
const string Zbil_tag="SAPVT";

//! number of Z
const size_t nZbil=5;

//! compute the projected bilinears
vector<djvec_t> compute_proj_bil(const vjprop_t &jprop_inv1,vjbil_vert_t &jverts,const vjprop_t &jprop_inv2);

//! compute a certain vertex, with a given pair of propagators
void build_jackknifed_vert_Gamma(vjbil_vert_t &jverts,const vprop_t &prop1,size_t iG,const vprop_t &prop2,size_t ijack);

//! compute all vertices
void build_all_mr_gBil_jackknifed_verts(bool use_QED,size_t ijack);

//! clusterize the vertexes
void clusterize_verts(vjbil_vert_t &jverts,size_t clust_size=1);

//! index of the (im,ir)(im,ir)(iZbil) prop combo, first is reverted
EXTERN_ZBIL index_t mr_Zbil_ind;
void set_mr_Zbil_ind(size_t nm,size_t nr);

//! index of the (im,ir)(im,ir)(ig) prop combo, first is reverted
EXTERN_ZBIL index_t mr_gbil_ind;
void set_mr_gbil_ind(size_t nm,size_t nr);

EXTERN_ZBIL size_t nmr_Zbil; //!< total number of m and r pair times nZbil
EXTERN_ZBIL size_t nmr_gbil; //!< total number of m and r pair times nGamma

EXTERN_ZBIL vector<vjbil_vert_t> jverts_0; //!< jackkniffed vertex
EXTERN_ZBIL vector<vjbil_vert_t> jverts_em; //!< jackkniffed vertex with all em corrs
EXTERN_ZBIL vector<vjbil_vert_t> jverts_P; //!< jackkniffed P vertex
EXTERN_ZBIL vector<vjbil_vert_t> jverts_S; //!< jackkniffed S vertex

//! set the vertex combo
void set_jbil_verts(bool set_QED);

//! build all vertices
void build_all_mr_gbil_jackknifed_verts(bool use_QED,size_t ijack);

//! clusterize all vertices
void clusterize_all_mr_gbil_verts(bool use_QED,size_t clust_size);

#undef EXTERN_ZBIL
#undef INIT_TO

#endif
