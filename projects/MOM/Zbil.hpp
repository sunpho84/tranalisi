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

//! bilinears Z
enum iZbil{iZS,iZA,iZP,iZV,iZT};
const string Zbil_tag="SAPVT";

//! number of Z
const size_t nZbil=5;

//! holds jackkniffed vertex for an mr combo and for a given mom
class jbil_vert_t
{
  //! return a list of pointers to all internal data
  vector<jprop_t*> get_all_ptrs(bool use_QED)
  {
    vector<jprop_t*> out={&jvert_0};
    if(use_QED)
      for(auto &p : {&jvert_em,&jvert_P,&jvert_S})
	out.push_back(p);
    return out;
  }
  
public:
  jprop_t jvert_0; //!< jackkniffed vertex
  
  jprop_t jvert_em; //!< jackkniffed vertex with all em corrs
  jprop_t jvert_P; //!< jackkniffed P vertex
  jprop_t jvert_S; //!< jackkniffed S vertex
  
  //! clusterize
  void clusterize_all(bool use_QED,size_t clust_size)
  {
    cout<<"Clusterizing all verts, clust_size="<<clust_size<<endl;
    for(auto &p : get_all_ptrs(use_QED))
      clusterize(*p,clust_size);
  }
};

//! compute the projected bilinears
//vector<djvec_t> compute_proj_bil(const vjprop_t &jprop_inv1,vjbil_vert_t &jverts,const vjprop_t &jprop_inv2);

//! compute all vertices for a certain conf
void build_all_mr_gbil_jackknifed_verts(vector<jbil_vert_t> &jbil,const vector<m_r_mom_conf_props_t> &props,
					const index_t &im_r_im_r_igam_ind,const index_t &im_r_ijack_ind,bool use_QED);

//! clusterize all vertices
void clusterize_all_mr_gbil_verts(vector<jbil_vert_t> &jbils,bool use_QED,size_t clust_size);

#undef EXTERN_ZBIL
#undef INIT_TO

#endif
