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
  vector<vector<jprop_t>*> get_all_ptrs(bool use_QED)
  {
    vector<vector<jprop_t>*> out={&LO};
    if(use_QED)
      for(auto &p : {&EM,&P,&S})
	out.push_back(p);
    return out;
  }
  
public:
  vector<jprop_t> LO; //!< jackkniffed vertex
  
  vector<jprop_t> EM; //!< jackkniffed vertex with all em corrs
  vector<jprop_t> P; //!< jackkniffed P vertex
  vector<jprop_t> S; //!< jackkniffed S vertex
  
  jbil_vert_t(size_t size,bool use_QED)
  {for(auto &p : get_all_ptrs(use_QED)) p->resize(size);}
  
  //! clusterize
  void clusterize_all(bool use_QED,size_t clust_size)
  {
    //cout<<"Clusterizing all verts, clust_size="<<clust_size<<endl;
    for(auto &p : get_all_ptrs(use_QED))
      for(auto &pi : *p)
	clusterize(pi,clust_size);
  }
};

//! compute the projected bilinears
djvec_t compute_proj_bil(const vjprop_t &jprop_inv1,const vector<jprop_t> &jverts,const vjprop_t &jprop_inv2,const index_t &im_r_ind);

//! compute all vertices for a certain conf
void build_all_mr_gbil_jackknifed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_props_t> &props,
					const index_t &im_r_im_r_igam_ind,const index_t &im_r_ijack_ind,bool use_QED);

//! finish computing the em verts
void finish_jverts_EM(jbil_vert_t &jverts,const djack_t &deltam_cr);

#undef EXTERN_ZBIL
#undef INIT_TO

#endif
