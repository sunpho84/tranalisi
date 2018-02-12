#ifndef _ZMESLEP_HPP
#define _ZMESLEP_HPP

#include <Dirac.hpp>

#include <prop.hpp>

#ifndef EXTERN_MESLEP
 #define EXTERN_MESLEP extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) =__VA_ARGS__
#endif

EXTERN_MESLEP bool compute_meslep;

namespace meslep
{
  EXTERN_MESLEP const vector<vector<size_t>> iGq_of_iop        INIT_TO({{ 1, 2, 3, 4}, { 1, 2, 3, 4}, { 0}, { 0}, {10,11,12,13,14,15}});
  EXTERN_MESLEP const vector<int>            g5_sign_of_iop    INIT_TO({ -1,            +1,            -1,   +1 ,  +1});
  EXTERN_MESLEP const vector<vector<size_t>> &iGl_of_iop       INIT_TO(iGq_of_iop);
  EXTERN_MESLEP const vector<vector<size_t>> &ipGl_of_iproj    INIT_TO(iGl_of_iop);
  EXTERN_MESLEP const vector<int>            &g5_sign_of_iproj INIT_TO(g5_sign_of_iop);
}

//! holds jackkniffed vertex for an mr combo and for a given mom
class jmeslep_vert_t
{
  //! return a list of pointers to all internal data
  vector<vector<jqprop_t>*> get_all_ptrs()
  {
    vector<vector<jqprop_t>*> out={&ML1,&ML2};
    return out;
  }
  
public:
  vector<jqprop_t> ML1; //!< jackkniffed mesolep vertex, photon going out from line 1
  vector<jqprop_t> ML2; //!< jackkniffed mesolep vertex, photon going out from line 2
  
  jmeslep_vert_t(size_t size)
  {
    for(auto &p : get_all_ptrs()) p->resize(size);
  }
  
  //! clusterize
  void clusterize_all(size_t clust_size)
  {
    vector<vector<jqprop_t>*> ps=get_all_ptrs();
    index_t ind({{"type",ps.size()},{"icombo",ps[0]->size()}});
    
#pragma omp parallel for
    for(size_t i=0;i<ind.max();i++)
      {
	const vector<size_t> comp=ind(i);
	const size_t itype=comp[0];
	const size_t icombo=comp[1];
	
	(*ps[itype])[icombo].clusterize(clust_size);
      }
  }
};

void build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
					   const vector<mom_conf_lprops_t> &props_lep,
					   const index_t &im_r_ind);

djvec_t compute_proj_measlep(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind);

#undef EXTERN_MESLEP
#undef INIT_TO

#endif
