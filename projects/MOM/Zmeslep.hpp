#ifndef _ZMESLEP_HPP
#define _ZMESLEP_HPP

#include <Dirac.hpp>

#include <prop.hpp>

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
  {for(auto &p : get_all_ptrs()) p->resize(size);}
  
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

void build_all_mr_gmeslep_jackknifed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
					   const vector<mom_conf_lprops_t> &props_lep,
					   const index_t &im_r_im_r_igam_iprojlep_ind,const index_t &im_r_iclust_ind,const double q1,const double q2);

#endif