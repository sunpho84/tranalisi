#ifndef _SIGMA_HPP
#define _SIGMA_HPP

#include <tranalisi.hpp>

//! holds jackkniffed sigma for an mr and all moms
class jsigma_t
{
  //! return a list of pointers to all internal data
  vector<djvec_t*> get_all_ptrs(bool use_QED)
  {
    vector<djvec_t*> out={&LO};
    if(use_QED)
      for(auto &p : {&PH,&CR_CT,&TM_CT})
	out.push_back(p);
    return out;
  }
  
public:
  djvec_t LO; //!< jackkniffed LO
  
  djvec_t PH;    //!< jackkniffed sigma with all photons
  djvec_t CR_CT; //!< jackkniffed sigma with critical counterterm
  djvec_t TM_CT; //!< jackkniffed sigma with twisted counterterm
  
  jsigma_t()
  {}
  
  void resize(size_t size,bool use_QED)
  {
    for(auto &p : get_all_ptrs(use_QED)) p->resize(size);
  }
  
  jsigma_t(size_t size,bool use_QED)
  {
    resize(size,use_QED);
  }
  
  //! clusterize
  void clusterize_all(bool use_QED,size_t clust_size,const index_t &im_r_ind)
  {
    vector<djvec_t*> ps=get_all_ptrs(use_QED);
    index_t ind({{"type",ps.size()},{"im_r",ps[0]->size()}});
    
#pragma omp parallel for
    for(size_t i=0;i<ind.max();i++)
      {
	const vector<size_t> comp=ind(i);
	const size_t itype=comp[0];
	const size_t im_r=comp[1];
	
	(*ps[itype])[im_r].clusterize(clust_size);
      }
  }
};
#endif
