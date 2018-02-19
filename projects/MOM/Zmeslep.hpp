#ifndef _ZMESLEP_HPP
#define _ZMESLEP_HPP

#include <Dirac.hpp>

#include <prop.hpp>

#ifndef EXTERN_MESLEP
 #define EXTERN_MESLEP extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) __VA_ARGS__
#endif

EXTERN_MESLEP bool compute_meslep;

namespace meslep
{
  //! holds info on constructing operators
  struct Zop_t
  {
    struct listGl_Gq_t
    {
      const size_t ilistGl; //!< insertions to be taken from the list
      const size_t Gq;      //!< basic gamma on the quark side of the operator...
    };
    
    const vector<listGl_Gq_t> contr;
    const int Qg5_sign;               //!< 1+sign*g5 on the quark side
  };
  
  EXTERN_MESLEP vector<Zop_t> zops INIT_TO({
      {{{1,1},{2,2},{3,3},{4,4}},-1},
      {{{1,1},{2,2},{3,3},{4,4}},+1},
      {{{ 0,0}},-1},
      {{{ 0,0}},+1},
      {{{{5,10},{6,11},{7,12},{8,13},{9,14},{10,15}}},+1}});
  
  EXTERN_MESLEP size_t nZop INIT_TO({5});
  
  EXTERN_MESLEP       vector<size_t>         listGl            INIT_TO(={{0,1,2,3,4,10,11,12,13,14,15}}); //!< list of Gamma to be inserted on lepton side of the operator
  EXTERN_MESLEP       vector<size_t>         &listpGl INIT_TO(=listGl);
}

//! holds jackkniffed vertex for an mr combo and for a given mom
class jmeslep_vert_t
{
  //! return a list of pointers to all internal data
  vector<vector<jqprop_t>*> get_all_ptrs()
  {
    vector<vector<jqprop_t>*> out={&QCD,&ML1,&ML2};
    return out;
  }
  
public:
  vector<jqprop_t> QCD; //!< jackkniffed mesolep vertex, no photon
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
