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
  struct mesloop_t
  {
    vector<dcompl_t> LO;
    vector<dcompl_t> F;
    
    mesloop_t(size_t size) : LO(size),F(size)
    {
    }
  };
  
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
    const size_t norm;
  };
  
  EXTERN_MESLEP vector<Zop_t> zops INIT_TO({
      {{{1,1},{2,2},{3,3},{4,4}},-1,4},
      {{{1,1},{2,2},{3,3},{4,4}},+1,4},
      {{{ 0,0}},-1,1},
      {{{ 0,0}},+1,1},
      {{{{5,10},{6,11},{7,12},{8,13},{9,14},{10,15}}},+1,12}});
  
  EXTERN_MESLEP size_t nZop INIT_TO({5});
  
  EXTERN_MESLEP       vector<size_t>         listGl   INIT_TO(={{ 0, 1, 2, 3, 4,10,11,12,13,14,15}}); //!< list of Gamma to be inserted on lepton side of the operator
  EXTERN_MESLEP       vector<size_t>         &listpGl INIT_TO(=listGl);
  EXTERN_MESLEP       vector<int>            Lg5_sign INIT_TO(={{+1,-1,-1,-1,-1,+1,+1,+1,+1,+1,+1}}); //!< 1+sign*g5 on the quark side
}

//! holds jackkniffed vertex for an mr combo and for a given mom
class jmeslep_vert_t
{
  //! return a list of pointers to all internal data
  vector<vector<jqprop_t>*> get_all_ptrs()
  {
    return {&LO,&PH,&CT1,&CT2,&S,&QED};
  }
  
public:
  vector<jqprop_t> LO; //!< jackkniffed mesolep vertex, no photon
  
  vector<jqprop_t> PH; //!< jackkniffed vertex, nasty + self + tad
  vector<jqprop_t> CT1; //!< jackkniffed meson-like vertex, counterterm on line 1
  vector<jqprop_t> CT2; //!< jackkniffed meson-like vertex, counterterm on line 2
  vector<jqprop_t> S;   //!< jackkniffed meson-like vertex, S insertion
  
  vector<jqprop_t> QED; //!< full correction
  
  jmeslep_vert_t(size_t size)
  {
    for(auto &p : get_all_ptrs()) p->resize(size);
  }
  
  //! clusterize
  void clusterize_all(const size_t clust_size,const index_t &im_r_im_r_iop_ilistpGl_ind,const djvec_t &deltam_cr)
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
    
#pragma omp parallel for
    for(size_t i=0;i<im_r_im_r_iop_ilistpGl_ind.max();i++)
	{
	  const vector<size_t> comps=im_r_im_r_iop_ilistpGl_ind(i);
	  const size_t im1=comps[0],im2=comps[2];
	  for(size_t ijack=0;ijack<=njacks;ijack++)
	    QED[i][ijack]=
	      PH[i][ijack]
	      // #warning not including counterterm
	      -deltam_cr[im1][ijack]*CT1[i][ijack]
	      -deltam_cr[im2][ijack]*CT2[i][ijack]
	      ;
	}
  }
};

void build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
					   const vector<mom_conf_lprops_t> &props_lep,
					   const index_t &im_r_ind);

djvec_t compute_proj_meslep(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind);

#undef EXTERN_MESLEP
#undef INIT_TO

#endif
