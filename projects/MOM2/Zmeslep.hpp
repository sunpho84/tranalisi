#ifndef _ZMESLEP_HPP
#define _ZMESLEP_HPP

#include <Dirac.hpp>

#include <MOM2/prop.hpp>

#ifndef EXTERN_MESLEP
 #define EXTERN_MESLEP extern
 #define INIT_MESLEP_TO(...)
#else
 #define INIT_MESLEP_TO(...) __VA_ARGS__
#endif

namespace meslep
{
  //these are the charges in the lagrangian
  const double ql=-1.0;        //!< the program simulates muon *antiparticle*
  const double q_in=+2.0/3.0;  //!< charge of the quark which comes into the vertex
  const double q_ou=-1.0/3.0;  //!< charge of the quark which comes out the vertex
  
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
    const size_t norm;                //!< norm when inserting the operator
    const size_t pnorm;               //!< norm when projecting
  };
  
  EXTERN_MESLEP vector<Zop_t> zops INIT_MESLEP_TO({
      {{{1,1},{2,2},{3,3},{4,4}},-1,1,4},
      {{{1,1},{2,2},{3,3},{4,4}},+1,1,4},
      {{{ 0,0}},-1,1,1},
      {{{ 0,0}},+1,1,1},
      {{{{5,10},{6,11},{7,12},{8,13},{9,14},{10,15}}},+1,2,24}});
  
  const size_t nZop=5;
  
  EXTERN_MESLEP       vector<size_t>         listGl   INIT_MESLEP_TO(={{ 0, 1, 2, 3, 4,10,11,12,13,14,15}}); //!< list of Gamma to be inserted on lepton side of the operator
  EXTERN_MESLEP       vector<size_t>         &listpGl INIT_MESLEP_TO(=listGl);
  EXTERN_MESLEP       vector<int>            Lg5_sign INIT_MESLEP_TO(={{+1,-1,-1,-1,-1,+1,+1,+1,+1,+1,+1}}); //!< 1+sign*g5 on the quark side
}

//! holds jackkniffed vertex for an mr combo and for a given mom
class jmeslep_vert_t
{
  //! return a list of pointers to all internal data
  vector<vector<jqprop_t>*> get_all_ptrs()
  {
    return {&LO,&PH,&CR_CT_IN,&CR_CT_OU,&TM_CT_IN,&TM_CT_OU,&QED};
  }
  
public:
  vector<jqprop_t> LO; //!< jackkniffed mesolep vertex, no photon
  
  vector<jqprop_t> PH;        //!< jackkniffed vertex, nasty + self + tad
  vector<jqprop_t> CR_CT_IN;  //!< jackkniffed meson-like vertex, critical counterterm on line IN
  vector<jqprop_t> CR_CT_OU;  //!< jackkniffed meson-like vertex, critical counterterm on line OUT
  vector<jqprop_t> TM_CT_IN;  //!< jackkniffed meson-like vertex, twisted counterterm on line IN
  vector<jqprop_t> TM_CT_OU;  //!< jackkniffed meson-like vertex, twisted counterterm on line OUT
  
  vector<jqprop_t> QED; //!< full correction
  
  jmeslep_vert_t(size_t size)
  {
    for(auto &p : get_all_ptrs()) p->resize(size);
  }
  
  //! clusterize
  void clusterize_all(const size_t clust_size,const index_t &im_r_im_r_iop_ilistpGl_ind,const djvec_t &deltam_cr,const djvec_t &deltam_tm)
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
    
    const size_t nm=im_r_im_r_iop_ilistpGl_ind.max(0);
    const size_t nr=im_r_im_r_iop_ilistpGl_ind.max(1);
    const index_t im_r_ind({{"nm",nm},{"nr",nr}});
    
#pragma omp parallel for
    for(size_t i=0;i<im_r_im_r_iop_ilistpGl_ind.max();i++)
	{
	  const vector<size_t> comps=im_r_im_r_iop_ilistpGl_ind(i);
	  const size_t im_r_in=im_r_ind({comps[0],comps[1]});
	  const size_t im_r_ou=im_r_ind({comps[2],comps[3]});
	  for(size_t ijack=0;ijack<=njacks;ijack++)
	    QED[i][ijack]=
	      PH[i][ijack]
	      +deltam_cr[im_r_in][ijack]*CR_CT_IN[i][ijack]
	      +deltam_cr[im_r_ou][ijack]*CR_CT_OU[i][ijack]
	      +deltam_tm[im_r_in][ijack]*TM_CT_IN[i][ijack]
	      +deltam_tm[im_r_ou][ijack]*TM_CT_OU[i][ijack]
	      ;
	}
  }
};

#undef EXTERN_MESLEP
#undef INIT_MESLEP_TO

#endif
