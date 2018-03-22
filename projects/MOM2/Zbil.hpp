#ifndef _ZBIL_HPP
#define _ZBIL_HPP

#include <Dirac.hpp>
#include <meas_vec.hpp>

#include <MOM2/geometry.hpp>
#include <MOM2/prop.hpp>

#ifndef EXTERN_ZBIL
 #define EXTERN_ZBIL extern
 #define INIT_ZBIL_TO(...)
#else
 #define INIT_ZBIL_TO(...) = __VA_ARGS__
#endif

//! bilinears Z
enum iZbil_t{iZS,iZA,iZP,iZV,iZT};
const string Zbil_tag="SAPVT";

//! number of Z
const size_t nZbil=5;

//! holds all iZbil_t
const iZbil_t iZbil_t_list[nZbil]={iZS,iZA,iZP,iZV,iZT};
const vector<vector<size_t>> iG_of_Zbil={{0},{6,7,8,9},{5},{1,2,3,4},{10,11,12,13,14,15}};

//! holds jackkniffed vertex for an mr combo and for a given mom
class jbil_vert_t
{
  //! return a list of pointers to all internal data
  vector<vector<jqprop_t>*> get_all_ptrs(bool use_QED)
  {
    vector<vector<jqprop_t>*> out={&LO};
    if(use_QED)
      for(auto &p : {&PH,&CT_in,&CT_ou,&S,&QED})
	out.push_back(p);
    return out;
  }
  
public:
  vector<jqprop_t> LO; //!< jackkniffed vertex
  
  vector<jqprop_t> PH;    //!< jackkniffed vertex with all photons
  vector<jqprop_t> CT_in; //!< jackkniffed vertex with counterterm on line in
  vector<jqprop_t> CT_ou; //!< jackkniffed vertex with counterterm on line out
  vector<jqprop_t> S;     //!< jackkniffed S vertex
  vector<jqprop_t> QED;   //!< full correction
  
  jbil_vert_t(size_t size,bool use_QED)
  {for(auto &p : get_all_ptrs(use_QED)) p->resize(size);}
  
  //! clusterize
  void clusterize_all(bool use_QED,size_t clust_size,const index_t &im_r_im_r_igam_ind,const djvec_t &deltam_cr)
  {
    vector<vector<jqprop_t>*> ps=get_all_ptrs(use_QED);
    index_t ind({{"type",ps.size()},{"icombo",ps[0]->size()}});
    
#pragma omp parallel for
    for(size_t i=0;i<ind.max();i++)
      {
	const vector<size_t> comp=ind(i);
	const size_t itype=comp[0];
	const size_t icombo=comp[1];
	
	(*ps[itype])[icombo].clusterize(clust_size);
      }
    
    const size_t nm=im_r_im_r_igam_ind.max(0);
    const size_t nr=im_r_im_r_igam_ind.max(1);
    const index_t im_r_ind({{"nm",nm},{"nr",nr}});
    
    if(use_QED)
#pragma omp parallel for
      for(size_t i=0;i<im_r_im_r_igam_ind.max();i++)
	{
	  const vector<size_t> comps=im_r_im_r_igam_ind(i);
	  const size_t im_r_in=im_r_ind({comps[0],comps[1]});
	  const size_t im_r_ou=im_r_ind({comps[2],comps[3]});
	  for(size_t ijack=0;ijack<=njacks;ijack++)
	    QED[i][ijack]=
	      PH[i][ijack]
	      +deltam_cr[im_r_in][ijack]*CT_in[i][ijack]
	      +deltam_cr[im_r_ou][ijack]*CT_ou[i][ijack]
	      ;
	}
  }
};

#undef EXTERN_ZBIL
#undef INIT_ZBIL_TO

#endif
