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
enum iZbil_t{iZS,iZA,iZP,iZV,iZT};
const string Zbil_tag="SAPVT";

//! number of Z
const size_t nZbil=5;

//! holds all iZbil_t
const iZbil_t iZbil_t_list[nZbil]={iZS,iZA,iZP,iZV,iZT};

//! holds jackkniffed vertex for an mr combo and for a given mom
class jbil_vert_t
{
  //! return a list of pointers to all internal data
  vector<vector<jqprop_t>*> get_all_ptrs(bool use_QED)
  {
    vector<vector<jqprop_t>*> out={&LO};
    if(use_QED)
      for(auto &p : {&PH,&CT1,&CT2,&S,&QED})
	out.push_back(p);
    return out;
  }
  
public:
  vector<jqprop_t> LO; //!< jackkniffed vertex
  
  vector<jqprop_t> PH;  //!< jackkniffed vertex with all photons
  vector<jqprop_t> CT1; //!< jackkniffed vertex with counterterm on line 1
  vector<jqprop_t> CT2; //!< jackkniffed vertex with counterterm on line 2
  vector<jqprop_t> S;   //!< jackkniffed S vertex
  vector<jqprop_t> QED; //!< full correction
  
  jbil_vert_t(size_t size,bool use_QED)
  {for(auto &p : get_all_ptrs(use_QED)) p->resize(size);}
  
  //! clusterize
  void clusterize_all(bool use_QED,size_t clust_size,const index_t &im_r_im_r_igam_ind,const djvec_t &deltam_cr)
  {
    vector<vector<jqprop_t>*> ps=get_all_ptrs(use_QED);
    index_t ind({{"type",ps.size()},{"icombo",ps[0]->size()}});
    
    //cout<<"Clusterizing all verts, clust_size="<<clust_size<<endl;
#pragma omp parallel for
    for(size_t i=0;i<ind.max();i++)
      {
	const vector<size_t> comp=ind(i);
	const size_t itype=comp[0];
	const size_t icombo=comp[1];
	
	(*ps[itype])[icombo].clusterize(clust_size);
      }
    
    if(use_QED)
#pragma omp parallel for
      for(size_t i=0;i<im_r_im_r_igam_ind.max();i++)
	{
	  const vector<size_t> comps=im_r_im_r_igam_ind(i);
	  const size_t im1=comps[0],im2=comps[2];
	  for(size_t ijack=0;ijack<=njacks;ijack++)
	    QED[i][ijack]=
	      PH[i][ijack]-deltam_cr[im1][ijack]*CT1[i][ijack]-deltam_cr[im2][ijack]*CT2[i][ijack];
	}
  }
};

//! compute the projected bilinears
djvec_t compute_proj_bil(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind);

//! compute all vertices for a certain conf
void build_all_mr_gbil_jackkniffed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
					const index_t &im_r_im_r_igam_ind,const index_t &im_r_ijack_ind,bool use_QED);

#undef EXTERN_ZBIL
#undef INIT_TO

#endif
