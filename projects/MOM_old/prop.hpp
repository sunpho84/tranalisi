#ifndef _PROP_HPP
#define _PROP_HPP

#include <geometry.hpp>

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using vprop_t=vector<qprop_t>; //! vector of propagators
using vjqprop_t=vector<jqprop_t>; //! vector of jackkniffed props

namespace glb
{
  EXTERN_PROP size_t nm; //!< number of masses
  EXTERN_PROP vector<double> am; //! quark mass
  EXTERN_PROP size_t nr; //!< number of r
  EXTERN_PROP size_t nmr; //!< total number of m and r
  
  EXTERN_PROP index_t im_r_ind; //!< index of im,r
}

EXTERN_PROP index_t i_in_clust_ihit_ind; //!< index of i_in_clust,ihit
EXTERN_PROP index_t conf_ind; //!< index of a conf given ijack and i_in_clust

EXTERN_PROP string suff_hit INIT_TO(""); //!< suffix for each source

EXTERN_PROP bool use_QED; //!< perform or not the QED analysis

const double tau3[2]={-1.0,+1.0}; //!< tau entering the propagator

//! read a quark propagator
void read_qprop(qprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom,const int r_si,const int r_so);

// read a lepton propagator
void read_lprop(lprop_t *prop,raw_file_t &file,const dcompl_t &fact,const size_t imom,const int r_si,const int r_so);

//! holds a given m, r and momentum, for a fixed conf, for quarks
class m_r_mom_conf_qprops_t
{
public:
  qprop_t LO; //!< propagator with no insertion
  
  qprop_t FF; //!< propagator with 2 photon insertions
  qprop_t F;  //!< propagator with 1 photon insertion
  qprop_t T;  //!< propagator with Tadpole insertion
  qprop_t S;  //!< propagator with Scalar insertion
  qprop_t P;  //!< propagator with Pseudoscalar insertion
  
  //! number of all kinds
  static const size_t NPROP_WITH_QED=6;
  //! tag to read
  static const char tag[NPROP_WITH_QED][3];
  //! kind of propagator
  array<qprop_t*,NPROP_WITH_QED> kind{&LO,&FF,&F,&T,&S,&P};
  
  //! number of propagator kind
  static size_t nprop_kind()
  {
    if(use_QED) return NPROP_WITH_QED;
    else        return 1;
  }
  
  //! return the coefficient for reading the kind
  static dcompl_t coeff_to_read(const size_t ikind,const size_t r)
  {
    const size_t nmax=m_r_mom_conf_qprops_t::nprop_kind();
    if(ikind>=nmax) CRASH("cannot ask for coeff of kind %zu, maximum is %zu",ikind,nmax);
    
    if(ikind==5) //P
      return dcompl_t(0,tau3[r]);
    else
      if(ikind==4) //S
	return -1.0;
      else //others
	return 1.0;
  }
};

class mom_conf_lprops_t
{
public:
  lprop_t LO; //!< propagator with no insertion
  
  lprop_t F;  //!< propagator with 1 photon insertions, amputated
  
  //! number of all kinds
  static const size_t NPROP_WITH_QED=2;
  //! tag to read
  static const char tag[NPROP_WITH_QED][3];
  //! kind of propagator
  array<lprop_t*,NPROP_WITH_QED> kind{&LO,&F};
  
  //! number of propagator kind
  static size_t nprop_kind()
  {
    return NPROP_WITH_QED;
  }
};

//! return the tag of a quark prop
string get_qprop_tag(const size_t im,const size_t ir,const size_t ikind);

//! return the tag of a lepton prop
string get_lprop_tag(const size_t ikind);

//! holds a given m and r, all jackks, for quarks
class jm_r_mom_qprops_t
{
  //! return a list of pointers to all internal data
  vector<jqprop_t*> get_all_ptrs(bool use_QED)
  {
    vector<jqprop_t*> out={&LO};
    if(use_QED)
      for(auto &p : {&PH,&CT,&S,&QED})
	out.push_back(p);
    return out;
  }
  
public:
  jqprop_t LO; //!< propagator with no photon, given momentum, all m and r
  
  jqprop_t PH;  //!< propagator with 2 photons, Tadpole
  jqprop_t CT;  //!< propagator with Pseudoscalar inserion
  jqprop_t S;   //!< propagator with Scalar insertion
  jqprop_t QED; //!< propagator with 2 photons, Tadpole, Pseudoscalar and possibly Scalar insertions
  
  //! clusterize the propagators
  void clusterize_all_mr_props(bool use_QED,size_t clust_size)
  {
    //cout<<"Clusterizing all props, clust_size="<<clust_size<<endl;
    for(auto &p : get_all_ptrs(use_QED))
      p->clusterize(clust_size);
  }
};

//! add the charge to the propagators of quarks
void incorporate_charge(vector<m_r_mom_conf_qprops_t> &props,const double ch);

//! add the charge to the propagators of leptons
void incorporate_charge(vector<mom_conf_lprops_t> &props,const double ch);

//! add the prop of a given conf on the jackknife
void build_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,const vector<m_r_mom_conf_qprops_t> &props,bool set_QED,const index_t &im_r_ind);

//! clusterize all props
void clusterize_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,bool use_QED,size_t clust_size,const index_t &im_r_ind,const djvec_t &deltam_cr);

//! return the propagator at maximal twist
lprop_t free_prop(const imom_t &pi,double mu,double kappa,size_t r);

//! returns the inverse propagators
void get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
			     const vector<jm_r_mom_qprops_t> &jprops,
			     const index_t &im_r_ijackp1_ind);

#undef EXTERN_PROP
#undef INIT_TO

#endif
