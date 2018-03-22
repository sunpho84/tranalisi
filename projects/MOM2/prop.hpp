#ifndef _PROP_HPP
#define _PROP_HPP

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
 #define INIT_PROP_TO(...)
#else
 #define INIT_PROP_TO(...) = __VA_ARGS__
#endif

#include <MOM2/geometry.hpp>
#include <MOM2/pars.hpp>

const double tau3[2]={-1.0,+1.0}; //!< tau entering the propagator

using vprop_t=vector<qprop_t>; //! vector of propagators
using vjqprop_t=vector<jqprop_t>; //! vector of jackkniffed props

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
    if(pars::use_QED) return NPROP_WITH_QED;
    else        return 1;
  }
  
  //! return the coefficient for reading the kind
  static dcompl_t coeff_to_read(const size_t ikind,const size_t r)
  {
    const size_t nmax=m_r_mom_conf_qprops_t::nprop_kind();
    if(ikind>=nmax) CRASH("cannot ask for coeff of kind %zu, maximum is %zu",ikind,nmax);
    
    if(ikind==5) //P
      return dcompl_t(0,-1// tau3[r]
		      );
    else
      if(ikind==4) //S
	return -1.0;
      else //others
	return 1.0;
  }
};

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
  
  void build_all_mr_jackkniffed_qprops(vector<jm_r_mom_qprops_t> &jprops,const vector<m_r_mom_conf_qprops_t> &props);
  
  //! clusterize the propagators
  void clusterize_all_mr_props(bool use_QED,size_t clust_size)
  {
    //cout<<"Clusterizing all props, clust_size="<<clust_size<<endl;
    for(auto &p : get_all_ptrs(use_QED))
      p->clusterize(clust_size);
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
    if(pars::use_QED) return NPROP_WITH_QED;
    else        return 1;
  }
  
  //! returns the inverse propagators
  void get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
			       const vector<jm_r_mom_qprops_t> &jprops,
			       const index_t &im_r_ijackp1_ind) const;
};

//! incorporate charge into the quark propagators
void incorporate_charge(vector<m_r_mom_conf_qprops_t> &props,const double ch);

//! incorporates charge into the lepton propagator
void incorporate_charge(vector<mom_conf_lprops_t> &props,const double ch);

#undef EXTERN_PROP
#undef INIT_PROP_TO

#endif
