#ifndef _PR_BIL_HPP
#define _PR_BIL_HPP

#include <Dirac.hpp>
#include <meas_vec.hpp>

#include <MOM2/geometry.hpp>
#include <MOM2/prop.hpp>

#ifndef EXTERN_PR_BIL
 #define EXTERN_PR_BIL extern
 #define INIT_PR_BIL_TO(...)
#else
 #define INIT_PR_BIL_TO(...) = __VA_ARGS__
#endif

//! bilinears
enum ibil_t{iS,iA,iP,iV,iT};
const string bil_tag="SAPVT";

//! number of bilinears
const size_t nbil=5;

//! holds all ibil_t
const ibil_t ibil_t_list[nbil]={iS,iA,iP,iV,iT};
const vector<vector<size_t>> iG_of_bil={{0},{6,7,8,9},{5},{1,2,3,4},{10,11,12,13,14,15}};

//! holds jackkniffed vertex for an mr combo and for a given mom
class jbil_vert_t
{
  //! return a list of pointers to all internal data
  vector<jqprop_t*> get_all_ptrs()
  {
    vector<jqprop_t*> out={&LO};
    if(pars::use_QED)
      for(auto &p : {&PH,&CR_CT_in,&CR_CT_ou,&TM_CT_in,&TM_CT_ou,&QED})
	out.push_back(p);
    return out;
  }
  
public:
  jqprop_t LO; //!< jackkniffed vertex
  
  jqprop_t PH;       //!< jackkniffed vertex with all photons
  jqprop_t CR_CT_in; //!< jackkniffed vertex with critical counterterm on line in
  jqprop_t CR_CT_ou; //!< jackkniffed vertex with critical counterterm on line out
  jqprop_t TM_CT_in; //!< jackkniffed vertex with twisted counterterm on line in
  jqprop_t TM_CT_ou; //!< jackkniffed vertex with twisted counterterm on line out
  jqprop_t QED;      //!< full correction
  
  //! clusterize
  void clusterize_all(size_t clust_size)
  {
    for(auto &p : get_all_ptrs())
      p->clusterize(clust_size);
  }
};

#undef EXTERN_PR_BIL
#undef INIT_PR_BIL_TO

#endif
