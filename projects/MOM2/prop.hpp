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

namespace qprop
{
  enum ins{LO ,FF , F , T , S ,P };
  EXTERN_PROP vector<ins>    ins_list;
  EXTERN_PROP vector<string> ins_tag;
  EXTERN_PROP size_t nins;
}

namespace lprop
{
  enum ins{LO , F};
  EXTERN_PROP vector<ins>    ins_list;
  EXTERN_PROP vector<string> ins_tag;
  EXTERN_PROP size_t nins;
}

//! return the coefficient for reading the kind
dcompl_t coeff_to_read(const lprop::ins ikind,const size_t r);

namespace jqprop
{
  enum ins{LO , PH , CR_CT , TM_CT };
  EXTERN_PROP vector<ins>    ins_list;
  EXTERN_PROP vector<string> ins_tag;
  EXTERN_PROP size_t nins;
}

//! returns the inverse propagators
void get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
			     const vector<jqprop_t> &jprops,
			     const index_t &im_r_ijackp1_ind);

#undef EXTERN_PROP
#undef INIT_PROP_TO

#endif
