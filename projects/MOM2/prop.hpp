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
  enum ins{LO ,F , FF , T , S ,P, QED, RI, RI_QED};
  EXTERN_PROP vector<ins>    ins_list;
  EXTERN_PROP vector<size_t> iins_of_ins;
  EXTERN_PROP vector<string> ins_tag INIT_PROP_TO({"0" ,"F","FF","T","S","P","QED","RI","RI_QED"});
  
  EXTERN_PROP size_t nins;
  
  void set_ins();
}

namespace lprop
{
  enum ins{LO , F};
  EXTERN_PROP vector<ins>    ins_list;
  EXTERN_PROP vector<string> ins_tag INIT_PROP_TO({"0" ,"F"});
  
  EXTERN_PROP size_t nins;
  
  void set_ins();
}

//! return the coefficient for reading the kind
dcompl_t coeff_to_read(const lprop::ins ikind,const size_t r);

namespace jqprop
{
  enum ins{LO , PH , CR , TM , QED, RI, RI_QED};
  EXTERN_PROP vector<ins>    ins_list;
  EXTERN_PROP vector<size_t> iins_of_ins;
  EXTERN_PROP vector<string> ins_tag INIT_PROP_TO({"0" ,"PH","CR","TM","QED","RI","RI_QED"});
  EXTERN_PROP size_t nins;
  
  void set_ins();
}

//! returns the coefficient to insert for each insertion
dcompl_t coeff_to_read(const qprop::ins ins,const size_t r);

#undef EXTERN_PROP
#undef INIT_PROP_TO

#endif
