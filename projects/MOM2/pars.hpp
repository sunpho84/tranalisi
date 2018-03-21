#ifndef _PARS_HPP
#define _PARS_HPP

#ifndef EXTERN_PARS
 #define EXTERN_PARS extern
 #define INIT_PARS_TO(...)
#else
 #define INIT_PARS_TO(...) __VA_ARGS__
#endif

#include <tranalisi.hpp>

#include <MOM2/types.hpp>

namespace gaz
{
  enum type_t{PLAQ,TLSYM,IWA};
  const map<string,type_t> decr={{"Plaq",PLAQ},{"tlSym",TLSYM},{"Iwa",IWA}};
  const size_t n=3;
  PROVIDE_DECRYPTER;
}

namespace gf
{
  enum type_t{FEYNMAN,LANDAU};
  const map<string,type_t> decr={{"Feynman",FEYNMAN},{"Landau",LANDAU}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
  const double lambda[2]={1.0,0.0};
}

//! chiral extrapolation
namespace chir_extr
{
  enum type_t{MQUARK,MMES};
  const map<string,type_t> decr={{"MQuark",MQUARK},{"MMes",MMES}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
}

//! kind of scheme supported
namespace reno_scheme
{
  enum type_t{RI_MOM,SMOM};
  const map<string,type_t> decr={{"RI_MOM",RI_MOM},{"SMOM",SMOM}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
}

namespace pars
{
  //! gauge action
  EXTERN_PARS gaz::type_t act INIT_PARS_TO({gaz::PLAQ});
  //! extrapolation method
  EXTERN_PARS chir_extr::type_t chir_extr_method INIT_PARS_TO({chir_extr::MQUARK});
  //! compute mesoleptonic vertex
  EXTERN_PARS bool compute_meslep INIT_PARS_TO({false});
  //! compute deltam_cr counterterm and subtract
  EXTERN_PARS bool deltam_cr_ct INIT_PARS_TO({true});
  //! Filter for democracy
  EXTERN_PARS double filter_thresh INIT_PARS_TO({1});
  //! number of momenta between each print
  EXTERN_PARS size_t print_each_mom INIT_PARS_TO({100});
  //! reno scheme used
  EXTERN_PARS reno_scheme::type_t scheme INIT_PARS_TO({reno_scheme::RI_MOM});
  //! use QED
  EXTERN_PARS bool use_QED INIT_PARS_TO({false});
  
  //! list of ensembles
  EXTERN_PARS vector<string> ens;
}

#undef INIT_PARS_TO
#undef EXTERN_PARS

#endif
