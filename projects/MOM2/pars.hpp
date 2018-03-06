#ifndef _PARS_HPP
#define _PARS_HPP

#ifndef EXTERN_PARS
 #define EXTERN_PARS extern
 #define INIT_PARS_TO(...)
#else
 #define INIT_PARS_TO(...) __VA_ARGS__
#endif

#include <tranalisi.hpp>

#include <types.hpp>

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

namespace pars
{
  //! number of momenta between each print
  EXTERN_PARS size_t print_each_mom INIT_PARS_TO({100});
  //! gauge action
  EXTERN_PARS gaz::type_t act INIT_PARS_TO({gaz::PLAQ});
  //! use QED
  EXTERN_PARS bool use_QED INIT_PARS_TO({false});
  //! list of ensembles
  EXTERN_PARS vector<string> ens;
}

#undef INIT_PARS_TO
#undef EXTERN_PARS

#endif
