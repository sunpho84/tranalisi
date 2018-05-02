#ifndef _SIGMA_HPP
#define _SIGMA_HPP

#include <tranalisi.hpp>

#ifndef EXTERN_SIGMA
 #define EXTERN_SIGMA extern
 #define INIT_SIGMA_TO(...)
#else
 #define INIT_SIGMA_TO(...) __VA_ARGS__
#endif

namespace sigma
{
  enum SIG_ins{LO,  CR,  TM,  PH};
  EXTERN_SIGMA vector<SIG_ins>  SIG_ins_list;
  EXTERN_SIGMA vector<string>   SIG_ins_tag;
  EXTERN_SIGMA size_t nsig_ins;
}

//! set all sigma insertions
void set_sigma_ins();

#undef EXTERN_SIGMA
#undef INIT_SIGMA_TO

#endif
