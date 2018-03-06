#ifndef _DELTAM_CR_HPP
#define _DELTAM_CR_HPP

#include <tranalisi.hpp>

#ifndef EXTERN_DELTAM_CR
 #define EXTERN_DELTAM_CR extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

EXTERN_DELTAM_CR djvec_t deltam_cr;

//! read or compute and write deltam_cr
void get_deltam_cr();

#undef EXTERN_DELTAM_CR
#undef INIT_TO

#endif
