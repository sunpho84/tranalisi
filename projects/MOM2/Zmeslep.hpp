#ifndef _ZMESLEP_HPP
#define _ZMESLEP_HPP

#include <Dirac.hpp>

#include <prop.hpp>

#ifndef EXTERN_MESLEP
 #define EXTERN_MESLEP extern
 #define INIT_MESLEP_TO(...)
#else
 #define INIT_MESLEP_TO(...) __VA_ARGS__
#endif

#undef EXTERN_MESLEP
#undef INIT_MESLEP_TO

#endif
