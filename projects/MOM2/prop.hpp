#ifndef _PROP_HPP
#define _PROP_HPP

#include <geometry.hpp>

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
 #define INIT_PROP_TO(...)
#else
 #define INIT_PROP_TO(...) = __VA_ARGS__
#endif

using vprop_t=vector<qprop_t>; //! vector of propagators
using vjqprop_t=vector<jqprop_t>; //! vector of jackkniffed props

#undef EXTERN_PROP
#undef INIT_PROP_TO

#endif
