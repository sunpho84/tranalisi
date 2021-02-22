#ifndef _GLOBAL_HPP
#define _GLOBAL_HPP

#include <stdlib.h>

#ifndef EXTERN_GLOBAL
 #define EXTERN_GLOBAL extern
#endif

constexpr int q1=2,q2=-1;

/// Number of X
EXTERN_GLOBAL size_t nX;

/// Subtract or not the zero momentum
EXTERN_GLOBAL size_t subtractZeroV;

#undef EXTERN_GLOBAL

#endif
