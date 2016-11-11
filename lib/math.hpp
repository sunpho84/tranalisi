#ifndef _MATH_HPP
#define _MATH_HPP

#include <type_traits>

using namespace std;

//! helper for unary minus
template <class T> T uminus(const T &in)
{return -in;}

//! define the square of a double, float or integer
template <class T,class=typename enable_if<is_arithmetic<T>::value>::type> T sqr(T x)
{return x*x;}

#endif
