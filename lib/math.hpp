#ifndef _MATH_HPP
#define _MATH_HPP

#include <cmath>

#include <type_traits>

using namespace std;

//! helper for unary minus
template <class T>
T uminus(const T &in)
{return -in;}

//! define the square of a double, float or integer
template <class T>
T sqr(const T &x)
{return x*x;}

//! define the cube
template <class T>
T cube(const T &x)
{return x*x*x;}

//! returns the normalized gaussian
template <class T>
T gauss_distr(double x,double x0,const T &sigma)
{return 1/(sqrt(2*M_PI)*sigma)*exp(-sqr((x-x0)/sigma)/2.0);}

#endif
