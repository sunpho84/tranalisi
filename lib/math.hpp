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
constexpr T sqr(const T &x)
{return x*x;}

//! define the cube
template <class T>
T cube(const T &x)
{return x*x*x;}

//! returns the normalized gaussian
template <class T>
T gauss_distr(double x,double x0,const T &sigma)
{return 1/(sqrt(2*M_PI)*sigma)*exp(-sqr((x-x0)/sigma)/2.0);}

//////////////////////////////// static powers

namespace Impl
{
  //! internal implementation of power
  template <int N>
  struct PowN
  {
    template <class T>
    static constexpr T get(T x)
    {
      return
	PowN<N/2>::get(x)*
	PowN<N-N/2>::get(x);
    }
  };
  
  //! special case, N=0, terminator for the recursion
  template <>
  struct PowN<0>
  {
    template <class T>
    static constexpr T get(T x)
    {return 1;}
  };
  
  //! special case, N=0, terminator for the recursion
  template <>
  struct PowN<1>
  {
    template <class T>
    static constexpr T get(T x)
    {return x;}
  };
}

//! static version of power
template <int N,class T>
constexpr T Pow(T x)
{return Impl::PowN<N>::get(x);}

//! static Square
template <class T>
constexpr T Sqr(T x)
{return Pow<2>(x);}

//! static Cube
template <class T>
constexpr T Cube(T x)
{return Pow<3>(x);}

/// Factorial of n
constexpr int fact(const int& n)
{
  int r=1;
  
  for(int i=2;i<=n;i++)
    r*=i;
  
  return r;
}

/// Gaussian centered in x0 with standard deviation sigma, normalized to 1
template <typename Real>
Real gauss(const Real& x0,
	   const Real& sigma,
	   const Real& x)
{
  return exp(-sqr((x-x0)/sigma)/2)/(sqrt(2*M_PI)*sigma);
}

#endif
