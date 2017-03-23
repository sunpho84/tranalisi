#ifndef _SOLVE_HPP
#define _SOLVE_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <tools.hpp>

#include <iostream>

using namespace std;

//! minimal difference
constexpr double epsilon=2*numeric_limits<double>::epsilon();

//! crash if the two points do not bracket a zero
inline void check_not_same_sign(double fa,double fb)
{if(same_sign(fa,fb)) CRASH("f(a) and f(b) do not have opposite sign: %lg %lg",fa,fb);}

//! check if |a-b|<=|a+b|
inline bool are_close(double a,double b,double rel_tol)
{return fabs(a-b)<=rel_tol*(fabs(a+b));}

//! solve using bisection
template <typename func_t> double bisect_solve(func_t fun,double a,double b,double rel_tol=epsilon*10)
{
  double fa=fun(a),fb=fun(b);
  check_not_same_sign(fa,fb);
  
  while(fabs(fa)>epsilon and !are_close(a,b,rel_tol))
    {
      //compute midpoint
      double x=(a+b)/2;
      double f=fun(x);
      cout<<a<<" "<<x<<" "<<b<<endl;
      
      if(same_sign(fa,f))
	{
	  a=x;
	  fa=f;
	}
      else
	{
	  b=x;
	  fb=f;
	}
    }
  
  return a;
}

//! solve using Brent method
double Brent_solve(function<double(double)> fun,double a,double b,double rel_tol=epsilon*10);

//! Brent solve for type
template <class TS,class=enable_if_t<has_method_size<TS>::value>> TS Brent_solve(function<double(double,size_t)> fun,double a,double b,double rel_tol=epsilon*10)
{
  TS out;
  for(size_t iel=0;iel<out.size();iel++) out[iel]=Brent_solve(bind(fun,placeholders::_1,iel),a,b,rel_tol);
  return out;
}

#endif
