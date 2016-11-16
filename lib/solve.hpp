#ifndef _SOLVE_HPP
#define _SOLVE_HPP

#include <algorithm>
#include <cmath>
#include <tools.hpp>

#include <iostream>

using namespace std;

//! minimal difference
constexpr double epsilon=2*numeric_limits<double>::epsilon();

//! crash if the two points do not bracket a zero
inline void check_not_same_sign(double fa,double fb)
{if(same_sign(fa,fb)) CRASH("f(a) and f(b) do not have opposite sign: %lg %lg",fa,fb);}

//! solve using bisection
template <typename func_t> double bisect_solve(func_t fun,double a,double b)
{
  double fa=fun(a),fb=fun(b);
  check_not_same_sign(fa,fb);
  
  while(fabs(fa)>epsilon and fabs(a-b)>epsilon)
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
template <typename func_t> double Brent_solve(func_t fun,double a,double b)
{
  double fa=fun(a),fb=fun(b);
  check_not_same_sign(fa,fb);
  double d=epsilon,s=0;
  
  double c=a,fc=fa;
  bool mflag=true;
  int iter=0;
  
#ifdef DEBUG
  cout.precision(16);
  cout<<"Fa: "<<fa<<", Fb: "<<fb<<endl;
#endif
  
  while(fabs(fa)>epsilon and fabs(a-b)>epsilon*(fabs(a)+fabs(b)))
    {
      if(fa!=fc and fb!=fc)
	{
#ifdef DEBUG
	  cout<<iter+1<<" Inverse quadratic interpolation"<<endl;
#endif
	  s=a*fb*fc/((fa-fb)*(fa-fc))+
	    b*fc*fa/((fb-fc)*(fb-fa))+
	    c*fa*fb/((fc-fa)*(fc-fb));
	}
      else
	{
#ifdef DEBUG
	  cout<<iter+1<<" Secant Rule"<<endl;
#endif
	  s=b-fb*(b-a)/(fb-fa);
	}
      
      double tmp=(3*a+b)/4;
      double what_to_comp=mflag?fabs(b-c):fabs(c-d);
      if((!((s>tmp and s<b) or
	    (s<tmp and s>b))) or
	 fabs(s-b)>=what_to_comp/2 or
	  what_to_comp<epsilon)
	{
#ifdef DEBUG
	  cout<<iter+1<<" Bisection"<<endl;
#endif
	  s=(a+b)/2;
	  mflag=1;
	}
      else mflag = 0;
      
      double fs=fun(s);
      d=c;
      c=b;
      fc=fb;
      
      //replace the point with opposite sign
      if(same_sign(fb,fs))
	{
	  b=s;
	  fb=fs;
	}
      else
	{
	  a=s;
	  fa=fs;
	}
      
      //makes again a the larger
      if(fabs(fa)<fabs(fb))
	{
	  swap(a,b);
	  swap(fa,fb);
	}
      iter++;
      
      if(iter>100000) CRASH("Error");
#ifdef DEBUG
      cout<<" current interval: ["<<a<<", "<<b<<"], values: "<<fa<<", "<<fb<<endl;
#endif
    }
  
#ifdef DEBUG
  cout<<"Number of iterations: "<<iter<<endl;
#endif
  
  return b;
}

#endif
