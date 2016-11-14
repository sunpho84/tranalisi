#ifndef _SOLVE_HPP
#define _SOLVE_HPP

#include <algorithm>
#include <cmath>

#include <iostream>

using namespace std;

template <typename func_t> double get_solve_err(func_t fun,double guess,double err)
{
  //enlarge until opposite sign bracketed
  while(fun(guess-err)*fun(guess+err)>=0) err*=2;
  
  return err;
}

//! get the guess with fb[1]>fb[0]
template <typename func_t> double get_guess(double xb[2],double fb[2],func_t fun,double guess,double err)
{
  double x=guess;
  err=get_solve_err(fun,guess,err);
  
  //bisect
  xb[0]=x-err;
  xb[1]=x+err;
  fb[0]=fun(xb[0]);
  fb[1]=fun(xb[1]);
  
  //swap to have fb[1]>0
  if(fb[1]<fb[0])
    {
      swap(fb[0],fb[1]);
      swap(xb[0],xb[1]);
    }
  
  return x;
}

//! solve using bisection
template <typename func_t> double bisect_solve(func_t fun,double guess=1,double err=1e-10)
{
  double xb[2],fb[2];
  double x=get_guess(xb,fb,fun,guess,err);
  
  do
    {
      //compute midpoint
      x=(xb[0]+xb[1])/2;
      double f=fun(x);
      cout<<xb[1]<<" "<<x<<" "<<xb[0]<<endl;
      //change into fb[1] if f>0
      bool w=(f>0);
      xb[w]=x;
      fb[w]=f;
    }
  while(xb[0]!=xb[1] && fb[1]>0 && fb[0]<0);
  
  return x;
}

//! solve using Brent method
template <typename func_t> double brent_solve(func_t fun,double guess=1,double err=1e-10)
{
  double x=guess;
  err=get_solve_err(fun,guess,err);
  double a=x-err,b=x+err;
  double fa=fun(a),fb=fun(b);
  double s,fs;
  const double tol=1e-16;
  
  do
    {
      double c=(a+b)/2;
      double fc=fun(c);
      if(fabs(fa-fc)>tol && fabs(fb-fc)>tol)
	{
	  cout<<"Inverse quadratic interpolation"<<endl;
	  s=a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb));
      }
      else
	{
	  cout<<"Secant rule"<<endl;
	  s=b-fb*(b-a)/(fb-fa);
	}
      fs=fun(s);
      if(c>s)
	{
	  swap(s,c);
	  swap(fs,fc);
	}
      if(fc*fs<0)
	{
	  a=s;
	  fa=fs;
	  b=c;
	  fb=fc;
	}
      else
	if(fs*fb<0)
	  {
	    a=c;
	    fa=fc;
	  }
	else
	  {
	    b=s;
	    fb=fs;
	  }
      
      cout<<a<<" "<<s<<" "<<b<<" "<<fs<<" "<<tol<<" "<<(fs<tol)<<" "<<(fabs(b-a)<tol)<<endl;
    }
  while(fs>tol||fabs(b-a)>tol);	// Convergence conditions
  
  return s;
}

#endif
