#ifndef _SOLVE_HPP
#define _SOLVE_HPP

#include <iostream>

using namespace std;

template <typename func_t> double get_solve_err(func_t fun,double guess,double err)
{
  //enlarge until opposite sign bracketed
  while(fun(guess-err)*fun(guess+err)>=0) err*=2;
  
  return err;
}

//! solve using bisection
template <typename func_t> double bisect_solve(func_t fun,double guess=1,double err=1e-10)
{
  //get original guess
  double x=guess;
  err=get_solve_err(fun,guess,err);
  
  //bisect
  double xl=x-err,xr=x+err,yl=fun(xl),yr=fun(xr);
  do
    {
      x=(xl+xr)/2;
      double y=fun(x);
      
      if((yl>0 && y>0) || (yl<0 && y<0))
	{
	  xl=x;
	  yl=y;
	}
      else
	{
	  xr=x;
	  yr=y;
	}
      cout<<xl<<" "<<x<<" "<<xr<<endl;
      }
  while((xl!=xr) && ((yl>0 && yr<0) || (yl<0 && yr>0)));
  
  return x;
}

#endif
