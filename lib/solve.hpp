#ifndef _SOLVE_HPP
#define _SOLVE_HPP

#include <algorithm>
#include <cmath>
#include <tools.hpp>

#include <iostream>

using namespace std;

constexpr double epsilon=numeric_limits<double>::epsilon();

template <typename func_t> double get_bracket(func_t fun,double guess)
{
  //enlarge until opposite sign bracketed
  double err=epsilon;
  while(fun(guess-err)*fun(guess+err)>=0) err*=2;
  
  return err;
}

//! get the guess with fb[1]>fb[0]
template <typename func_t> double get_guess(double xb[2],double fb[2],func_t fun,double guess,double err)
{
  double x=guess;
  err=get_bracket(fun,guess);
  
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
template <typename func_t> double Brent_solve(func_t fun,double guess=1)
{
  double x=guess,err=get_bracket(fun,guess);
  double a=x-err,b=x+err,d=epsilon;
  double fa=fun(a),fb=fun(b);
  double s=0;
  
  if(same_sign(fa,fb)) CRASH("f(a) and f(b) do not have opposite sign: %lg %lg",fa,fb);
  
  if(fabs(fa)<fabs(fb))
    {
      swap(a,b);
      swap(fa,fb);
    }
  
  double c=a,fc=fa;
  bool mflag=true;
  int iter=0;
   
   while(!(fb==0) and (fabs(a-b)>epsilon))
     {
       if(fa!=fc and fb!=fc)
	 //{
	 //cout<<iter+1<<" Inverse quadratic interpolation"<<endl;
	 s=a*fb*fc/((fa-fb)*(fa-fc))+
	   b*fc*fa/((fb-fc)*(fb-fa))+
	   c*fa*fb/((fc-fa)*(fc-fb));
       //}
       else
	 //{
	 //cout<<i+1<<" Secant Rule"<<endl;
	 s=b-fb*(b-a)/(fb-fa);
       //}
       
       double tmp=(3*a+b)/4;
       double what_to_comp=mflag?fabs(b-c):fabs(c-d);
       if((!((s>tmp and s<b) or
	     (s<tmp and s>b))) or
	  fabs(s-b)>=what_to_comp/2 or
	  what_to_comp<epsilon)
	 {
	   //cout<<iter+1<<" Bisection"<<endl;
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
       
       if(fabs(fa)<fabs(fb))
	 {
	   swap(a,b);
	   swap(fa,fb);
	 }
       iter++;
       
       if(iter>100000) CRASH("Error");
       //cout<<s<<endl;
     }
   
   //cout<<"Number of iterations: "<<iter<<endl;
   
   return b;
}

#endif
