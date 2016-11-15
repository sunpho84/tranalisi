#ifndef _SOLVE_HPP
#define _SOLVE_HPP

#include <algorithm>
#include <cmath>
#include <tools.hpp>

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
template <typename func_t> double brent_solve2(func_t fun,double guess=1,double err=1e-10)
{
  double x=guess;
  err=get_solve_err(fun,guess,err);
  double a=x-err,b=x+err;
  a=-4;
  b=4./3;
  //double fa=fun(a),fb=fun(b);
  
  if(fun(a)*fun(b)>=0) CRASH("Something went wrong");
  if(fabs(fun(a))<fabs(fun(b)))
    {
      swap(a,b);
      //swap(fun(a),fun(b));
    }
  double c=a;//,fc=fun(a);
  double d=1e300;
  double s=0;//,fs=1e300;
  bool mflag=true;
  double tol=numeric_limits<double>::epsilon();
  int iter=0;
  while(fun(b)!=0 && fun(s)!=0 && b!=a)
    {
      if(fun(a)!=fun(c) and fun(b)!=fun(c))
	{
	  cout<<iter+1<<" Inverse quadratic interpolation"<<endl;
	  s=a*fun(b)*fun(c)/((fun(a)-fun(b))*(fun(a)-fun(c)))+
	    b*fun(c)*fun(a)/((fun(b)-fun(c))*(fun(b)-fun(a)))+
	    c*fun(a)*fun(b)/((fun(c)-fun(a))*(fun(c)-fun(b)));
	}
      else
	{
	  cout<<iter+1<<" Secant method"<<endl;
	  s=b-fun(b)*(b-a)/(fun(b)-fun(a));
	}
      
      double tmp=(3*a+b)*0.25;
      if(!((tmp<s and s<b) or
	   (tmp>s and s>b)) or
	 ( mflag and fabs(s-b)>=fabs(b-c)/2) or
	 (!mflag and fabs(s-b)>=fabs(c-d)/2) or
	 ( mflag and fabs(b-c)<tol) or
	 (!mflag and fabs(c-d)<tol))
	{
	  cout<<iter+1<<" Bisection method"<<endl;
	  s=(a+b)/2;
	  mflag=true;
	}
      else mflag=false;
      
      //fs=fun(s);
      d=c;
      c=b;
      //fc=fb;
      
      if(fun(a)*fun(s)<0)
	{
	  b=s;
	  cout<<" changing b to s"<<endl;
	  //fb=fs;
	}
      else
	{
	  a=s;
	  cout<<" changing a to s"<<endl;
	  //fa=fs;
	}
      
      if(fabs(fun(a))<fabs(fun(b)))
	{
	  swap(a,b);
	  cout<<" swapping a,b"<<endl;
	  //swap(fun(a),fun(b));
	}
      
      iter++;
      cout<<s<<endl;
    }
  
  return s;
}

//! solve using Brent method
template <typename func_t> double brent_solve3(func_t fun,double guess=1,double err=1e-10)
{
  double x=guess;
  err=get_solve_err(fun,guess,err);
  double a=x-err,b=x+err;
  a=-4;
  b=4./3;
  double c = 0;
   double d = 1e-16;
   
   double fa = fun(a);
   double fb = fun(b);
   
   double fc = 0;
   double s = 0;
   double fs = 0;
   
   // if f(a) f(b) >= 0 then error-exit
   if (fa * fb >= 0)
   {
      if (fa < fb)
        return a;
      else
        return b;
   }
   
   // if |f(a)| < |f(b)| then swap (a,b) end if
   if (fabs(fa) < fabs(fb))
   { double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp; }
   
   c = a;
   fc = fa;
   int mflag = 1;
   int i = 0;
   
   while (!(fb==0) && (fabs(a-b) > 1e-16))
   {
      if ((fa != fc) && (fb != fc))
	{
	  cout<<i+1<<" Inverse quadratic interpolation"<<endl;
	  s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) /
	    (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
	}
      else
	{
	  cout<<i+1<<" Secant Rule"<<endl;
	  s = b - fb * (b - a) / (fb - fa);
	}
      
      double tmp2 = (3 * a + b) / 4;
      if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) ||
          (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) ||
          (!mflag && (fabs(s - b) >= (fabs(c - d) / 2))))
      {
	cout<<i+1<<" Bisection"<<endl;
         s = (a + b) / 2;
         mflag = 1;
      }
      else
      {
         if ((mflag && (fabs(b - c) < 1e-16)) ||
             (!mflag && (fabs(c - d) < 1e-16)))
         {
	cout<<i+1<<" Bisection"<<endl;
            s = (a + b) / 2;
            mflag = 1;
         }
         else
           mflag = 0;
      }
      fs = fun(s);
      d = c;
      c = b;
      fc = fb;
      if (fa * fs < 0) { b = s; fb = fs; }
      else { a = s; fa = fs; }
      
      // if |f(a)| < |f(b)| then swap (a,b) end if
      if (fabs(fa) < fabs(fb))
      { double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp; }
      i++;
      if (i > 100000)
      {
         printf("Error is %f \n", fb);
         break;
      }
      cout<<s<<endl;
   }
   
   printf("Number of iterations : %d\n",i);
   return b;
}

#endif
