#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <solve.hpp>

double Brent_solve(function<double(double)> fun,double a,double b,double rel_tol)
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
  
  while(fabs(fa)>epsilon and fabs(fb)>epsilon and !are_close(a,b,rel_tol))
    {
      if(same_sign(fa,fb)) CRASH("Something went wrong");
      
      if(!are_close(fa,fc,rel_tol) and !are_close(fb,fc,rel_tol))
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
