#include <tranalisi.hpp>

#include <luscherZeta.hpp>

#include <iostream>

using namespace std;

const int lambdaMax=30;
static DistanceMultiplicity distanceMultiplicity(lambdaMax);

const LuscherZetaCalculator luscherZetaCalculator(lambdaMax);

double q2=0.3;
int lambda=4;


double aaaa(const double& t)
{
  /// Previous iteration result
  double oldRes=0.0;
  
  /// Result
  double res=0.0;
  
  /// Point multiplicity
  int mult;
  
  /// Iteration distance
  int distance2=lambda*lambda;
  do
    {
      mult=distanceMultiplicity(distance2);
      
      if(mult)
	{
	  /// Current contribution
	  const double contr=
	    mult*exp((q2-distance2)*t);
	  
	  oldRes=res;
	  res=oldRes+contr;
	}
      
      distance2++;
    }
  while(distance2<lambdaMax*lambdaMax and (fabs(oldRes!=res or mult==0)));
  
  const double tot=
    res*pow(2*M_PI,-3.0)/sqrt(4*M_PI)-pow(4*M_PI,-2)*pow(t,-1.5);
  
  return tot;
}

double bbbb(const double& t)
{
  
  /// Previous iteration result
  double oldRes=0.0;
  
  /// Result
  double res=0.0;
  
  /// Point multiplicity
  int mult;
  
  /// Iteration distance
  int distance2=0;
  do
    {
      mult=distanceMultiplicity(distance2);
      
      if(mult)
	{
	  /// Current contribution
	  const double contr=
	    mult*exp(q2*t-M_PI*M_PI*distance2/t);
	  
	  cout<<"small: "<<distance2<<" "<<contr<<endl;
	  oldRes=res;
	  res=oldRes+contr;
	}
      
      distance2++;
    }
  while(distance2<lambdaMax*lambdaMax and (fabs(oldRes!=res or mult==0)));
  
  res*=pow(4*M_PI*t,-1.5);
  
  double sub=0.0;
  for(int distance2=0;distance2<lambda*lambda;distance2++)
    {
      const int mult=distanceMultiplicity(distance2);
      
      /// Current contribution
      const double contr=
	exp(-t*(distance2-q2));
      
      cout<<"large: "<<distance2<<" "<<contr<<endl;
      
      sub+=contr*mult;
    }
  sub*=pow(2*M_PI,-3.0);
  
  const double tot=(res-sub)/sqrt(4*M_PI)-pow(4*M_PI,-2)*pow(t,-1.5);
  cout<<" resa: "<<res/sqrt(4*M_PI)<<" "<<sub/sqrt(4*M_PI)<<endl<<endl;
  
  return tot;
}

//! kernel of eq.10
double lzc(double t,
	   void *params)
{
  return ((t>1)?aaaa:bbbb)(t);
}

//! compute tilde for double
double luscherZeta()
{
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  //! function structure
  gsl_function f;
  f.function=lzc;
  f.params=nullptr;
  
  //integrate
  double // part1=0.0,
    part2;
  double abserr;
  const double epsabs=0,epsrel=1e-8;
  //gsl_integration_qags(&f,0.0,1.0,epsabs,epsrel,workspace_size,workspace,&part1,&abserr);
  gsl_integration_qagiu(&f,0.0,epsabs,epsrel,workspace_size,workspace,&part2,&abserr);
  gsl_integration_workspace_free(workspace);
  
  /////////////////////////////////////////////////////////////////
  
  double add=0.0;
  for(int distance2=0;distance2<lambda*lambda;distance2++)
    {
      const int mult=distanceMultiplicity(distance2);
      
      /// Current contribution
      const double contr=
	1.0/(distance2-q2);
      
      add+=contr*mult;
    }
  add/=sqrt(4.0*M_PI);
  
  return add+part2*pow(2.0*M_PI,3.0);
}

int main()
{
  cout.precision(16);
  
  q2=3.0;
  lambda=3.0;
  const double o=lzc(0.4,nullptr);
  const double n=luscherZetaCalculator.heatKernelCalculator(0.4,q2);
  cout<<o<<" "<<n<<endl;
  return 0;
  
  
  double oldPhi=0.0;
  int iCorr=0;
  for(q2=0.081;q2<15.01;q2+=0.01)
    {
      lambda=floor(sqrt(q2))+2;
      const double l=luscherZeta();
      const double tanPhi=-pow(M_PI,1.5)*sqrt(q2)/l;
      const double phi=atan(tanPhi);//+floor(q2)*M_PI;
      if(oldPhi and fabs(phi-oldPhi)>0.1) iCorr++;
      //const int iCorr=round((phi-M_PI*q2)/M_PI);
      const double corr=phi+iCorr*M_PI;
      
      cout<<q2<<" "<<corr/(M_PI*q2)<<// phi/<<
	  endl;
      
      oldPhi=phi;
    }
  return 0;
  // 	lzc(2.0,nullptr);
  
  
  // return 0;
  for(double t=1-1e-6;t<1+1e-6;t*=1+1e-8)
    {
      double a=// (t<M_PI)?
	// luscherZetaCalculator.integrand(t,q2,lambda):
	lzc(t,nullptr);
      
      
      cout<<t<<" "<<a<<endl;
    }
    //cout<<t<<" "<<luscherZetaCalculator.heatKernelCalculator.heatKernelForLargeArg(t)<<endl;
  
  return 0;
}
