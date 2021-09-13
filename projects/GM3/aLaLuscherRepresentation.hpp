#ifndef _ALALUSCHERREPRESENTATION_HPP
#define _ALALUSCHERREPRESENTATION_HPP

#include <array>
#include <complex>
#include <gslWrap.hpp>

#include <math.hpp>

#include <phiFunction.hpp>

template <typename T>
std::array<T,2> hHPrime(const T& mPi,
			const T& g2,
			const T& k,
			const T& x)
{
  const T m2Pi=
    sqr(mPi);
  
  const T f1=
    g2/(6*M_PI)*
    sqr(k)/(M_PI*x);
  
  const T f2=
    log((x+2*k)/(2*mPi));
  
  const T h=
    f1*2*k*f2;
  
  const T hPrime=
    f1*(1+(1+2*m2Pi/sqr(x))*x/k*f2);
  
  return
    {h,hPrime};
}

inline double cotDelta11(const double& mPi,
			 const double& mRho,
			 const double& g2,
			 const double& kN)
{
  const double m2Pi=
    sqr(mPi);
  
  const double m2Rho=
    sqr(mRho);
  
  const double omega=
    2*sqrt(m2Pi+sqr(kN));
  
  const double omega2=
    sqr(omega);
  
  const double k=
      sqrt(omega2/4-m2Pi);
  
  const double Gamma=
    g2*k*k*k/(6*M_PI*omega2);
  
  const auto hHPrimeOmega=
    hHPrime(mPi,g2,k,omega);
  
  const double& hOmega=
    hHPrimeOmega[0];

  CRASH("%lg %lg\n",hOmega,g2/(6*M_PI)*k*k*k/omega*2/M_PI*log((omega+2*k)/(2*mPi)));
  
  // const double kRho=
  //   sqrt(m2Rho/4-m2Pi);
 
  const auto hHPrimeRho=
    hHPrime(mPi,g2,k,mRho);
  
 const double& hRho=
   hHPrimeRho[0];
 //g2/(6*M_PI)*kRho*kRho*kRho/mRho*2/M_PI*log((mRho+2*kRho)/(2*mPi));
  
 const double& hPrimeRho=
   hHPrimeRho[1];
 //g2/(6*M_PI)*kRho*kRho/(M_PI*mRho)*(1+(1+2*m2Pi/m2Rho)*mRho/kRho*log((mRho+2*kRho)/(2*mPi)));
  
  const double cotDelta11=
    (m2Rho-omega2-hRho-(omega2-m2Rho)*hPrimeRho/(2*mRho)+hOmega)/(omega*Gamma);
  
  return
    cotDelta11;
}

inline vector<double> momentumFinder(const double& mPi,
				     const double& mRho,
				     const double& g2,
				     const double& L,
				     const int& n,
				     const string& path="")
{
  TanPhiCalculator tanPhiCalculator;
  
  auto f=
    [&mPi,&mRho,&g2,&L,&tanPhiCalculator](const double& kN)
    {
      const double cotD11=
	cotDelta11(mPi,mRho,g2,kN);
      
      const double tanD11=
	1.0/cotD11;
      
      const double tanPhi=
	tanPhiCalculator(kN*L/(2*M_PI));
      
      /// tan(a+b)=(tan(a)+tan(b))/(1-tan(a)*tan(b))
      const double tanD11PlusPhi=
	(tanD11+tanPhi)/(1-tanD11*tanPhi);
	
      return
	tanD11PlusPhi;
    };
  
  auto g=
    [&f](const double& x)
    {
      return
	1/f(x);
    };
  
  /// Positions of discontinuities
  vector<double> discontinuities;
  
  /// Positions of zeroes
  vector<double> zero;
  
  /// Starting value to be searched
  double x=
    0.1;
  
  // plotLhs.set_no_line();
  for(int i=1;i<=n;i++)
    {
      x=NewtonSolve(g,x+0.01,1e-8);
       discontinuities.push_back(x);
      
      x=NewtonSolve(f,x+0.01);
      zero.push_back(x);
    }
  
  if(path!="")
    {
      grace_file_t plot(path);
      
      for(double x=0.1;x<0.7;x+=0.001)
	plot.write_xy(x,atan(f(x)));
      
      plot.set_no_symbol();
      plot.new_data_set();
      
      for(auto& m : discontinuities)
	plot.write_xy(m,M_PI/2.0);
      
      plot.new_data_set();
      plot.set_no_line();
      
      for(auto& z : zero)
	plot.write_xy(z,0.0);
    }
  
  return
    zero;
}

#endif
