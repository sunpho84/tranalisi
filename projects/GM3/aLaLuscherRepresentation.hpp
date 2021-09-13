#ifndef _ALALUSCHERREPRESENTATION_HPP
#define _ALALUSCHERREPRESENTATION_HPP

#include <array>
#include <complex>
#include <gslWrap.hpp>

#include <math.hpp>

#include <phiFunction.hpp>

struct ALaLuscherRepresentation
{
  const double mPi;
  
  const double m2Pi;
  
  const double L;
  
  const double g2;
  
  const double mRho;
  
  const double m2Rho;
  
  const double kRho;
  
  const TanPhiCalculator tanPhiCalculator;
  
  ALaLuscherRepresentation(const double& mPi,
			   const double& L,
			   const double& g2,
			   const double& mRho) :
    mPi(mPi),
    m2Pi(sqr(mPi)),
    L(L),
    g2(g2),
    mRho(mRho),
    m2Rho(sqr(mRho)),
    kRho(sqrt(m2Rho/4-m2Pi))
  {
  }
  
  std::array<double,2> hHPrime(const double& k,
			       const double& x) const
  {
    const double f1=
    g2/(6*M_PI)*
      sqr(k)/(M_PI*x);
    
    const double f2=
      log((x+2*k)/(2*mPi));
    
    const double h=
      f1*2*k*f2;
    
    const double hPrime=
      f1*(1+(1+2*m2Pi/sqr(x))*x/k*f2);
    
    return
      {h,hPrime};
  }
  
  double cotDelta11(const double& kN) const
  {
    const double omega=
      2*sqrt(m2Pi+sqr(kN));
    
    const double omega2=
      sqr(omega);
    
    const double k=
      sqrt(omega2/4-m2Pi);
    
    const double Gamma=
      g2*k*k*k/(6*M_PI*omega2);
    
    const auto hHPrimeOmega=
      hHPrime(k,omega);
    
    const double& hOmega=
      hHPrimeOmega[0];
    
    const auto hHPrimeRho=
      hHPrime(kRho,mRho);
    
    const double& hRho=
      hHPrimeRho[0];
    
    const double& hPrimeRho=
      hHPrimeRho[1];
    
    const double cotDelta11=
      (m2Rho-omega2-hRho-(omega2-m2Rho)*hPrimeRho/(2*mRho)+hOmega)/(omega*Gamma);
    
    return
      cotDelta11;
  }
  
  double levelsFinderLhs(const double& kN) const
  {
    const double cotD11=
      cotDelta11(kN);
    
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
  
  vector<double> findLevels(const int& n,
			    const string& path="")
    const
  {
    const auto zeroFinder=
      [this](const double& kN)
      {
	return
	  levelsFinderLhs(kN);
      };
    
    auto discontinuityFinder=
      [this](const double& kN)
      {
	return
	  1/levelsFinderLhs(kN);
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
	x=NewtonSolve(discontinuityFinder,x+0.01,1e-8);
	discontinuities.push_back(x);
	
	x=NewtonSolve(zeroFinder,x+0.01);
	zero.push_back(x);
      }
    
    if(path!="")
      {
	grace_file_t plot(path);
	
	for(double x=0.1;x<0.7;x+=0.001)
	  plot.write_xy(x,atan(zeroFinder(x)));
	
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
};

#endif
