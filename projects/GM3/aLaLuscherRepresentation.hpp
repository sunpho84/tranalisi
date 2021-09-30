#ifndef _ALALUSCHERREPRESENTATION_HPP
#define _ALALUSCHERREPRESENTATION_HPP

#include <array>
#include <complex>
#include <gslWrap.hpp>

#include <math.hpp>

#include "params.hpp"
#include <phiFunction.hpp>

/// Representation
struct ALaLuscherRepresentation
{
  /// Factor of the correction
  static constexpr double isoVectorCorrection=
    10.0/9;
  
  /// Single state contribution
  struct Contr
  {
    /// Weight
    double weight;
    
    /// Energy of the state
    double energy;
    
    /// Sets the parameters
    void set(const double& _weight,
	     const double& _energy)
    {
      weight=_weight;
      energy=_energy;
    }
    
    /// Evaluates for time t, including the isovector correction
    double operator()(const double& t) const
    {
      return
	isoVectorCorrection*weight*exp(-t*energy);
    }
  };
  
  /// Coefficient of the expansion
  vector<Contr> coeffs;
  
  /// Evaluates for time t
  double operator()(const double& t) const
  {
    double res=
      0.0;
    
    for(const auto& p : coeffs)
      res+=p(t);
    
    return
      res;
  }
  
  /// Set the state
  void setState(const int& i,
		const double& coeff,
		const double& energy)
  {
    coeffs[i].set(coeff,energy);
  }
  
  /// Initialize with n states
  ALaLuscherRepresentation(const int& n=0) :
    coeffs(n)
  {
  }
};

struct ALaLuscherRepresentationCalculator
{
  const double mPi;
  
  const double m2Pi;
  
  const double L;
  
  const double g2;
  
  const double mRho;
  
  const double m2Rho;
  
  const double kRho;
  
  const double App0;
  
  double Gamma(const double& omega2,
	       const double& k) const
  {
    return
      g2*k*k*k/(6*M_PI*omega2);
  }
  
  double getApp0() const
  {
    const auto [h,hPrime]=
      hHPrime(kRho,mRho);
    
    return
      h-mRho/2*hPrime+g2*m2Pi/(6*M_PI*M_PI);
  }
  
  complex<double> App(const double& omega) const
  {
    const double omega2=
      omega*omega;
    
    const double k=
      sqrt(omega2/4-m2Pi);
    
    const auto [hRho,hPrimeRho]=
      hHPrime(kRho,mRho);
    
    const auto [h,hPrime]=
      hHPrime(k,omega);
    
    return
      {hRho+(omega2-m2Rho)*hPrimeRho/(2*mRho)-h,omega*Gamma(omega2,k)};
  }
  
  complex<double> F(const double& omega) const
  {
    const double num=
      m2Rho-App0;
    
    const complex<double> den=
      m2Rho-omega*omega-App(omega);
    
    return
      num/den;
  }
  
  /// Computes omega given kN
  double omegaOf(const double& kN) const
  {
    return
      2*sqrt(m2Pi+sqr(kN));
  }
  
  /// Computes kN given omega
  double kNOf(const double& omega) const
  {
    return
      sqrt(omega*omega/4-m2Pi);
  }
  
  /// Computes the weight
  double contrWeight(const double& kN) const
  {
    /// Omega entering the definition
      const double omega=
	omegaOf(kN);
    
    /// d11'
    const double d11deriv=
      gslDeriv([this](const double& x)
      {
	return
	  delta11(x);
      },
	kN,
	1e-3);
    
    /// Argument of the derivative of thi
    const double arg=
      kN*L/(2*M_PI);
    
    /// Computes phi'
    const double phiDeriv=
      hashedPhiAndDerivCalculator.phiDeriv(arg);
    
    /// Norm2 of the form factor
    const double F2Norm=
      norm(F(omega));
    
    return
      2*pow(kN,5)/(3*M_PI*sqr(omega))*F2Norm/
      (kN*d11deriv+arg*phiDeriv);
  }
  
  /// Constructor
  ALaLuscherRepresentationCalculator(const double& mPi,
				     const double& L,
				     const double& g2,
				     const double& mRho) :
    mPi(mPi),
    m2Pi(sqr(mPi)),
    L(L),
    g2(g2),
    mRho(mRho),
    m2Rho(sqr(mRho)),
    kRho(sqrt(m2Rho/4-m2Pi)),
    App0(getApp0())
  {
  }
  
  /// Computes h and d'
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
  
  /// cot(delta11) function
  double cotDelta11(const double& kN) const
  {
    const double omega=
      omegaOf(kN);
    
    const double omega2=
      sqr(omega);
    
    const double k=
      sqrt(omega2/4-m2Pi);
    
    const auto [hOmega,hPrimeOmega]=
      hHPrime(k,omega);
    
    const auto [hRho,hPrimeRho]=
      hHPrime(kRho,mRho);
    
    const double cotDelta11=
      (m2Rho-omega2-hRho-(omega2-m2Rho)*hPrimeRho/(2*mRho)+hOmega)/(omega*Gamma(omega2,k));
    
    return
      cotDelta11;
  }
  
  /// Delta11 function
  double delta11(const double& kN) const
  {
    double delta11=
      atan(1/cotDelta11(kN));
    
    if(delta11<0)
      delta11+=M_PI;
    
    return
      delta11;
  }
  
  /// Derivative of the delta11 function
  double derDelta11(const double& kN) const
  {
    return
      gslDeriv([this](const double& x)
      {
	return
	  delta11(x);
      },
	kN,
	1e-3);
  }
  
  /// Computes the l.h.s of the quantization condition
  double levelsFinderLhs(const double& kN) const
  {
    const double cotD11=
      cotDelta11(kN);
    
    const double tanD11=
      1.0/cotD11;
    
    const double tanPhi=
      hashedPhiAndDerivCalculator.tanPhi(kN*L/(2*M_PI));
    
    /// tan(a+b)=(tan(a)+tan(b))/(1-tan(a)*tan(b))
    const double tanD11PlusPhi=
      (tanD11+tanPhi)/(1-tanD11*tanPhi);
    
    return
      tanD11PlusPhi;
  };
  
  /// Find n levels
  vector<double> findLevels(const int& n,
			    const string& path="")
    const
  {
    //cout<<"  searching "<<n<<" levels"<<endl;
    
    const auto zeroFinder=
      [this](const double& kN)
      {
	//cout<<"     searching for zero near "<<kN<<endl;
	
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
    
    /// Minimal value to search
    const double xMin=
      1e-4;
    
    /// Starting value to be searched
    double x=
      xMin;
    
    for(int i=1;i<=n;i++)
      {
	//cout<<"   "<<i<<" initial "<<x<<endl;
	x=NewtonSolve(discontinuityFinder,x+1e-4);
	discontinuities.push_back(x);
	//cout<<"   "<<i<<" found disco "<<x<<endl;
	
	x=NewtonSolve(zeroFinder,x+1e-4);
	zero.push_back(x);
	//cout<<"   "<<i<<" found zero "<<x<<endl;
      }
    
    /// Last point to plot
    const double xMax=
      x*1.1;
    
    if(path!="")
      {
	grace_file_t plot(path);
	
	for(double x=xMin;x<xMax;x+=0.001)
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
  
  /// Gets the parameters of the states
  ALaLuscherRepresentation getLevelsPars(const int& n,
					 const string& path="") const
  {
    /// Gets all levels momentum
    const vector<double> kNs=
      findLevels(n,path);
    
    /// Result
    ALaLuscherRepresentation res(n);
    
    for(int i=0;i<n;i++)
      {
	//cout<<"  adding solution for level "<<i<<endl;
	
	/// Single momentum
	const double& kN=
	  kNs[i];
	
	/// Weight of the state
	const double weight=
	  contrWeight(kN);
	
	const double omega=
	  omegaOf(kN);
	
	res.setState(i,weight,omega);
      }
    
    return
      res;
  }
  
  ALaLuscherRepresentation operator()(const int& n,
				      const string& path="") const
  {
    return
      getLevelsPars(n,path);
  }
};

/// Caches the calculation of the Luscher representation
struct ALaLuscherRepresentationCached
{
  /// Cached value of mPi
  mutable double cachedMPi;
  
  /// Number of levels to be searched
  const int n;
  
  /// Storage of the interpolation parameters
  mutable map<pair<double,double>,ALaLuscherRepresentation> cachedPars;
  
  /// Constructor
  ALaLuscherRepresentationCached(const int& n) :
    cachedMPi(0.0),
    n(n)
  {
  }
  
  /// Computes the representation
  const ALaLuscherRepresentation& operator()(const double& mPi,
					     const double& L,
					     const double& mRho,
					     const double& g2) const
  {
    if(cachedMPi!=mPi)
      {
	//cout<<"Changing mPi from "<<cachedMPi<<" to "<<mPi<<endl;
	
	cachedPars.clear();
	cachedMPi=mPi;
      }
    
    /// Search key
    const pair<double,double> key{mRho,g2};
    
    /// Result of the research
    auto find=
      cachedPars.find(key);
    
    if(find==cachedPars.end())
      {
	//cout<<"mPi "<<mPi<<" , L "<<L<<" , g2 "<<g2<<" , mRho "<<mRho<<endl;
	
	//cout<<" computing for values: mPi="<<mPi<<" , g2="<<g2<<" , mRho="<<mRho<<endl;
	const ALaLuscherRepresentationCalculator interacting(mPi,L,g2,mRho);
	const ALaLuscherRepresentation pars=
	  interacting.getLevelsPars(n);
	
	const auto insertRes=
	  cachedPars.insert({key,pars});
	
	find=insertRes.first;
      }
    
    return
      find->second;
  }
};

static const double qMax=
    4.3;
  
const int nIntervals=
  1000;

struct ALaLuscherRepresentationInfVol
{
  gsl_integration_workspace *workspace;
  
  /// Absolute error
  static constexpr double epsAbs=0;
  
  /// Size of the workspace
  static constexpr int workspaceSize=
    1000;
  
  const ALaLuscherRepresentationCalculator aLaLusch;
  
  const double mPi;
  
  double integrand(const double& omega,
		   const double& t) const
  {
    return
      10.0/9*(1.0/(48.0*sqr(M_PI)))*sqr(omega)*pow(1.0-sqr(2.0*mPi/omega),1.5)*exp(-omega*t)*norm(aLaLusch.F(omega));
  };
  
  double operator()(const double& t) const
  {
    /// Result
    double res;
    
    /// Absolute error
    double absErr;
    
    const double& eps=1e-8;
    
    gsl_integration_qagiu(GslFunction([this,t](const double& omega)
      {
	return
	  integrand(omega,t);
      }).getWrapper(),2*mPi,epsAbs,eps,workspaceSize,workspace,&res,&absErr);
    
    return
      res;
  }
  
  ALaLuscherRepresentationInfVol(const double& mPi,
				 const double& g2,
				 const double& mRho) :
    workspace(gsl_integration_workspace_alloc(workspaceSize)),
    aLaLusch(mPi,1e30,g2,mRho),
    mPi(mPi)
  {
  }
  
  ALaLuscherRepresentationInfVol(const ALaLuscherRepresentationInfVol& oth) :
    workspace(gsl_integration_workspace_alloc(workspaceSize)),
    aLaLusch(oth.aLaLusch),
    mPi(oth.mPi)
  {
  }
  
  ALaLuscherRepresentationInfVol(ALaLuscherRepresentationInfVol&& oth) :
    aLaLusch(std::move(oth.aLaLusch)),
    mPi(std::move(oth.mPi))
  {
    workspace=
      oth.workspace;
    
    oth.workspace=
      nullptr;
  }
  
  ~ALaLuscherRepresentationInfVol()
  {
    if(workspace)
      gsl_integration_workspace_free(workspace);
  }
};

#endif
