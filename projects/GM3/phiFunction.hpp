#ifndef _PHIFUNCTION_HPP
#define _PHIFUNCTION_HPP

#include <gslWrap.hpp>
#include <hashedFunction.hpp>
#include <fit.hpp>

#include <luscherZeta.hpp>

/// Computes the r.h.s of eq.A.3 of Nuclear Physics B364 (1991) 237-251
struct TanPhiCalculator
{
  /// Computes the Zeta function
  const LuscherZetaCalculator luscherZetaCalculator;
  
  /// Approximation according to eq. A.4
  double approximation(const double& q2) const
  {
    /// Equation A.5
    const double s=
      q2*(8.91363-
	  q2*(16.5323-
	      q2*(8.402-
		  q2*6.95)));
    
    /// Numerator
    const double num=
      2*M_PI*M_PI*pow(q2,1.5);
    
    /// Denominator
    const double den=
      1+s;
    
    return
      num/den;
  }
  
  /// Full evalation
  double full(const double& q2) const
  {
    /// Square root of q2
    const double q=
      sqrt(q2);
    
    /// Pi^(3/2)
    static const double pi32=
      pow(M_PI,1.5);
    
    /// Numerator
    const double num=
      pi32*q;
    
    /// Denominator
    const double den=
      luscherZetaCalculator(q2);
    
    return
      -num/den;
  }
  
  /// Calculate tan(phi)
  double operator()(const double& q) const
  {
    /// Square of q
    const double q2=
      q*q;
    
    return
      (q2<=0.01)?
      approximation(q2):
      full(q2);
  }
};

/// Structure to compute phi
struct PhiCalculator
{
  /// Calculator of tan
  TanPhiCalculator tanPhiCalculator;
  
  /// Superior limit
  const double q2Max;
  
  /// Numer of known ones (without shifting)
  const vector<double> knownPhiOnes;
  
  /// Raw value of phi without imposing continuity
  double rawPhi(const double& q) const
  {
    return
      atan(tanPhiCalculator(q));
  }
  
  /// Finds all points where the function is approximately one in modulo
  vector<double> findPhiOnes() const
  {
    /// List of ones positions
    vector<double> res;
    
    /// Previous value
    double prev=
      0.0;
    
    /// Step to be used to find ones
    const double step=
      0.1;
    
    /// Starting point
    const double start=
      step;
    
    /// Ending point of the search
    const double end=
      q2Max+2.0;
    
    for(double q2=start;q2<end;q2+=step)
      {
	/// Value of q
	const double q=
	  sqrt(q2);
	
	///Current value
	const double cur=
	  rawPhi(q);
	
	if((fabs(prev)-1.0)*(fabs(cur)-1.0)<0)
	  {
	    /// Position of the one
	    const double pos=
	      Brent_solve([this](const double& q)
	      {
		return
		  fabs(this->rawPhi(q))-1.0;
	      },
		q,sqrt(q2-step));
	    
	    res.push_back(pos);
	  }
	
	prev=cur;
      }
    
    return
      res;
  }
  
  /// Computes the shift needed
  int neededShift(const double& q,
		  const double& rawPhi) const
  {
    /// Find the number of ones before
    int nSmallerOnes=0;
    while(nSmallerOnes<(int)knownPhiOnes.size()-1 and knownPhiOnes[nSmallerOnes+1]<q)
      nSmallerOnes++;
    
    /// Computes the number of pairs of ones
    const int numberOfSmallerPairs=
      (nSmallerOnes+1)/2;
    
    /// Check if it is an even one, in which case the megative part must be shifted
    const int isEvenOne=
      nSmallerOnes%2==0;
    
    /// Check if the raw value of phi is negative
    const int rawPhiIsNegative=
      rawPhi<0;
    
    /// Computes the shift to be added
    const int shift=
      numberOfSmallerPairs+(isEvenOne and rawPhiIsNegative);
    
    return
      shift;
  }
  
  /// Computes putting the correct shift to provide continuity
  double operator()(const double& q) const
  {
    if(q*q>=q2Max)
      CRASH("Cannot compute beyond xMax");
    
    /// Compute rawPhi
    const double rawPhi=
      atan(tanPhiCalculator(q));
    
    /// Shift to be added
    const int shift=
      neededShift(q,rawPhi);
    
    /// Shifted value of phi
    const double shiftedPhi=
      rawPhi+M_PI*shift;
    
    return
      shiftedPhi;
  }
  
  /// Constructor
  PhiCalculator(const double& q2Max) :
    q2Max(q2Max),
    knownPhiOnes(findPhiOnes())
  {
  }
};

/// Provides an hashed version of the tanPhi function
struct HashedTanPhiAndDerivFunction
{
  const PhiCalculator pc;
  
  /// Interpolates for phi(q)/(pi*q^2)
  HashedFunction<> hashedPhiOverPiQ2Function;
  
  /// Interpolates for phi(q)/(pi*q^2)
  HashedFunction<> hashedPhiDeriv;
  
  /// Constructor
  HashedTanPhiAndDerivFunction(const string& pathPhi,
			       const string& pathPhiDeriv,
			       const double& qMax,
			       const int& nIntervals) :
    pc(qMax*qMax),
    hashedPhiOverPiQ2Function(pathPhi,
			      [this](const double& q)
			      {
				return
				  pc(q)/(M_PI*q*q);
			      },
			      0.1,
			      qMax,
			      nIntervals),
      hashedPhiDeriv(pathPhiDeriv,
		     [this](const double& q)
		     {
		       return
			 gslDeriv(pc,q);
		     },
		     0.1,
		     qMax,
		     nIntervals)
  {
  }
  
  /// Evaluate tan(phi(q)) for the required value of q
  double tanPhi(const double& q) const
  {
    return
      tan(hashedPhiOverPiQ2Function(q)*M_PI*q*q);
  }
  
  /// Evaluate phi'(q) for the required value of q
  double phiDeriv(const double& q) const
  {
    return
      hashedPhiDeriv(q);
  }
  
  void plot(const char* phiPath,
	    const char* phiDerivPath) const
  {
    for(auto& p : {make_pair(phiPath,hashedPhiOverPiQ2Function),
		   {phiDerivPath,hashedPhiDeriv}})
      {
	grace_file_t plot(p.first);
	
	plot.set_no_symbol();
	for(double q=p.second.xMin;q<p.second.xMax;q+=0.001)
	  plot.write_xy(q,p.second(q));
      }
  }
};

#endif
