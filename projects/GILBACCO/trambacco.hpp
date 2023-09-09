#include <tranalisi.hpp>

template <typename Real>
struct Trambacco
{
  /// Time extension of the correlator
  const size_t T;
  
  /// First time from which to reconstruct
  const size_t& tMin;
  
  /// Sign of the backward signal
  const double periodicitySign;
  
  /// Basis function
  ///
  /// The value of "it" is the element of the function, so we need to
  /// increase by tmin to get the time
  Real b(const size_t& it,
	 const Real& E) const
  {
    /// Refered time
    const Real t=
      it+tMin;
    
    return exp(-t*E)+exp(-(T-t)*E);
  }
  
  Trambacco(const size_t& T,
	    const size_t& tMin,
	    const double& periodicitySign=+1) :
    T(T),tMin(tMin),periodicitySign(periodicitySign)
  {
  }
};
