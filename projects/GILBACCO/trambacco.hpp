#include <tranalisi.hpp>

#include <PrecFloat.hpp>
#include <hpIntegrator.hpp>

const map<char,double> aMlist{{'Z',0.00077},{'B',0.00072},{'C',0.0006},{'D',0.00054},{'E',0.00044}};

/// Physical pion mass
constexpr double MPiPhys=0.1350;

/// Physical pion decay constant
constexpr double fPiPhys=0.1304;

/// Estimates the lattice spacing
inline djack_t estimateA(const double& am,
			 const djvec_t& cP5P5,
			 const size_t& T,
			 const string& plotPath)
{
  const auto [Z2Pi,aMPi]=
    two_pts_fit(cP5P5,T/2,30,40,plotPath);
  
  const djack_t aFPi=sqrt(Z2Pi)*2*am/(aMPi*sinh(aMPi));
  
  cout<<"am: "<<am<<endl;
  cout<<"aMPi: "<<aMPi.ave_err()<<endl;
  cout<<"aFPi: "<<aFPi.ave_err()<<endl;
  
  const djack_t xi=aMPi/aFPi;
  cout<<"Xi: "<<xi.ave_err()<<" phys: "<<MPiPhys/fPiPhys<<endl;
  
  const djack_t aInv=
    fPiPhys/aFPi;
  
  return aInv;
}

/////////////////////////////////////////////////////////////////

inline const double& get(const double& x)
{
  return x;
}

inline double get(const PrecFloat& x)
{
  return x.get();
}

template <typename Real,
	  typename F>
Real integrateUpToInfinity(F&& f,
			   const double& lower=0,
			   const double& tol=1e-8)
{
  if constexpr(std::is_same_v<Real,double>)
    return gslIntegrateUpToInfinity(f,lower,tol);
  else
    return precIntegrateUpToInfinity(f,lower);
}

/////////////////////////////////////////////////////////////////

template <typename Real>
struct Basis
{
  /// Time extension of the correlator
  const size_t T;
  
  /// First time from which to reconstruct
  const size_t tMin;
  
  /// Number of basis elements
  const size_t nT;
  
  /// Sign of the backward signal
  const double periodicitySign;
  
  /// Basis function
  ///
  /// The value of "it" is the element of the function, so we need to
  /// increase by tmin to get the time
  Real operator()(const size_t& it,
		  const Real& E) const
  {
    /// Refered time
    const Real t=
      it+tMin;
    
    return exp(-t*E)+periodicitySign*exp(-(T-t)*E);
  }
  
  Basis(const size_t& T,
	const size_t& tMin,
	const size_t& tMax,
	const double& periodicitySign) :
    T(T),tMin(tMin),nT(tMax-tMin),periodicitySign(periodicitySign)
  {
  }
};

/// Functional of the trambacco
template <typename Real>
struct TrambaccoFunctional
{
  VectorX<Real> L;
  
  MatrixX<Real> QStat;
  
  MatrixX<Real> QSyst;
  
  TrambaccoFunctional(const size_t& nT) :
    L(nT),QStat(nT,nT),QSyst(nT,nT)
  {
  }
  
  VectorX<Real> getG(const vector<Real>& preco,
		     const double& statOverSyst=1) const
  {
    const size_t& nT=
      preco.size();
    
    const MatrixX<Real> Q=
      QSyst+QStat/statOverSyst;
    
    const VectorX<Real> gPreco=
      Q.inverse()*L;
    
    VectorX<Real> g(nT);
    
    for(size_t iT=0;iT<preco.size();iT++)
      g[iT]=gPreco[iT]/preco[iT];
    
    return g;
  }
};

template <typename Real,
	  typename SpecAns,
	  typename Targ>
TrambaccoFunctional<Real> getTrambaccoFunctional(const Basis<Real>& basis,
						 const SpecAns& specAns,
						 const Targ& targ,
						 const Real& EMin,
						 const Real& EMax,
						 const MatrixX<Real>& corrCov,
						 const vector<Real>& preco)
{
  static constexpr double tol=1e-8;
  
  using Vect=Vector<Real,Dynamic>;
  
  using Matr=Matrix<Real,Dynamic,Dynamic>;
  
  const size_t& nT=basis.nT;
  
  Vect e(nT);
  
  Matr f(nT,nT);
  
  for(size_t iT=0;iT<nT;iT++)
    {
      e(iT)=
	gslIntegrateFromTo([&targ,&specAns,&basis,iT](const Real& E)
	{
	  return get(targ(E)*sqr(specAns(E))*basis(iT,E));
	},get(EMin),get(EMax),tol);
      
      for(size_t iS=iT;iS<nT;iS++)
	f(iT,iS)=f(iS,iT)=
	  gslIntegrateFromTo([iT,iS,&basis,&specAns](const Real& E)
	  {
	    return get(sqr(specAns(E))*basis(iT,E)*basis(iS,E));
	  },get(EMin),get(EMax),tol);
    }
  
  TrambaccoFunctional<Real> tf(nT);
  
  /// Prepares the functional components
  for(size_t iT=0;iT<nT;iT++)
    {
      tf.L(iT)=e(iT)/preco[iT];
      for(size_t iS=0;iS<nT;iS++)
	{
	  constexpr double corrFading=1;
	  const size_t dist=abs((int)iS-(int)iT);
	  tf.QSyst(iT,iS)=f(iT,iS)/(preco[iT]*preco[iS]);
	  tf.QStat(iT,iS)=corrCov(iT,iS)*pow(corrFading,dist)/(preco[iT]*preco[iS]);
	}
    }
  
  return tf;
}

