#ifndef _TANTALOBACCO_HPP
#define _TANTALOBACCO_HPP

#include <grace.hpp>
#include <highPrec.hpp>
#include <invert.hpp>
#include <math.hpp>
#include <meas_vec.hpp>
#include <solve.hpp>

// https://arxiv.org/pdf/1903.06476.pdf

namespace Bacco
{
  /// Holds the parameters related to the correlation function to be analyzed
  struct CorrelatorPars
  {
    /// Time extents
    const int T;
    
    /// Take note if the correlation function has a backward signal
    const bool hasBwSignal;
    
    /// Minimal energy declared for the correlator
    const double E0;
    
    /// Correlation function to be analyzed
    const djvec_t corr;
    
    template <typename F>
    void loopOverIrIt(F&& f,const int tMin,const int nT) const
    {
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  {
	    f(iR,iT,iR+iT+2*tMin);
	    
	    if(hasBwSignal)
	      {
		f(iR,iT,T-iR+iT);
		f(iR,iT,T+iR-iT);
		f(iR,iT,2*T-iR-iT-2*tMin);
	      }
	  }
    }
    
    template <typename F>
    void loopOverIt(F&& f,const int tMin,const int nT) const
    {
	for(int iT=0;iT<nT;iT++)
	  {
	    f(iT,iT+tMin);
	    
	    if(hasBwSignal)
	      f(iT,T-iT-tMin);
	  }
    }
    
    /// Constructor
    CorrelatorPars(const int& T,
		   const bool& hasBwSignal,
		   const double& E0,
		   const djvec_t& corr) :
      T(T),
      hasBwSignal(hasBwSignal),
      E0(E0),
      corr(corr)
    {
    }
    
    /// Basis to be used
    PrecFloat bT(const PrecFloat& t,
		 const PrecFloat& E) const
    {
      PrecFloat res=
	exp(-t*E);
      
      if(hasBwSignal)
	res+=exp(-(T-t)*E);
      
      return res;
    };
    
#define PROVIDE_MOMENT(NAME,SPEC...)		\
    						\
    PrecFloat NAME(const PrecFloat& t) const	\
    {						\
      auto c=					\
	[this](const PrecFloat& t)		\
	{					\
	  return exp(-t*E0)*(SPEC);		\
	};					\
						\
      PrecFloat res=				\
	c(t);					\
						\
      if(hasBwSignal)				\
	res+=c(T-t);				\
						\
      return res;				\
    }
    
    PROVIDE_MOMENT(norm,1/t)
    PROVIDE_MOMENT(firstMoment,(1+E0*t)/(t*t))
    PROVIDE_MOMENT(secondMoment,(2+E0*t*(2+E0*t))/(t*t*t))
    
#undef PROVIDE_MOMENT
    
    /// Integral of E*bT(E) from E0 to infinity
    PrecFloat basisFirstMoment(const PrecFloat& t) const
    {
      auto c=
	[this](const PrecFloat& t)
	{
	  return exp(-t*E0)*(1/t+E0)/t;
	};
      
      PrecFloat res=
	c(t);
      
      if(hasBwSignal)
	res+=c(T-t);
      
      return res;
    };
  };
  
  /// Implements the plottability function
  struct Plottable
  {
    /// Plots the smearing function
    void plot(const string& path,
	      const double& Emin,
	      const double& Emax,
	      const int& nPoints=101) const
    {
      grace_file_t plot(path);
      
      plot.write_line([this](const double& E)
      {
	return this->plotFunction(E);
      },Emin,Emax,nPoints);
      
      plot.set_xaxis_label("E");
      plot.set_yaxis_label("\xD");
      plot.set_no_symbol();
    }
    
    virtual double plotFunction(const double&E) const =0;
  };
  
  /// Holds a reconstrucion of the smearing function
  struct Reconstruction :
    public Plottable
  {
    /// Basis size
    const int nT;
    
    /// Minimal value of t
    const int tMin;
    
    /// Holds the correlation function
    const CorrelatorPars correlatorPars;
    
    /// Parameters of the reconstrucion
    const PrecVect g;
    
    /// Norm of the square within the support
    PrecFloat squareNorm() const
    {
      PrecFloat a=0.0;
      
      correlatorPars.loopOverIrIt([&a,this](const int& iR,
					    const int& iT,
					    const PrecFloat& i)
      {
	a+=g[iT]*g[iR]*exp(-correlatorPars.E0*i)/i;
      },tMin,nT);
      
      return a;
    }
    
    double widthAssumingGaussianAround(const double& E) const
    {
      PrecFloat f=0.0;
      
      correlatorPars.loopOverIrIt([&f,this,&E](const int& iR,
				     const int& iT,
				     const PrecFloat& i)
      {
	f+=
	  g[iR]*g[iT]*((1/i-E)*2/i+E*E)/i;
      },tMin,nT);
      
      const double functional=
	f.get();
      
      const double est=
	functional*4*sqrt(M_PI);
      
      const double res=
	// est;
	NewtonSolve([functional,E](const double& s)
	{
	  const double x=
	    E/s;
	  
	  const double num=
	    2*erf(x)+2-4*x*exp(-x*x)/sqrt(M_PI);
	  
	  const double den=
	    sqr(erf(x/sqrt(2))+1);
	  
	  const double xi=
	    num/den;
	  
	  const double guess=
	    s/(4*sqrt(M_PI));
	  
	  const double fOfS=
	    guess*xi;
	  
	  return fOfS-
	    functional;
	},est);
      
      return res;
    }
    
    /// Computes the reconstrucion function at a given energy
    PrecFloat smearingFunction(const PrecFloat& E) const
    {
      PrecFloat s=0;
      
      for(int iT=0;iT<nT;iT++)
	s+=g[iT]*correlatorPars.bT(iT+tMin,E);
      
      return s;
    }
    
    double plotFunction(const double&E) const
    {
      return smearingFunction(E).get();
    }
    
    /// Integer power, to be used for mean and variance
    static int pow(const int a,const int& n)
    {
      int res=1;
      for(int i=1;i<=n;i++)
	res*=a;
      
      return res;
    }
    
// #define PROVIDE_MOMENT(NAME,SPEC...)					
//     PrecFloat NAME() const						
//     {									
//       PrecFloat m=0.0;							
// 									
//       correlatorPars.loopOverIt([this,&m](const int& iT,		
// 					  const PrecFloat& i)		
//       {									
// 	m+=g[iT]*(SPEC);						
//       },tMin,nT);							
//       									
//       return m;								
//     }
    
//     PROVIDE_MOMENT(norm,correlatorPars.norm(i))
//     PROVIDE_MOMENT(mean,correlatorPars.firstMoment(i))
//     PROVIDE_MOMENT(meanOfE2,correlatorPars.secondMoment(i))
    
    /// Computes the estimator of the width as <(E-<E>)^2>
    ///
    /// Notice that since the smearing function is not garanteed to be
    /// definite positive, the width can be undefinite
    // PrecFloat width() const
    // {
    //   return sqrt(meanOfE2()-sqr(mean()));
    // }
    
    /// Reconstruct the density for the fixed energy
    djack_t recoDensity() const
    {
      djack_t s;
      
      for(size_t ijack=0;ijack<=njacks;ijack++)
	{
	  PrecFloat temp=0.0;
	  for(int iT=0;iT<nT;iT++)
	    temp+=g[iT]*correlatorPars.corr[iT+tMin][ijack];
	  s[ijack]=temp.get();
	}
      
      return s;
    }
    
    /// Default constructor
    Reconstruction(const int& nT,
		   const int& tMin,
		   const CorrelatorPars& correlatorPars,
		   const PrecVect& g) :
      nT(nT),
      tMin(tMin),
      correlatorPars(correlatorPars),
      g(g)
    {
    }
  };
  
  struct Reconstructor
  {
    const CorrelatorPars correlatorPars;
    
    const djvec_t& corr=
      correlatorPars.corr;
    
    const int& T=
      correlatorPars.T;
    
    const bool& hasBwSignal=
      correlatorPars.hasBwSignal;
    
    const int tMin;
    
    const int nT;
    
    const double Estar;
    
    const double lambda;
    
    /// Basis normalization
    PrecVect R;
    
    virtual PrecVect prepareWVect(const PrecMatr& Winv) const =0;
    
    virtual PrecFloat aFun(const PrecFloat& i) const =0;
    
    Reconstruction getReco()
    {
      // Fill R
      
      R.resize(nT,1);
      
      for(int iT=0;iT<nT;iT++)
	{
	  R[iT]=
	    1.0/PrecFloat(iT+tMin);
	  
	  if(hasBwSignal)
	    R[iT]+=
	      1.0/PrecFloat(T-iT-tMin);
	}
      
      grace_file_t RFile("/tmp/R.xmg");
      RFile.new_data_set();
      for(int iT=0;iT<nT;iT++)
	RFile.write_xy(iT,R[iT].get());
      
      /////////////////////////////////////////////////////////////////
      
      // fill A
      
      /// Orthonormalization
      PrecMatr A(nT,nT);
      A.setZero();
      
      correlatorPars.loopOverIrIt([&A,this](const int& iR,
					    const int& iT,
					    const int& i)
      {
	A(iR,iT)+=
	  aFun(i);
      },tMin,nT);
      
      grace_file_t AFile("/tmp/A.xmg");
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  AFile.write_xy(iT+nT*iR,A(iR,iT).get());
      
      /////////////////////////////////////////////////////////////////
      
      // fill B
      
      PrecMatr B(nT,nT);
      
      for(int iR=0;iR<nT;iR++)
	{
	  for(int iT=0;iT<nT;iT++)
	    if(lambda!=0)
	      B(iR,iT)=sqr(corr[iR+tMin].err()/corr[1].ave())*(iR==iT);
	    else
	      B(iR,iT)=0.0;
	}
      
      grace_file_t BFile("/tmp/B"+to_string(Estar)+".xmg");
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  BFile.write_xy(iT+nT*iR,B(iR,iT).get());
      
      /////////////////////////////////////////////////////////////////
      
      // fill W
      
      PrecMatr W(nT,nT);
      
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  W(iR,iT)=
	    (1-lambda)*A(iR,iT)+
	    lambda*B(iR,iT);
      
      grace_file_t WFile("/tmp/W"+to_string(Estar)+".xmg");
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  WFile.write_xy(iT+nT*iR,W(iR,iT).get());
      
      /////////////////////////////////////////////////////////////////
      
      PrecMatr Winv=W.inverse();
      
      // cout<<"Check inversion: "<<(Winv*W-PrecMatr::Identity(nT,nT)).squaredNorm().get()<<" "<<(W*Winv-PrecMatr::Identity(nT,nT)).squaredNorm().get()<<endl;
      
      grace_file_t WinvFile("/tmp/Winv"+to_string(Estar)+".xmg");
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  WinvFile.write_xy(iT+nT*iR,Winv(iR,iT).get());
      
      const PrecVect v=
	prepareWVect(Winv);
      
      /////////////////////////////////////////////////////////////////
      
      PrecVect g(nT);
      
      for(int iR=0;iR<nT;iR++)
	{
	  g[iR]=0.0;
	  for(int iC=0;iC<nT;iC++)
	    g[iR]+=
	      Winv(iR,iC)*v[iC];
	}
      
      grace_file_t gFile("/tmp/g"+to_string(Estar)+".xmg");
      for(int iT=0;iT<nT;iT++)
	gFile.write_xy(iT,g[iT].get());
      
      Reconstruction reco(nT,tMin,correlatorPars,g);
      
      return reco;
    }
    
    Reconstructor(const CorrelatorPars& correlatorPars,
			 const int& tMin,
			 const int& tMax,
			 const double& Estar,
			 const double& lambda) :
    correlatorPars(correlatorPars),
    tMin(tMin),
    nT(tMax-tMin),
    Estar(Estar),
    lambda(lambda)
    {
    }
  };
  
  struct TargetedReconstructor :
    public Reconstructor,
    public Plottable
  {
    /// Normalize the weight with E^m
    const int mFact;
    
    const double E0;
    
    const PrecFloat alpha=0;
    
    /// In the targeted reconstruction, a is set to the orthonormality of the basis vector
    ///
    /// See eq.32 of Nazario's paper
    PrecFloat aFun(const PrecFloat& i) const
    {
      if(mFact!=0 and mFact!=2)
	CRASH("Not implemented");
      
      if(mFact==0)
	return
	  exp(-i*E0)/i;
      else
	return
	  (24+E0*i*(24+E0*i*(12+E0*i*(4+E0*i))))/(exp(E0*i)*pow(i,5));
      //(gamma(1+m)/pow(i,m)+(pow(E0,m)*(-(m*gamma(m))+gamma(1+m,E0*i)))/pow(E0*i,m))/i; wrong, check, maybe is just the ordinary gamma
    }
    
    // virtual PrecFloat normalization() const =0;
    
    virtual PrecFloat fFun(const PrecFloat& i) const =0;
    
    virtual double targetFunction(const double& E) const =0;
    
    double plotFunction(const double&E) const
    {
      return targetFunction(E);
    }
    
    /// Returns the actual squared norm, in the support
    virtual PrecFloat squareNorm() const =0;
    
    PrecVect prepareWVect(const PrecMatr& Winv) const
    {
      PrecVect f(nT);
      
      for(int iT=0;iT<nT;iT++)
	{
	  f[iT]=fFun(iT+tMin);
	  
	  if(hasBwSignal)
	    f[iT]+=
	      fFun(T-iT-tMin);
	}
      
      grace_file_t fFile("/tmp/f"+to_string(Estar)+".xmg");
      for(int i=0;i<nT;i++)
	fFile.write_xy(i,f[i].get());
      
      /////////////////////////////////////////////////////////////////
      
      PrecVect res(nT);
      
      const PrecFloat c=
	(1-(PrecFloat)lambda);
      
      // const PrecFloat num=
      // 	normalization()-c*R.transpose()*Winv*f;
      
      // const PrecFloat den=
      // 	R.transpose()*Winv*R;
      
      for(int iC=0;iC<nT;iC++)
	res[iC]=
	  c*f[iC]// +R[iC]*num/den
	  ;
      
      return res;
    }
    
    double recoErr(const double& E,
		   const Reconstruction& reco) const
    {
      return targetFunction(E)-reco.smearingFunction(E).get();
    }
    
    PrecFloat projectionWithReco(const Reconstruction& reco) const
    {
      PrecFloat res=0.0;
      
      for(int iT=0;iT<nT;iT++)
	{
	  PrecFloat f=
	    fFun(iT+tMin);
	  
	  if(hasBwSignal)
	    f+=
	      fFun(T-iT-tMin);
	  
	  res+=
	    reco.g[iT]*f;
	}
      
      return res;
    }
    
    double deviation(const Reconstruction& reco) const
    {
      return (reco.squareNorm()-2*projectionWithReco(reco)+squareNorm()).get();
    }
    
    TargetedReconstructor(const int mFact,
			  const CorrelatorPars& correlatorPars,
			  const int& tMin,
			  const int& tMax,
			  const double& Estar,
			  const double& lambda,
			  const double& E0) :
      Reconstructor(correlatorPars,tMin,tMax,Estar,lambda),
      mFact(mFact),
      E0(E0)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  struct GaussReconstructor :
    public TargetedReconstructor
  {
    const double sigma;
    
    PrecFloat fFun(const PrecFloat& t) const
    {
      return exp(-t*Estar+sqr(t*sigma)/2)*erfc((t*sigma*sigma+E0-Estar)/(sigma*sqrt((PrecFloat)2)))/2;
    }
    
    PrecFloat squareNorm() const
    {
      const PrecFloat x=
	((PrecFloat)E0-Estar)/sigma;
      
      return erfc(x)/(4*sqrt(precPi())*sigma);
    }
    
    double targetFunction(const double& E) const
    {
      return
	exp(-sqr((E-Estar)/sigma)/2)/(sqrt(2*M_PI)*sigma);
    };
    
    GaussReconstructor(const CorrelatorPars& correlatorPars,
		       const int& tMin,
		       const int& tMax,
		       const double& Estar,
		       const double& lambda,
		       const double& E0,
		       const double& sigma) :
      TargetedReconstructor(0,correlatorPars,tMin,tMax,Estar,lambda,E0),
      sigma(sigma)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  struct NumericalReconstructor :
    public TargetedReconstructor
  {
    PrecFloat squareNorm() const
    {
      return integrateUpToInfinite([this](const PrecFloat& E)
      {
	return sqr(preciseTargetFunction(E));
      },E0);
    }
    
    virtual PrecFloat preciseTargetFunction(const PrecFloat& E) const=0;
    
    double targetFunction(const double& E) const
    {
      return preciseTargetFunction(E).get();
    }
    
    double deviation2(const Reconstruction& reco) const
    {
      return integrateUpToInfinite([this,&reco](const PrecFloat& E)
      {
	return sqr(preciseTargetFunction(E)-reco.smearingFunction(E));
      },E0).get()/
	integrateUpToInfinite([this](const PrecFloat& E)
      {
	return sqr(preciseTargetFunction(E));
      },E0).get();
    }
    
    PrecFloat fFun(const PrecFloat& t) const
    {
      return integrateUpToInfinite([t,this](const PrecFloat& E)
      {
	return pow(E,2*mFact)*preciseTargetFunction(E)*exp(-t*E);//not using bT since bw signal taken explicitly into account correlatorPars.bT(t,E);
      },E0);
    }
    
    NumericalReconstructor(const int& mFact,
			   const CorrelatorPars& correlatorPars,
			   const int& tMin,
			   const int& tMax,
			   const double& Estar,
			   const double& lambda,
			   const double& E0) :
      TargetedReconstructor(mFact,correlatorPars,tMin,tMax,Estar,lambda,E0)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  struct GaussDivE2Reconstructor :
    public NumericalReconstructor
  {
    const double sigma;
    
    PrecFloat preciseTargetFunction(const PrecFloat& E) const
    {
      return
	exp(-sqr((E-Estar)/sigma)/2)/(sqrt(2*precPi())*sigma)
	/sqr(E);
    }
    
    GaussDivE2Reconstructor(const CorrelatorPars& correlatorPars,
			    const int& tMin,
			    const int& tMax,
			    const double& Estar,
			    const double& lambda,
			    const double& E0,
			    const double sigma) :
      NumericalReconstructor(0,correlatorPars,tMin,tMax,Estar,lambda,E0),
      sigma(sigma)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  struct LegoReconstructor :
    public TargetedReconstructor
  {
    const double width;
    
    PrecFloat fFun(const PrecFloat& t) const
    {
      return
	(exp(-(Estar-width/2)*t)-exp(-(Estar+width/2)*t))/(width*t);
    }
    
    double targetFunction(const double& E) const
    {
      return
	((E>Estar-width/2) and (E<Estar+width/2))/width;
    };
    
    PrecFloat squareNorm() const
    {
      return PrecFloat(1)/width;
    }
    
    LegoReconstructor(const CorrelatorPars& correlatorPars,
			const int& tMin,
			const int& tMax,
			const double& Estar,
			const double& lambda,
			const double& E0,
			const double& width) :
      TargetedReconstructor(0,correlatorPars,tMin,tMax,Estar,lambda,E0),
      width(width)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  struct BGReconstructor :
    public Reconstructor
  {
    PrecVect prepareWVect(const PrecMatr& Winv) const
    {
      PrecVect res(nT);
      
      const PrecFloat den=
	R.transpose()*Winv*R;
      
      for(int iC=0;iC<nT;iC++)
	res[iC]=
	  R[iC]/den;
      
      return res;
    }
    
    PrecFloat aFun(const PrecFloat& i) const
    {
      return
	((1/i-Estar)*2/i+Estar*Estar)/i;
    }
    
    BGReconstructor(const CorrelatorPars &correlatorPars,const int &tMin,const int &tMax,const double &Estar,const double &lambda) :
      Reconstructor(correlatorPars,tMin,tMax,Estar,lambda)
    {
    }
  };
}

// struct TantaloBaccoRecoEngine :
//   TantaloBaccoPars
// {
  
  
//   TantaloBaccoRecoEngine(const TantaloBaccoPars& pars,
// 			 const PrecFloat& Estar) :
//     TantaloBaccoPars(pars)
//   {
//     fillA(Estar);
//     fillR();
//     fillM();
//   }
// };

// struct TantaloBaccoReco :
//   TantaloBaccoRecoEngine
// {
//   enum TargetFunction{GAUSS,LEGO};
  
//   TargetFunction targetFunction{GAUSS};
  
//   const PrecFloat Estar;
  
//   const djvec_t corr;
  
//   /////////////////////////////////////////////////////////////////
  
//   PrecFloat recoErrDensity() const
//   {
//     PrecFloat s=0.0;
//     for(int iT=0;iT<nT;iT++)
//       s+=sqr(g[iT]*corr[iT+tMin].err());
    
//     return sqrt(s);
//   }
  
//   TantaloBaccoReco(const TantaloBaccoRecoEngine& engine,
// 		   const PrecFloat& Estar,
// 		   const djvec_t& corr) :
//     TantaloBaccoRecoEngine(engine),
//     Estar(Estar),
//     corr(corr)
//   {
//   }
// };

#endif
