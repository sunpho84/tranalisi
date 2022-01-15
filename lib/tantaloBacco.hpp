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
    PrecFloat bT(const int& t,
		 const PrecFloat& E) const
    {
      PrecFloat res=
	exp(-t*E);
      
      if(hasBwSignal)
	res+=exp(-(T-t)*E);
      
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
	NewtonSolve([functional,E](const double& s)
	{
	  return (-2*exp(-sqr(E/s))+sqrt(M_PI)*s*(1+erf(E/s)))/
	    (2*M_PI*sqr(1+erf(E/(sqrt(2)*s))))-
	    functional;
	},est);
      
      return res;
    }
    
    /// Computes the reconstrucion function at a given energy
    PrecFloat smearingFunction(const double& E) const
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
    
    /// Compute \sum_it g[it]/(it+tMin)^n with periodicity
    PrecFloat gWeightedWithTToMinusN(const int n) const
    {
      PrecFloat m=0.0;
      
      for(int iT=0;iT<nT;iT++)
	{
	  m+=g[iT]/pow(iT+tMin,n);
	  if(correlatorPars.hasBwSignal)
	    m+=g[iT]/pow(correlatorPars.T-iT-tMin,n);
	}
      
      return m;
    }
    
    /// Returns the normalization, which should be 1
    PrecFloat norm() const
    {
      return gWeightedWithTToMinusN(1);
    }
    
    /// Returns the mean
    PrecFloat mean() const
    {
      return gWeightedWithTToMinusN(2);
    }
    
    /// Returns the expectation value of E^2
    PrecFloat meanOfE2() const
    {
      return 2*gWeightedWithTToMinusN(3);
    }
    
    /// Computes the estimator of the width as <(E-<E>)^2>
    ///
    /// Notice that since the smearing function is not garanteed to be
    /// definite positive, the width can be undefinite
    PrecFloat width() const
    {
      return sqrt(meanOfE2()-sqr(mean()));
    }
    
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
      PrecMatr A(nT,nT,0.0);
      
      correlatorPars.loopOverIrIt([&A,this](const int& iR,
					    const int& iT,
					    const int& i)
      {
	A(iR,iT)+=
	  aFun(iR+iT+2*tMin);
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
      
      PrecMatr Winv(nT,nT);
      invert_matrix(W,Winv);
      
      //cout<<"Check inversion: "<<(Winv*W-PrecMatr::Identity(nT,nT)).squaredNorm().get()<<" "<<(W*Winv-PrecMatr::Identity(nT,nT)).squaredNorm().get()<<endl;
      
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
    const double E0;
    
    const PrecFloat alpha=0;
    
    /// In the targeted reconstruction, a is set to the orthonormality of the basis vector
    ///
    /// See eq.32 of Nazario's paper
    PrecFloat aFun(const PrecFloat& i) const
    {
      return
	exp(-i*E0)/i;
    }
    
    virtual PrecFloat fFun(const PrecFloat& i) const =0;
    
    virtual double targetFunction(const double& E) const =0;
    
    double plotFunction(const double&E) const
    {
      return targetFunction(E);
    }
    
    virtual PrecFloat squareNorm() const =0;
    
    PrecVect prepareWVect(const PrecMatr& Winv) const
    {
      PrecVect f(nT,1);
      
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
      
      const PrecFloat num=
	(1-Winv.formWith(R,f));
      
      const PrecFloat den=
	Winv.formWith(R,R);
      
      for(int iC=0;iC<nT;iC++)
	res[iC]=
	  f[iC]+R[iC]*num/den;
      
      return res;
    }
    
    double recoErr(const double& E,
		   const Reconstruction& reco) const
    {
      return targetFunction(E)-reco.smearingFunction(E).get();
    }
    
    double deviation(const Reconstruction& reco) const
    {
      PrecFloat a=0.0;
      correlatorPars.loopOverIrIt([&a,reco](const int& iR,
					    const int& iT,
					    const int& i)
      {
	      a+=reco.g[iT]*reco.g[iR]/i;
      },tMin,nT);
      
      
      for(int iT=0;iT<nT;iT++)
	{
	  PrecFloat f=
	    fFun(iT+tMin);
	  
	  if(hasBwSignal)
	    f+=
	      fFun(T-iT-tMin);
	  
	  a-=
	    2*reco.g[iT]*f;
	}
      
      return sqrt(a+squareNorm()).get();
    }
    
    TargetedReconstructor(const CorrelatorPars& correlatorPars,
			   const int& tMin,
			   const int& tMax,
			   const double& Estar,
			   const double& lambda,
			   const double& E0) :
      Reconstructor(correlatorPars,tMin,tMax,Estar,lambda),
      E0(E0)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  struct GaussReconstructor :
    public TargetedReconstructor
  {
    const double sigma;
    
    PrecFloat zFun() const
    {
      return
	(PrecFloat(1)+erf(Estar/(sqrt(PrecFloat(2))*sigma)))/2;
    }
    
    PrecFloat fFun(const PrecFloat& t) const
    {
      const PrecFloat Z=
	zFun();
      
      const PrecFloat N=
	(1-lambda)/(2*Z)*exp((alpha-t)*((alpha-t)*sigma*sigma+2*Estar)/2);
      
      const PrecFloat F=
	1+erf(((alpha-t)*sigma*sigma+Estar-E0)/(sqrt(PrecFloat(2))*sigma));
      
      return N*F;
    }
    
    PrecFloat squareNorm() const
    {
      const PrecFloat x=
	((PrecFloat)E0-Estar)/sigma;
      
      return erfc(x)/(sqrt(precPi())*sigma*sqr(erfc(x/sqrt(PrecFloat(2)))));
    }
    
    double targetFunction(const double& E) const
    {
      const double Z=
	zFun().get();
      
      return
	exp(-sqr(E-Estar)/(2*sigma*sigma))/(M_SQRT2*sqrt(M_PI)*sigma*Z);
    };
    
    GaussReconstructor(const CorrelatorPars& correlatorPars,
			const int& tMin,
			const int& tMax,
			const double& Estar,
			const double& lambda,
			const double& E0,
			const double& sigma) :
      TargetedReconstructor(correlatorPars,tMin,tMax,Estar,lambda,E0),
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
      TargetedReconstructor(correlatorPars,tMin,tMax,Estar,lambda,E0),
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
	Winv.formWith(R,R);
      
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
