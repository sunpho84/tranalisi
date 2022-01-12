#ifndef _TANTALOBACCO_HPP
#define _TANTALOBACCO_HPP

#include <grace.hpp>
#include <highPrec.hpp>
#include <invert.hpp>
#include <math.hpp>
#include <meas_vec.hpp>

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
		 const double& E) const
    {
      PrecFloat res=
	exp((PrecFloat)(-t*E));
      
      if(hasBwSignal)
	res+=exp((PrecFloat)(-(T-t)*E));
      
      return res;
    };
  };
  
  struct Reconstruction
  {
    const int nT;
    
    const int tMin;
    
    const CorrelatorPars correlatorPars;
    
    const PrecVect g;
    
    PrecFloat smearingFunction(const double& E) const
    {
      PrecFloat s=0;
      
      for(int iT=0;iT<nT;iT++)
	s+=g[iT]*correlatorPars.bT(iT+tMin,E);
      
      return s;
    }
    
    /// Integer power, to be used for future reference
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
    
    PrecFloat mean() const
    {
      return gWeightedWithTToMinusN(2);
    }
    
    PrecFloat meanOfSquare() const
    {
      return 2*gWeightedWithTToMinusN(3);
    }
    
    PrecFloat widthOfSquare() const
    {
      PrecFloat w2=0.0;
      const int& T=
	correlatorPars.T;
      
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  {
	    auto add=[&w2,gtr=g[iT]*g[iR]](const int i){w2+=gtr/i;};
	    
	    add(iR+iT+2*tMin);
	    
	    if(correlatorPars.hasBwSignal)
	      for(const int i : {T-iR+iT,
				 T+iR-iT,
				 2*T-iR-iT-2*tMin})
		add(i);
	  }
      
      return sqrt(w2);
    }
    
    PrecFloat width() const
    {
      return sqrt(meanOfSquare()-sqr(mean()));
    }
    
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
  
  struct ReconstructionEngine
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
      
      for(int iR=0;iR<nT;iR++)
	for(int iT=0;iT<nT;iT++)
	  {
	    A(iR,iT)=
	      aFun(iR+iT+2*tMin);
	    
	    if(hasBwSignal)
	      {
		A(iR,iT)+=
		  aFun(T-iR+iT)+
		  aFun(T+iR-iT)+
		  aFun(2*T-iR-iT-2*tMin);
	      }
	  }
      
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
	  WFile.write_xy(iT+nT*iR,Winv(iR,iT).get());
      
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
    
    ReconstructionEngine(const CorrelatorPars& correlatorPars,
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
  
  struct TargettedReconstructor :
    public ReconstructionEngine
  {
    const PrecFloat alpha=0;
    
    const double sigma;
    
    const PrecVect f;
    
    virtual PrecFloat fFun(const PrecFloat& i) const =0;
    
    virtual PrecFloat targetFunction(const double& E) const =0;
    
    PrecVect prepareWVect(const PrecMatr& Winv) const
    {
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
    
    PrecVect fillF()
    {
      PrecVect f(nT,1);
      
      for(int iT=0;iT<nT;iT++)
	{
	  f[iT]=fFun(iT+tMin);
	  
	  if(hasBwSignal)
	    f[iT]+=
	      fFun(T-T-iT-tMin);
	}
      
      grace_file_t fFile("/tmp/f"+to_string(Estar)+".xmg");
      for(int i=0;i<nT;i++)
	fFile.write_xy(i,f[i].get());
      
      return f;
    }
    
    TargettedReconstructor(const CorrelatorPars& correlatorPars,
			   const int& tMin,
			   const int& tMax,
			   const double& Estar,
			   const double& lambda,
			   const double& sigma) :
      ReconstructionEngine(correlatorPars,tMin,tMax,Estar,lambda),
      sigma(sigma),
      f(fillF())
    {
    }
  };
  
  struct GaussReconstructor :
    public TargettedReconstructor
  {
    const double E0;
    
    const PrecFloat Z;
    
    static PrecFloat Zfun(const double& Estar,
			  const double& sigma)
    {
      return
	(PrecFloat(1)+erf(Estar/(sqrt(PrecFloat(2))*sigma)))/2;
    };
    
    PrecFloat aFun(const PrecFloat& i) const
    {
      return
	exp(-i*E0)/i;
    }
    
    PrecFloat N(const PrecFloat& Estar,
		const PrecFloat& k) const
    {
      return
	(1-lambda)/(2*Z)*exp((alpha-k)*((alpha-k)*sigma*sigma+2*Estar)/2);
    }
    
    PrecFloat F(const PrecFloat& Estar,
		const PrecFloat& k) const
    {
      return
	1+erf(((alpha-k)*sigma*sigma+Estar-E0)/(M_SQRT2*sigma));
    }
    
    PrecFloat fFun(const PrecFloat& t) const
    {
      return
	N(Estar,t)*F(Estar,t);
    }
    
    PrecFloat targetFunction(const double& E) const
    {
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
      TargettedReconstructor(correlatorPars,tMin,tMax,Estar,lambda,sigma),
      E0(E0),
      Z(Zfun(Estar,sigma))
    {
    }
  };
  
  struct BGReconstructor :
    public ReconstructionEngine
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
      ReconstructionEngine(correlatorPars,tMin,tMax,Estar,lambda)
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
