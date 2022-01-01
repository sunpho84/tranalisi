#ifndef _TANTALOBACCO_HPP
#define _TANTALOBACCO_HPP

#include <grace.hpp>
#include <highPrec.hpp>
#include <math.hpp>
#include <meas_vec.hpp>

// https://arxiv.org/pdf/1903.06476.pdf

struct TantaloBaccoPars
{
  const bool useTantalo;
  
  const int T;
  
  const int tMin;
  
  const int nT;
  
  const PrecFloat E0;
  
  const PrecFloat lambda;
  
  const PrecFloat sigma;
  
  const PrecFloat alpha=0;
  
  const int useBw;
  
  TantaloBaccoPars(const bool& useTantalo,
		   const int& T,
		   const int& tMin,
		   const int& nT,
		   const PrecFloat& E0,
		   const PrecFloat& lambda,
		   const PrecFloat& sigma,
		   const int& useBw) :
    useTantalo(useTantalo),
    T(T),
    tMin(tMin),
    nT(nT),
    E0(E0),
    lambda(lambda),
    sigma(sigma),
    useBw(useBw)
  {
  }
  
  PrecFloat bT(const PrecFloat& t,
	       const PrecFloat& E) const
    {
      return
	exp(-t*E)+useBw*
	exp(-((PrecFloat)T-t)*E);
    };
  
  PrecFloat Z(const PrecFloat& Estar) const
  {
    return
      (PrecFloat(1)+erf(Estar/(M_SQRT2*sigma)))/2;
  };
  
  PrecFloat N(const PrecFloat& Estar,
	      const PrecFloat& k) const
    {
      return
	(1-lambda)/(2*Z(Estar))*exp((alpha-k)*((alpha-k)*sigma*sigma+2*Estar)/2);
    }
  
  PrecFloat F(const PrecFloat& Estar,
	      const PrecFloat& k) const
    {
      return
	1+erf(((alpha-k)*sigma*sigma+Estar-E0)/(M_SQRT2*sigma));
    };
  
  PrecFloat Delta(const PrecFloat& Estar,
		  const PrecFloat& E) const
  {
    return
      exp(-sqr(E-Estar)/(2*sigma*sigma))/(M_SQRT2*sqrt(M_PI)*sigma*Z(Estar));
  };
};

struct TantaloBaccoRecoEngine :
  TantaloBaccoPars
{
  PrecVect R;
  
  PrecMatr A;
  
  void fillR()
  {
    R.resize(nT,1);
    
    for(int iT=0;iT<nT;iT++)
      R[iT]=
	1.0/PrecFloat(iT+tMin)+useBw*
	1.0/PrecFloat(T-iT-tMin);
    
    grace_file_t RFile("/tmp/R.xmg");
    RFile.new_data_set();
    for(int iT=0;iT<nT;iT++)
      RFile.write_xy(iT,R[iT].get());
  }
  
  static PrecFloat tantalArt(const PrecFloat& Estar,
			     const PrecFloat& E0,
			     const PrecFloat& i)
  {
    return
      exp(-i*E0)/i;
  }
  
  static PrecFloat baccArt(const PrecFloat& Estar,
			   const PrecFloat& E0,
			   const PrecFloat& i)
  {
    return
      ((1/i-Estar)*2/i+Estar*Estar)/i;
  }
  
  void fillA(const PrecFloat& Estar)
  {
    auto art=[this,&Estar](const PrecFloat& i)
    {
      return
	(useTantalo?tantalArt:baccArt)(Estar,E0,i);
    };
    
    A.resize(nT,nT);
    
    for(int iR=0;iR<nT;iR++)
      for(int iT=0;iT<nT;iT++)
	{
	  const PrecFloat a=(PrecFloat)iR+iT+2*tMin;
	  const PrecFloat b=(PrecFloat)T-iR+iT;
	  const PrecFloat c=(PrecFloat)T+iR-iT;
	  const PrecFloat d=(PrecFloat)2*T-iR-iT-2*tMin;
	
	  A(iR,iT)=
	    art(a)+
	    useBw*(art(b)+
		   art(c)+
		   art(d));
	}
    
    grace_file_t AFile("/tmp/A.xmg");
    for(int iR=0;iR<nT;iR++)
      for(int iT=0;iT<nT;iT++)
	AFile.write_xy(iT+nT*iR,A(iR,iT).get());
  }
  
  TantaloBaccoRecoEngine(const TantaloBaccoPars& pars,
			 const PrecFloat& Estar) :
    TantaloBaccoPars(pars)
  {
    fillA(Estar);
    fillR();
  }
};

struct TantaloBaccoReco :
  TantaloBaccoRecoEngine
{
  const PrecFloat Estar;
  
  const djvec_t corr;
  
  PrecVect f;
  
  PrecMatr B;
  
  PrecMatr W;
  
  PrecVect g;
  
  void fillF()
  {
    f.resize(nT,1);
    
    for(int iT=0;iT<nT;iT++)
      f(iT)=
	N(Estar,iT+tMin)*F(Estar,iT+tMin)+useBw*
	N(Estar,T-iT-tMin)*F(Estar,T-iT-tMin);
    
    grace_file_t fFile("/tmp/f"+to_string(Estar.get())+".xmg");
    for(int i=0;i<nT;i++)
      fFile.write_xy(i,f(i).get());
  }
  
  void fillB()
  {
    B.resize(nT,nT);
    
    for(int iR=0;iR<nT;iR++)
      {
	for(int iT=0;iT<nT;iT++)
	  if(lambda!=0)
	    if(useTantalo)
	      B(iR,iT)=sqr(corr[iR+tMin].err()/corr[1].ave())*(iR==iT);
	    else
	      B(iR,iT)=sqr(corr[iR+tMin].err()/corr[1].ave())*(iR==iT);
	  else
	    B(iR,iT)=0.0;
      }
    
    grace_file_t BFile("/tmp/B"+to_string(Estar.get())+".xmg");
    for(int iR=0;iR<nT;iR++)
      for(int iT=0;iT<nT;iT++)
	BFile.write_xy(iT+nT*iR,B(iR,iT).get());
  }
  
  void fillW()
  {
    W.resize(nT,nT);
    
    for(int iR=0;iR<nT;iR++)
      for(int iT=0;iT<nT;iT++)
	W(iR,iT)=
	  (1-lambda)*A(iR,iT)+
	  lambda*B(iR,iT);
    
    grace_file_t WFile("/tmp/W"+to_string(Estar.get())+".xmg");
    for(int iR=0;iR<nT;iR++)
      for(int iT=0;iT<nT;iT++)
	WFile.write_xy(iT+nT*iR,W(iR,iT).get());
  }
  
  void computeG()
  {
    PrecMatr Winv=
      W.inverse();
    
    grace_file_t WFile("/tmp/Winv"+to_string(Estar.get())+".xmg");
    for(int iR=0;iR<nT;iR++)
      for(int iT=0;iT<nT;iT++)
	WFile.write_xy(iT+nT*iR,Winv(iR,iT).get());
    
    const PrecFloat den=
      R.transpose()*Winv*R;
    
    if(useTantalo)
      {
	const PrecFloat num=
	  1-R.transpose()*Winv*f;
	
	g=
	  Winv*f+
	  Winv*R*num/den;
      }
    else
      g=Winv*R/den;
    
    grace_file_t gFile("/tmp/g"+to_string(Estar.get())+".xmg");
    for(int iT=0;iT<nT;iT++)
      gFile.write_xy(iT,g[iT].get());
  }
  
  PrecFloat recoDelta(const PrecFloat& E) const
  {
    PrecFloat s=0;
    
    for(int iT=0;iT<nT;iT++)
      s+=g(iT)*bT(iT+tMin,E);
    
    return s;
  }
  
  djack_t recoDensity() const
  {
    djack_t s;
    
    for(size_t ijack=0;ijack<=njacks;ijack++)
      {
	PrecFloat temp=0.0;
	for(int iT=0;iT<nT;iT++)
	  temp+=g[iT]*corr[iT+tMin][ijack];
	s[ijack]=temp.get();
      }
    
    return s;
  }
  
  PrecFloat recoErrDensity() const
  {
    PrecFloat s=0.0;
    for(int iT=0;iT<nT;iT++)
      s+=sqr(g[iT]*corr[iT+tMin].err());
    
    return sqrt(s);
  }
  
  TantaloBaccoReco(const TantaloBaccoRecoEngine& engine,
		   const PrecFloat& Estar,
		   const djvec_t& corr) :
    TantaloBaccoRecoEngine(engine),
    Estar(Estar),
    corr(corr)
  {
    fillF();
    
    fillB();
    
    fillW();
    
    computeG();
  }
};

#endif
