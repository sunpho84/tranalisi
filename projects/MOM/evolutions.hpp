#ifndef _EVOLUTIONS_HPP
#define _EVOLUTIONS_HPP

#include <tranalisi.hpp>

namespace evolution
{
  //! tag the bilinears, or Zq
  enum bil_t{Q=0,P=1,S=1,T=2};
  
  //! tag the number of quarks
  enum Nf_t{NF0=0,NF2=1,NF4=2};
  
  //! forward declaration of coefficient functions
  template <int Iz,int Nf,int Ord>
  constexpr double C();
  
  //! instantiate the table of C
#define DEFC(IZ,NF,ORD,VAL)			\
  template<>					\
  constexpr double C<IZ,NF,ORD>()		\
  {return VAL;}
  
  //Q,Nf0
  DEFC(Q,NF0,1, 2.0303);
  DEFC(Q,NF0,2, 42.1268);
  DEFC(Q,NF0,3, 1728.43);
  //Q,Nf2
  DEFC(Q,NF2,1, 2.03448);
  DEFC(Q,NF2,2, 35.9579);
  DEFC(Q,NF2,3, 1199.16);
  //Q,Nf4
  DEFC(Q,NF4,1, 2.4000);
  DEFC(Q,NF4,2, 29.6724);
  DEFC(Q,NF4,3, 719.141);
  //P or S,Nf0
  DEFC(P,NF0,1, -8.08264);
  DEFC(P,NF0,2, -151.012);
  DEFC(P,NF0,3, -5247.93);
  //P or S,Nf2
  DEFC(P,NF2,1, -8.55727);
  DEFC(P,NF2,2, -125.423);
  DEFC(P,NF2,3, -3797.71);
  //P or S,Nf4
  DEFC(P,NF4,1, -9.38987);
  DEFC(P,NF4,2, -96.2883);
  DEFC(P,NF4,3, -2403.82);
  //T,Nf0
  DEFC(T,NF0,1, 2.53260);
  DEFC(T,NF0,2, 57.8740);
  DEFC(T,NF0,3, 0.0);
  //T,Nf2
  DEFC(T,NF2,1, 2.66852);
  DEFC(T,NF2,2, 47.9701);
  DEFC(T,NF2,3, 0.0);
  //T,Nf4
  DEFC(T,NF4,1, 2.91662);
  DEFC(T,NF4,2, 37.9471);
  DEFC(T,NF4,3, 0.0);
  
  //! forward declaration of anomalous dimension numerator
  template<int Iz>
  constexpr double AdimNum();
  
  //! instantiate the table of Adim num
#define DEFA_NUM(IZ,VAL)			\
  template<>					\
  constexpr double AdimNum<IZ>()		\
  {return VAL;}
  
  DEFA_NUM(Q, 0.0);
  DEFA_NUM(P, -12.0);
  DEFA_NUM(T, +4.0);
  
  //! forward declaration of anomalous dimension denominator
  template<int Nf>
  constexpr double AdimDen();
  
  //! instantiate the table of Adim den
#define DEFA_DEN(NF,VAL)			\
  template<>					\
  constexpr double AdimDen<NF>()		\
  {return VAL;}
  
  DEFA_DEN(NF0, 33.0);
  DEFA_DEN(NF2, 29.0);
  DEFA_DEN(NF4, 25.0);
  
  //! anomalous dimensiom
  template<int Iz,int Nf>
  constexpr double Adim()
  {
    return AdimNum<Iz>()/AdimDen<Nf>();
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //! strong alpha at three loops, ainv is the inverse lattice spacing, mu2 is a2p2
  template <int Nf,int Nc=3>
  double alphas(double ainv,double mu2,const double lambdaQCD=0.250)
  {
    constexpr int CF=(Nc*Nc-1.0)/(2.0*Nc);
    constexpr double Z3=1.2020569031595942;
    constexpr double beta_0=(11.0*Nc-2.0*Nf)/3.0;
    constexpr double beta_1=34.0/3.0*Sqr(Nc)-10.0/3.0*Nc*Nf-2.0*CF*Nf;
    constexpr double beta_2=(2857.0/54.0)*Cube(Nc)+Sqr(CF)*Nf-205.0/18.0*CF*Nc*Nf -1415.0/54.0*Sqr(Nc)*Nf+11.0/9.0*CF*Sqr(Nf)+79.0/54.0*Nc*Sqr(Nf);
    constexpr double beta_3=(150653.0/486.0-44.0/9.0*Z3)*Pow<4>(Nc) +
      (-39143.0/162.0+68.0/3.0*Z3)*Cube(Nc)*Nf+
      (7073.0/486.0-328.0/9.0*Z3)*CF*Sqr(Nc)*Nf+
      (-2102.0/27.0+176.0/9.0*Z3)*Sqr(CF)*Nc*Nf+
      23.0*Cube(CF)*Nf+(3965.0/162.0+56.0/9.0*Z3)*Sqr(Nc)*Sqr(Nf)+
      (338.0/27.0-176.0/9.0*Z3)*Sqr(CF)*Sqr(Nf)+
      (4288.0/243.0+112.0/9.0*Z3)*CF*Nc*Sqr(Nf)+53.0/243.0*Nc*Cube(Nf)+
      154.0/243.0*CF*Cube(Nf)+
      (-10.0/27.0+88.0/9.0*Z3)*Sqr(Nc)*(Sqr(Nc)+36.0)+
      (32.0/27.0-104.0/9.0*Z3)*Nc*(Sqr(Nc)+6)*Nf+
      (-22.0/27.0+16.0/9.0*Z3)*(Pow<4>(Nc)-6.0*Sqr(Nc)+18.0)/Sqr(Nc)*Sqr(Nf);
    
   constexpr double b1=beta_1/beta_0/4.0/M_PI;
   constexpr double b2=beta_2/beta_0/16.0/Sqr(M_PI);
   constexpr double b3=beta_3/beta_0/64.0/Cube(M_PI);
   
   const double lam0=lambdaQCD/ainv;
   
   const double L2=log(mu2/(lam0*lam0));
   const double LL2=log(L2);
    
   const double als0=4.0*M_PI/beta_0/L2;
   const double als1=als0-pow(als0,2)*b1*LL2;
   const double als2=als1+pow(als0,3)*(pow(b1,2)*(pow(LL2,2)-LL2-1.0)+b2);
   const double als3=als2+pow(als0,4)*(pow(b1,3)*(-pow(LL2,3)+5.0/2.0*pow(LL2,2)+2*LL2-1.0/2.0)-3.0*b1*b2*LL2+b3/2.0);
   
   return als3;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////
  
  //! compute the evolution up to 1/a2 of a Z
  template<int Iz,int Nf,int Ord>
  double evolution_to_RIp(double ainv,double mu_2)
  {
    double mu[2]={mu_2,1.0}; //convert to scale 1/a2 in units of a2, that is 1
    double cmu[2];
    
    for(size_t i=0;i<2;i++)
      {
	double al=alphas<Nf>(ainv,mu[i])/(4*M_PI);
	cmu[i]=1.0;
	cmu[i]+=al*C<Iz,Nf,1>();
	cmu[i]+=al*al*C<Iz,Nf,2>();
	cmu[i]+=al*al*al*C<Iz,Nf,3>();
	cmu[i]*=pow(al,Adim<Iz,Nf>());
      }
    
    return cmu[0]/cmu[1];
  }
}

#endif
