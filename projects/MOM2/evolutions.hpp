#ifndef _EVOLUTIONS_HPP
#define _EVOLUTIONS_HPP

#include <pr_bil.hpp>

namespace evol
{
  extern double lambda_QCD;
  
  //! tag the bilinears, or Zq
  enum ev_bil_t{Q=0,P=1,S=1,T=2};
  
  //! tag the number of quarks
  enum Nf_t{NF0=0,NF2=2,NF4=4};
  
  //! map Nf to Nf_t
  inline Nf_t Nf_t_of_Nf(int Nf)
  {
    switch(Nf)
      {
      case 0: return Nf_t::NF0; break;
      case 2: return Nf_t::NF2; break;
      case 4: return Nf_t::NF4; break;
      default:
	CRASH("Unknwonw Nf %d",Nf);
	return Nf_t::NF0;
      }
  }
  
  //! strong alpha at three loops, mu2 in physical units
  template <int Nf,int Nc=3>
  double alphas(double mu2)
  {
    constexpr double CF=(Nc*Nc-1.0)/(2.0*Nc);
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
    
    const double L2=log(mu2/sqr(lambda_QCD));
    const double LL2=log(L2);
    
    const double als0=4.0*M_PI/beta_0/L2;
    const double als1=als0-pow(als0,2)*b1*LL2;
    const double als2=als1+pow(als0,3)*(pow(b1,2)*(pow(LL2,2)-LL2-1.0)+b2);
    const double als3=als2+pow(als0,4)*(pow(b1,3)*(-pow(LL2,3)+5.0/2.0*pow(LL2,2)+2*LL2-1.0/2.0)-3.0*b1*b2*LL2+b3/2.0);
    
    return als3;
  }
  
//! forward declaration of coefficient functions
template <int Iz,int Nf,int Ord>
constexpr double C();

//! instantiate the table of C
#define DEFC(IZ,NF,ORD,VAL)			\
  template<>					\
  constexpr double C<evol::IZ,evol::NF,ORD>()	\
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
  constexpr double AdimNum<evol::IZ>()		\
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
  constexpr double AdimDen<evol::NF>()		\
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
  //! compute the evolution up to 1/a2 of a Z
  template<int Iz,int Nf,int Ord=3>
  double evolution_to_RIp(double ainv,double mu_2)
  {
    double mu2[2]={mu_2,sqr(ainv)};
    double cmu2[2];
    
    for(size_t i=0;i<2;i++)
      {
	double al=alphas<Nf>(mu2[i])/(4*M_PI);
	cmu2[i]=1.0;
	if(Ord>=1) cmu2[i]+=al*C<Iz,Nf,1>();
	if(Ord>=2) cmu2[i]+=al*al*C<Iz,Nf,2>();
	if(Ord>=3) cmu2[i]+=al*al*al*C<Iz,Nf,3>();
	
	cmu2[i]*=pow(al,Adim<Iz,Nf>());
      }
    
    return cmu2[0]/cmu2[1];
  }
  
  //! parse Ord
  template<int Iz,int Nf>
  inline double evolution_to_RIp(int evo_ord,double ainv,double mu_2)
  {
    double e=0.0;
    
    switch(evo_ord)
      {
      case 0:e=evolution_to_RIp<Iz,Nf,0>(ainv,mu_2);break;
      case 1:e=evolution_to_RIp<Iz,Nf,1>(ainv,mu_2);break;
      case 2:e=evolution_to_RIp<Iz,Nf,2>(ainv,mu_2);break;
      case 3:e=evolution_to_RIp<Iz,Nf,3>(ainv,mu_2);break;
      default:
	CRASH("Not defined for ord %d",evo_ord);
      }
    
    return e;
  }
  
  //! parse Nf
  template<int Iz>
  inline double evolution_to_RIp(evol::Nf_t Nf,int ord,double ainv,double mu_2)
  {
    double e=0.0;
    
    switch(Nf)
      {
      case evol::NF0:e=evolution_to_RIp<Iz,evol::NF0>(ord,ainv,mu_2);break;
      case evol::NF2:e=evolution_to_RIp<Iz,evol::NF2>(ord,ainv,mu_2);break;
      case evol::NF4:e=evolution_to_RIp<Iz,evol::NF4>(ord,ainv,mu_2);break;
      default:
	CRASH("Not defined for Nf2 %d",Nf);
      }
    
    return e;
  }
  
  //! parse Iz
  inline double evolution_Zbil_to_RIp(ibil_t iZ,evol::Nf_t Nf,int ord,double ainv,double mu_2)
  {
    double e=0.0;
    
    switch(iZ)
      {
      case ibil_t::iS:
	//falltrhough
      case ibil_t::iP:
	e=evolution_to_RIp<evol::S>(Nf,ord,ainv,mu_2);
	break;
      case ibil_t::iT:
	e=evolution_to_RIp<evol::T>(Nf,ord,ainv,mu_2);
	break;
      case ibil_t::iV:
	//falltrhough
      case ibil_t::iA:
	e=1.0;
	break;
      default:
	CRASH("Not defined for iZ %d",iZ);
      }
    
    return e;
  }
  
  //! parse Zq
  inline double evolution_Zq_to_RIp(evol::Nf_t Nf,int ord,double ainv,double mu_2)
  {
    return evolution_to_RIp<evol::Q>(Nf,ord,ainv,mu_2);
  }
}

#endif
