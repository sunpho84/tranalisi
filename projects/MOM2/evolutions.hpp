#ifndef _EVOLUTIONS_HPP
#define _EVOLUTIONS_HPP

#include <tranalisi.hpp>

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
}

#endif
