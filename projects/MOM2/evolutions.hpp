#ifndef _EVOLUTIONS_HPP
#define _EVOLUTIONS_HPP

#include <tranalisi.hpp>

namespace evol
{
  //! tag the bilinears, or Zq
  enum ev_bil_t{Q=0,P=1,S=1,T=2};
  
  //! tag the number of quarks
  enum Nf_t{NF0=0,NF2=2,NF4=4};
  
  //map Nf to Nf_t
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
}

#endif
