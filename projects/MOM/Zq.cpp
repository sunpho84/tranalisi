#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <Zq.hpp>

double compute_Zq(const prop_t &prop_inv,const size_t glb_mom)
{
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  const double pt2=ptilde.norm2();
  const prop_t pslash=slash(ptilde);
  
  const double Zq=(prop_inv*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}

djack_t compute_Zq(const jprop_t &jprop_inv,const size_t glb_mom)
{
  djack_t Zq;
  
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  const double pt2=ptilde.norm2();
  const prop_t pslash=slash(ptilde);
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    Zq[ijack]=(jprop_inv[ijack]*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}
