#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <MOM2/Zq.hpp>

#include <MOM2/perens.hpp>

double perens_t::compute_Zq(const qprop_t &prop_inv,const size_t glb_mom)
{
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  const double pt2=ptilde.norm2();
  const qprop_t pslash=qua_slash(ptilde);
  
  const double Zq=(prop_inv*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}

djack_t perens_t::compute_Zq(const jqprop_t &jprop_inv,const size_t glb_mom)
{
  djack_t Zq;
  
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  const double pt2=ptilde.norm2();
  const qprop_t pslash=qua_slash(ptilde);
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    Zq[ijack]=(jprop_inv[ijack]*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}
