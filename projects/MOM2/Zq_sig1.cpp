#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <MOM2/perens.hpp>
#include <MOM2/Zq_sig1.hpp>

double perens_t::compute_Zq_sig1(const qprop_t &prop_inv,const size_t glb_mom)
{
  double Zq;
  
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  
  Zq=0.0;
  for(size_t mu=0;mu<NDIM;mu++)
    if(fabs(ptilde[mu])>1e-10)
      Zq+=
	(prop_inv*quaGamma[igmu[mu]]).trace().imag()/
	(12.0*ptilde[mu]*V*all_moms[glb_mom].Np());
  
  return Zq;
}

djack_t perens_t::compute_Zq_sig1(const jqprop_t &jprop_inv,const size_t glb_mom)
{
  djack_t Zq;
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    Zq[ijack]=compute_Zq_sig1(jprop_inv[ijack],glb_mom);
  
  return Zq;
}
