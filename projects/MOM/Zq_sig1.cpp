#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <Zq.hpp>

double compute_Zq_sig1(const prop_t &prop_inv,size_t imom)
{
  double Zq;
  
  p_t ptilde=glb_moms[imom].p(L).tilde();
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      Zq=0.0;
      for(size_t mu=0;mu<NDIM;mu++)
	if(fabs(ptilde[mu])>1e-10)
	  Zq+=
	    (prop_inv*Gamma[igmu[mu]]).trace().imag()/
	    (12.0*ptilde[mu]*V*glb_moms[imom].Np());
    }
  
  return Zq;
}

djack_t compute_Zq_sig1(const jprop_t &jprop_inv,size_t imom)
{
  djack_t Zq;
  
  p_t ptilde=glb_moms[imom].p(L).tilde();
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      Zq[ijack]=0.0;
      for(size_t mu=0;mu<NDIM;mu++)
	if(fabs(ptilde[mu])>1e-10)
	  Zq[ijack]+=
	    (jprop_inv[ijack]*Gamma[igmu[mu]]).trace().imag()/
	    (12.0*ptilde[mu]*V*glb_moms[imom].Np());
    }
  
  return Zq;
}

djvec_t compute_Zq_sig1(const vjprop_t &jprop_inv)
{
  djvec_t Zq(glb_moms.size()); //!< Z of quark field
  
#pragma omp parallel for
  for(size_t imom=0;imom<glb_moms.size();imom++) Zq[imom]=compute_Zq_sig1(jprop_inv[imom],imom);
  
  return Zq;
}
