#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <Zq.hpp>

djvec_t compute_sig3(const vjprop_t &jprop_inv)
{
  djvec_t sigma3(imoms.size());
  
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      sigma3[imom][ijack]=
	(jprop_inv[imom][ijack]*Gamma[5]).trace().imag()/
	(12.0*V);
  
  return sigma3;
}
