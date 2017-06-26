#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <Zq.hpp>

double compute_Zq(const prop_t &prop_inv,size_t imom)
{
  double Zq;
  
  p_t ptilde=imoms[imom].p(L).tilde();
  double pt2=ptilde.norm2();
  prop_t pslash=slash(ptilde);
  
  Zq=(prop_inv*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}

djack_t compute_Zq(const jprop_t &jprop_inv,size_t imom)
{
  djack_t Zq;
  
  p_t ptilde=imoms[imom].p(L).tilde();
  double pt2=ptilde.norm2();
  prop_t pslash=slash(ptilde);
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      prop_t prop_inv=get_from_jackknife(jprop_inv,ijack);
      Zq[ijack]=(prop_inv*pslash).trace().imag()/(12.0*pt2*V);
    }
  
  return Zq;
}

djvec_t compute_Zq(const vjprop_t &jprop_inv)
{
  djvec_t Zq(imoms.size()); //!< Z of quark field
  
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++) Zq[imom]=compute_Zq(jprop_inv[imom],imom);
  
  return Zq;
}
