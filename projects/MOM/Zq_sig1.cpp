#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <Zq.hpp>

double compute_Zq_sig1(const prop_t &prop_inv,const size_t glb_mom)
{
  double Zq;
  
  cout<<"Computing Zq_sig1 for mom "<<glb_mom<<endl;
  const p_t ptilde=glb_moms[glb_mom].p(L).tilde();
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      Zq=0.0;
      for(size_t mu=0;mu<NDIM;mu++)
	if(fabs(ptilde[mu])>1e-10)
	  Zq+=
	    (prop_inv*Gamma[igmu[mu]]).trace().imag()/
	    (12.0*ptilde[mu]*V*glb_moms[glb_mom].Np());
    }
  
  return Zq;
}

djack_t compute_Zq_sig1(const jprop_t &jprop_inv,const size_t glb_mom)
{
  djack_t Zq;
  
  cout<<"Computing Zq_sig1 for mom "<<glb_mom<<endl;
  const p_t ptilde=glb_moms[glb_mom].p(L).tilde();
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      Zq[ijack]=0.0;
      for(size_t mu=0;mu<NDIM;mu++)
	if(fabs(ptilde[mu])>1e-10)
	  Zq[ijack]+=
	    (jprop_inv[ijack]*Gamma[igmu[mu]]).trace().imag()/
	    (12.0*ptilde[mu]*V*glb_moms[glb_mom].Np());
    }
  
  return Zq;
}
