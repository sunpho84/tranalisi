#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <Zbil.hpp>

vector<djvec_t> compute_proj_bil(const vjprop_t &jprop_inv1,vector<jverts_t> &jverts,const vjprop_t &jprop_inv2)
{
  const size_t iZbil_of_iG[nGamma]={0,1,1,1,1,2,3,3,3,3,4,4,4,4,4,4};
  const double Zdeg[nZbil]={1,4,1,4,6};
  vector<djvec_t> pr(nZbil,djvec_t(imoms.size(),0.0)); //!< Five summed projectors
  
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      {
	prop_t prop_inv1=get_from_jackknife(jprop_inv1[imom],ijack);
	prop_t prop_inv2=get_from_jackknife(jprop_inv2[imom],ijack);
	
	for(size_t iG=0;iG<nGamma;iG++)
	  {
	    size_t iZbil=iZbil_of_iG[iG];
	    prop_t vert=get_from_jackknife(jverts[imom][iG],ijack) ;
	    prop_t amp_vert=prop_inv1*vert*Gamma[5]*prop_inv2.adjoint()*Gamma[5];
	    pr[iZbil][imom][ijack]+=(amp_vert*Gamma[iG].adjoint()).trace().real()/(12.0*Zdeg[iZbil]);
	  }
      }
  
  return pr;
}

