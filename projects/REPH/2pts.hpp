#ifndef _2PTS_HPP
#define _2PTS_HPP

#include <REPH/reph.hpp>

djvec_t perens_t::load2ptsPP(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt)
{
  const index_t ind({{"iks",nMass},
		     {"ikt",nMass},
		     {"moms",nMoms},
		     {"momt",nMoms},
		     {"gamma",1},
		     {"reim",2}});
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,0,0});
  
  return read_djvec("jacks/oPPo-ss",T,i).symmetrized()/spatVol;
}

djvec_t perens_t::load2ptsAP(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iGamma)
{
  const index_t ind({{"iks",nMass},
		     {"ikt",nMass},
		     {"moms",nMoms},
		     {"momt",nMoms},
		     {"gamma",4},
		     {"reim",2}});
  
  const size_t iReIm=(iGamma==0)?RE:IM;
  const int par=(iGamma==0)?-1:+1;
  const int sign=(iGamma==0)?-1:+1;
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,iGamma,iReIm});
  
  return sign*read_djvec("jacks/oAmuPo-ss",T,i).symmetrized(par)/spatVol;
}

#endif
