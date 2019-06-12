#ifndef _3PTS_HPP
#define _3PTS_HPP

#include <REPH/permes_combo.hpp>

djvec_t permes_combo_t::load3pts(const size_t iVA,const size_t iMoms,const size_t iMomt,const size_t iMom0)
{
  const size_t& T=ens.T;
  const size_t& L=ens.L;
  const size_t& nMass=ens.nMass;
  const size_t& nMoms=ens.nMoms;
  
  const double s[2][2][2]={{{-1,+1},{-1,-1}},
			   {{-1,-1},{+1,-1}}};
  
  djvec_t corr(T);
  corr=0.0;
  
  const index_t ind({{"iks",nMass},
		     {"ikt",nMass},
		     {"moms",nMoms},
		     {"momt",nMoms},
		     {"mom0",nMoms},
		     {"pol",2},
		     {"gamma",4},
		     {"reim",2}});
  
  const index_t ind_ave({{"iks",nMass},
			 {"ikt",nMass},
			 {"moms",nMoms},
			 {"momt",nMoms},
			 {"mom0",nMoms}});
  
  const size_t iReIm=(iVA==0)?1:0;
  const int par=(iVA==0)?-1:+1;
  for(size_t iPol=0;iPol<2;iPol++)
    for(size_t iGamma=1;iGamma<=2;iGamma++)
      {
	const size_t i=ind({iMs,iMt,iMoms,iMomt,iMom0,iPol,iGamma,iReIm});
	const djvec_t contr=read_djvec(combine("jacks/o%smuGPo-gs",VA_tag[iVA]),T,i);
	corr+=contr*s[iVA][iPol][iGamma-1];
	
	contr.ave_err().write(combine("plots/o%smuGPo-gs_%s.xmg",VA_tag[iVA],ind.descr(i).c_str()));
      }
  corr/=2*sqrt(2);
  
  corr/=-(L*L*L);
  
  corr.symmetrize(par);
  const size_t iave=ind_ave({iMs,iMt,iMoms,iMomt,iMom0});
  corr.ave_err().write(combine("plots/o%smuGPo-gs_%s.xmg",VA_tag[iVA],ind_ave.descr(iave).c_str()));
  
  return corr;
}

#endif
