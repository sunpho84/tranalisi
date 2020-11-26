#include <tranalisi.hpp>

#include "perens.hpp"

constexpr int q1=+2,q2=-1;

int main()
{
  set_njacks(10);
  
  perens_t ens;
  
  
  // for(size_t i=0;i<index3pts.max();i++)
  //   {
  //     const size_t isl=index3pts(i)[2];
      
  //     threePts[i].ave_err().write(combine("plots/threeptsRaw/%s.xmg",index3pts.descr(i).c_str()));
  //     const djvec_t tsymm=threePts[i].symmetrized(+1)/norm[isl];
  //     const djvec_t tasymm=threePts[i].symmetrized(-1)/norm[isl];
  //     tasymm.ave_err().write(combine("plots/threeptsRaw/antisymm_%s.xmg",index3pts.descr(i).c_str()));
  //     tsymm.ave_err().write(combine("plots/threeptsRaw/symm_%s.xmg",index3pts.descr(i).c_str()));
  //   }
  
  // for(size_t hAV=0;hAV<2;hAV++)
  //   for(size_t icomb=0;icomb<(size_t)f.ncomb;icomb++)
  //     for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
  // 	{
  // 	  const djvec_t out=((hAV==0)?load3ptsHA:load3ptsHV)(icomb,isl)/norm[isl];
  // 	  out.ave_err().write(combine("plots/threepts/H%s_comb%zu_isl%zu.xmg",havTag[hAV],icomb,isl));
	  
  // 	  const djvec_t outImpr=out-((hAV==0)?load3ptsHA:load3ptsHV)(0,isl)/norm[isl];
  // 	  outImpr.ave_err().write(combine("plots/threeptsImpr/H%s_comb%zu_isl%zu.xmg",havTag[hAV],icomb,isl));
	  
  // 	  for(size_t alpha=1; alpha<=2;alpha++)
  // 	    for(size_t r=1;r<=2;r++)
  // 	      {
  // 		const djvec_t out=load3ptsPol((HAV)hAV,icomb,isl,alpha,r)/norm[isl];
  // 		out.ave_err().write(combine("plots/threeptsPol/H%s_comb%zu_isl%zu_al%zu_r%zu.xmg",havTag[hAV],icomb,isl,alpha,r));
  // 	      }
  // 	}
  
  // const djvec_t with_without_sme_ratio=
  //   twoPts[index2pts(PP,is,il,1/*with*/,0)].symmetrized()/twoPts[index2pts(PP,is,il,0/*without*/,0)].symmetrized();
  // with_without_sme_ratio.ave_err().write("plots/with_without_PP.xmg");
  
  cout<<"m1: "<<ens.mu1<<" , q1: "<<q1<<endl;
  cout<<"m2: "<<ens.mu2<<" , q2: "<<q2<<endl;
  
  const djack_t aFp=(ens.mu1+ens.mu2)*ens.ZPloc/sqr(ens.mP);
  cout<<"a*fP: "<<aFp.ave_err()<<endl;
  
  auto g=[&ens](const size_t iAV,const size_t iX) -> djvec_t
	 {
	   return
	     ens.c3[ens.iC3({0,iAV,iX})]*q1+
	     ens.c3[ens.iC3({1,iAV,iX})]*q2;
	 };
  
  djvec_t Rtot[2],fA(ens.nX),fV(ens.nX);
  for(size_t iAV=0;iAV<2;iAV++)
    Rtot[iAV]=djvec_t(ens.nX);
  
  for(size_t iX=0;iX<ens.nX;iX++)
    {
      cout<<""<<ens.xG[iX].ave_err()<<endl;
      
      const djvec_t RA=g(HA,iX)/g(HA,0)-1.0;
      const djvec_t RV=g(HV,iX)/g(HA,0);
      
      for(size_t iIns=0;iIns<2;iIns++)
	for(size_t iAV=0;iAV<2;iAV++)
	  {
	    const djvec_t e=ens.c3[ens.iC3({iIns,iAV,iX})]/ens.norm;
	    const djack_t y=constant_fit(e,18,23,combine("plots/%s_ins%zu_x%zu.xmg",((iIns==HA)?"A":"V"),iAV,iX));
	  }
      
      Rtot[HA][iX]=constant_fit(RA,18,23,combine("plots/A%zu.xmg",iX));
      Rtot[HV][iX]=constant_fit(RV,18,23,combine("plots/V%zu.xmg",iX));
      
      fA[iX]=aFp*Rtot[HA][iX]/ens.xG[iX];
      fV[iX]=Rtot[HV][iX]/ens.kZ[iX];
    }
  
  grace_file_t fAplot("plots/fAFunX.xmg");
  grace_file_t fVplot("plots/fVFunX.xmg");
  
  fAplot.write_vec_ave_err(ens.xG.ave(),fA.ave_err());
  fVplot.write_vec_ave_err(ens.xG.ave(),fV.ave_err());
  
  return 0;
}
