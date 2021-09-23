#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <data.hpp>
#include <renoConstants.hpp>
#include <VKVKRepresentation.hpp>

jack_t<VKVKRep> fitVKVK(const int nLevels,
			const size_t tMin)
{
  const djvec_t cP5P5=
    getAveForRego(0,nSources,idP5P5,REGO_TM);
  
  /// Charge factor of the correlator
  const double chargeFactor=
    5.0/9;
  
  djvec_t cVKVK=
    getAveForRego(0,nSources,idVKVK,REGO_TM)*sqr(Z[regoZId[REGO_TM]])*chargeFactor;
  cVKVK.ave_err().write("plots/corr_for_Rep.xmg");
  cVKVK[0]=cVKVK[1];
  
  jack_fit_t fitter;
  
  djvec_t pars(4);
  pars[0].fill_gauss(1.15,0.01,235235);
  pars[1].fill_gauss(0.13,0.01,7342);
  pars[2].fill_gauss(0.22,0.01,23423);
  pars[3].fill_gauss(28.1,0.01,32235);
  
  const size_t iRDual=
    fitter.add_fit_par_limits(pars[0],"RDual",pars[0].ave(),pars[0].err(), 0.2,2.0);
  const size_t iEThr=
    fitter.add_fit_par_limits(pars[1],"EThr",pars[1].ave(),pars[1].err(), 0.01,1.0);
  const size_t iMRho=
    fitter.add_fit_par_limits(pars[2],"MRho",pars[2].ave(),pars[2].err(), 0.01,1.0);
  const size_t iG2=
    fitter.add_fit_par_limits(pars[3],"g2",pars[3].ave(),pars[3].err(), 25.0,35.0);
  
  const djack_t aMPi=
    constant_fit(effective_mass(cP5P5),tMinP5P5[0],tMaxP5P5[0]);
  
  //fit_debug=true;
  
  const size_t tMinFit=
    5;
  
  size_t tFit=
    tMinFit;
  
  ALaLuscherRepresentationCached<true> cachedLuschRepFinder(nLevels);
  
  while(tFit<THp1 and cVKVK[tFit].err()<0.1*fabs(cVKVK[tFit].ave()))
    {
      fitter.add_point(cVKVK[tFit],
		       [&aMPiAll=aMPi,
			&cachedLuschRepFinder,
			iRDual,
			iEThr,
			iMRho,
			iG2,
			tFit](const vector<double>& p,
			      const size_t ijack)->double
		       {
			 const double aMPi=
			   aMPiAll[ijack];
			 
			 const double& rDual=
			   p[iRDual];
			 const double eThr=
			   p[iEThr];
			 const double amRho=
			   p[iMRho];
			 const double& g2=
			   p[iG2];
			 
			 VKVKRep rep(cachedLuschRepFinder(aMPi,L,amRho,g2),rDual,eThr,amRho);
			 
			 return
			   rep(tFit)();
		       });
      tFit++;
    };
  
  const size_t tMaxFit=
    tFit;
  
  const djvec_t oldPars=
    pars;
  
  fitter.fit();
  
  const djack_t& rDual=
    pars[iRDual];
  const djack_t eThr=
    pars[iEThr];
  const djack_t aMRho=
    pars[iMRho];
  const djack_t g2=
    pars[iG2];
  const djack_t g=
    sqrt(g2);
  
  vector<vector<double>> parsPerJack(njacks+1,vector<double>(4));
  for(int i=0;i<4;i++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      parsPerJack[ijack][i]=pars[i][ijack];
  
  const jack_t<ALaLuscherRepresentationCalculator> jLuschRepConstructor(aMPi,L,g2,aMRho);
  const jack_t<ALaLuscherRepresentation<true>> LuschRep=
    jLuschRepConstructor(nLevels);
  const jack_t<VKVKRep> rep(LuschRep,rDual,eThr,aMRho);
  
  djvec_t cVKVKDual(THp1);
  djvec_t cVKVKLuscher(THp1);
  djvec_t cVKVKFull(THp1);
  for(size_t t=0;t<=TH;t++)
    {
      const auto r=
	rep(t);
      
      for(size_t ijack=0;ijack<=njacks;ijack++)
	{
	  cVKVKDual[t][ijack]=r[ijack].dual;
	  cVKVKLuscher[t][ijack]=r[ijack].Lusch;
	}
    }
  
  cVKVKFull=
    cVKVKLuscher+cVKVKDual;
  
  djvec_t weights(nLevels),energy(nLevels);
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      const auto c=
	rep[ijack].LuschRep.coeffs;
      
      for(int ilev=0;ilev<nLevels;ilev++)
	{
	  weights[ilev][ijack]=c[ilev].weight;
	  energy[ilev][ijack]=c[ilev].energy;
	}
    }
  
  console<<"Levels (i,w,e):"<<endl;
  console<<"-------"<<endl;
  for(int ilev=0;ilev<nLevels;ilev++)
    console<<" "<<ilev<<" "<<weights[ilev].ave_err()<<" "<<energy[ilev].ave_err()<<endl;
  
  auto effMassFun=
    [&rep](const double& tPlot)
    {
      const double dTPlot=
	1;
      
      const djack_t y=
	rep(tPlot)();
      
      const djack_t yPlusDt=
	rep(tPlot+dTPlot)();
      
      const djack_t e=
	effective_mass(y,yPlusDt,tPlot,TH,1,0,dTPlot);
      
      return
	e;
    };
  
  cout<<"mPi: "<<aMPi.ave_err()<<endl;
  cout<<"rDual: "<<rDual.ave_err()<<endl;
  cout<<"eThr: "<<eThr.ave_err()<<endl;
  cout<<"amRho: "<<aMRho.ave_err()<<endl;
  cout<<"g: "<<g.ave_err()<<endl;
  
  grace_file_t plotFit("plots/aLaLuscherFit.xmg");
  plotFit.write_vec_ave_err(effective_mass(cVKVK).ave_err());
  plotFit.set_no_line();
  plotFit.write_polygon(effMassFun,tMinFit,tMaxFit);
  plotFit.write_polygon(effMassFun,tMaxFit,TH-1,grace::GREEN4);
  
  return
    rep;
}
