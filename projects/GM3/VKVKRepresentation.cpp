#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <data.hpp>
#include <renoConstants.hpp>
#include <VKVKRepresentation.hpp>

double VKVKInfVol(const VKVKRepInfiniteVol<>& rep,
		  const ALaLuscherRepresentationInfVol<true>& aLaLusch,
		  const double& t)
{
  const double dualPart=
    rep.dualPartFun(t);
  
  const double LuschPart=
    aLaLusch(t);
  
  return
    dualPart+LuschPart;
}

pair<jack_t<VKVKRepFiniteVol<>>,
     jack_t<VKVKRepInfiniteVol<>>> fitVKVK(const RegoType& rego,
					   const int& nLevels)
{
  const djvec_t cP5P5=
    getAveForRego(0,nSources,idP5P5,REGO_TM);
  
  /// Charge factor of the correlator
  const double chargeFactor=
    5.0/9;
  
  djvec_t cVKVK=
    getAveForRego(0,nSources,idVKVK,rego)*sqr(Z[regoZId[rego]])*chargeFactor;
  cVKVK.ave_err().write("plots/corr_for_Rep"+regoTag[rego]+".xmg");
    
  jack_fit_t fitter;
  
  djvec_t pars(4);
  pars[0].fill_gauss((rego==REGO_TM)?1.05:1.20,0.04,235235);
  pars[1].fill_gauss(((rego==REGO_TM)?0.416:0.5)*aAve,0.05,7342);
  pars[2].fill_gauss(((rego==REGO_TM)?0.760:0.7900)*aAve,0.05,23423);
  pars[3].fill_gauss(((rego==REGO_TM)?28:26),0.5,32235);
  
  const size_t iRDual=
    fitter.add_fit_par_limits(pars[0],"RDual",pars[0].ave(),pars[0].err(), 0.2,2.0);
  const size_t iEThr=
    fitter.add_fit_par_limits(pars[1],"EThr",pars[1].ave(),pars[1].err(), 0.01,1.0);
  const size_t iMRho=
    fitter.add_fit_par_limits(pars[2],"MRho",pars[2].ave(),pars[2].err(), 0.01,1.0);
  const size_t iG2=
    fitter.add_fit_par_limits(pars[3],"g2",pars[3].ave(),pars[3].err(), 20.0,35.0);
  
  //fitter.fix_par(iEThr);
  //fitter.fix_par(iMRho);
  //fitter.fix_par(iG2);
  
  const djack_t aMPi=
    constant_fit(effective_mass(cP5P5),tMinP5P5[REGO_TM],tMaxP5P5[REGO_TM],"plots/eff_mass_P5P5"+regoTag[REGO_TM]+"_for_rep.xmg");
  
  // fit_debug=true;
  
  const size_t tMinFit=
    (rego==REGO_TM)?
    5:
    6;
  
  ALaLuscherRepresentationCached<true> cachedLuschRepFinder(nLevels);
  
  size_t tFit=
    tMinFit;
  
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
  
  cout<<"Fitting in the range: ["<<tMinFit<<":"<<tFit-1<<"]"<<endl;
  
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
  const jack_t<VKVKRepFiniteVol<>> rep(LuschRep,rDual,eThr,aMRho);
  
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
  
  cVKVKFull.ave_err().write("plots/VKVK_"+regoTag[rego]+"_representation.xmg");
  
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
  
  const djack_t mRho=aMRho/(*a);
  
  cout<<"mPi: "<<aMPi.ave_err()<<endl;
  cout<<"rDual: "<<rDual.ave_err()<<endl;
  cout<<"eThr: "<<eThr.ave_err()<<endl;
  cout<<"amRho: "<<aMRho.ave_err()<<endl;
  cout<<"mRho: "<<mRho.ave_err()<<endl;
  cout<<"g: "<<g.ave_err()<<endl;
  
  grace_file_t plotFit("plots/aLaLuscherFit"+regoTag[rego]+".xmg");
  plotFit.write_vec_ave_err(effective_mass(cVKVK).ave_err());
  plotFit.set_no_line();
  plotFit.write_polygon(effMassFun,tMinFit,tMaxFit);
  plotFit.write_polygon(effMassFun,tMaxFit,TH-1,grace::GREEN4);
  
  djvec_t cVKVKInfVol=
    cVKVK;
  
  jack_t<ALaLuscherRepresentationInfVol<true>> aLaLusch(aMPi,g2,aMRho);
  jack_t<VKVKRepInfiniteVol<>> repInfVol(aLaLusch,rDual,eThr,aMRho);
  

  // for(size_t t=0;t<TH;t++)
  //   cVKVKInfVol[t]=jackCall(VKVKInfVol,rep,aLaLusch,t);
  // cVKVKInfVol.ave_err().write("plots/cVKVK_"+regoTag[rego]+"_InfVol.xmg");
  // effective_mass(cVKVKInfVol,TH,0).ave_err().write("plots/eff_mass_VKVK_InfVol_"+regoTag[rego]+".xmg");
 
  return
    {rep,repInfVol};
}
