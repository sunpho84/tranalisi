#ifndef _2PTS_HPP
#define _2PTS_HPP

#include <REPH/permes_combo.hpp>

std::string kinTag(const size_t iMom1,const size_t iMom2)
{
  return combine("imom1_%zu_imom2_%zu",iMom1,iMom2);
}

template <typename TV>
djvec_t permes_combo_t<TV>::load2ptsPP(const size_t iMoms,const size_t iMomt)
{
  const index_t ind({{"iks",ens.nMass},
		     {"ikt",ens.nMass},
		     {"moms",ens.nMoms},
		     {"momt",ens.nMoms},
		     {"gamma",1},
		     {"reim",2}});
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,0,0});
  
  return ens.loadCorr(combine("%s/jacks/oPPo-ss",ens.dirPath.c_str()),i,"2").symmetrized()/ens.spatVol;
}

template <typename TV>
djvec_t permes_combo_t<TV>::load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma)
{
  const index_t ind({{"iks",ens.nMass},
		     {"ikt",ens.nMass},
		     {"moms",ens.nMoms},
		     {"momt",ens.nMoms},
		     {"gamma",4},
		     {"reim",2}});
  
  const size_t iReIm=(iGamma==0)?RE:IM;
  const int par=(iGamma==0)?-1:+1;
  const int sign=(iGamma==0)?-1:+1;
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,iGamma,iReIm});
  
  return sign*ens.loadCorr(combine("%s/jacks/oAmuPo-ss",ens.dirPath.c_str()),i,"2").symmetrized(par)/ens.spatVol;
}

template <typename TV>
void permes_combo_t<TV>::load2pts(const bool forceLoad)
{
  const string dataPath=milledPath()+"/data2pts.dat";
  
  if(forceLoad or not file_exists(dataPath))
    {
      //cout<<"Loading 2pts correlators from scratch"<<endl;
      
      const string plotsPath=combine("%s/plots/%s/2pts_corr/",ens.dirPath.c_str(),mesTag.c_str());
      
      for(auto dir : {"PP","A0P","AP3"})
	{
	  const std::string path=combine("%s/%s/",plotsPath.c_str(),dir);
	  if(not dir_exists(path))
	    mkdir(path);
	}
      
      for(size_t iMesKin=0;iMesKin<ens.indMesKin.max();iMesKin++)
	{
	  const vector<size_t> c=ens.indMesKin(iMesKin);
	  const size_t iMom1=c[0];
	  const size_t iMom2=c[1];
	  
	  const std::string kinTag=combine("imom1_%zu_imom2_%zu",iMom1,iMom2);
	  
	  const bool loadSymm=false;
	  
	  if(loadSymm)
	    corrPP[iMesKin]=(load2ptsPP(iMom1,iMom2)+load2ptsPP(iMom2,iMom1))/2.0;
	  else
	    corrPP[iMesKin]=load2ptsPP(iMom1,iMom2);
	  corrPP[iMesKin].ave_err().write(combine("%s/PP/%s.xmg",plotsPath.c_str(),kinTag.c_str()));
	  
	  if(loadSymm)
	    corrA0P[iMesKin]=(load2ptsAP(iMom1,iMom2,0)+load2ptsAP(iMom2,iMom1,0))/2.0;
	  else
	    corrA0P[iMesKin]=load2ptsAP(iMom1,iMom2,0);
	  corrA0P[iMesKin].ave_err().write(combine("%s/A0P/%s.xmg",plotsPath.c_str(),kinTag.c_str()));
	  
	  if(loadSymm)
	    corrA3P[iMesKin]=(load2ptsAP(iMom1,iMom2,3)-load2ptsAP(iMom2,iMom1,3))/2.0;
	  else
	    corrA3P[iMesKin]=load2ptsAP(iMom1,iMom2,3);
	  corrA3P[iMesKin].ave_err().write(combine("%s/AP3/%s.xmg",plotsPath.c_str(),kinTag.c_str()));
	}
      
      //cout<<"Storing compacted 2pts"<<endl;
      raw_file_t dataFile(dataPath,"w");
      for(auto& c : {&corrPP,&corrA0P,&corrA3P})
	dataFile.bin_write(*c);
    }
  else
    {
      //cout<<"Loading compacted 2pts"<<endl;
      
      raw_file_t dataFile(dataPath,"r");
      for(auto& c : {&corrPP,&corrA0P,&corrA3P})
	dataFile.bin_read(*c);
    }
  
  for(size_t iMesKin=0;iMesKin<ens.indMesKin.max();iMesKin++)
    {
      eEff[iMesKin]=effective_mass(corrPP[iMesKin]);
      zEff[iMesKin]=sqrt(effective_squared_coupling(corrPP[iMesKin],eEff[iMesKin],ens.T/2));
    }
}

template <typename TV>
void permes_combo_t<TV>::plotDispRel(const string& mesPlotsPath) const
{
  grace_file_t dispRel(combine("%s/dispRel.xmg",mesPlotsPath.c_str()));
  dispRel.write_vec_ave_err(ens.pMes,E.ave_err());
  
  const djack_t &M=E[0];
  
  dispRel.set_legend("Lattice");
  dispRel.write_polygon([M](double x){return latt_en_1D(M,x);},0,ens.pMesMax,grace::GREEN);
  dispRel.set_legend("Continuum");
  dispRel.write_polygon([M](double x){return cont_en_1D(M,x);},0,ens.pMesMax,grace::VIOLET);
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::choose2ptsTint(const string& mesPlotsPath,const bool forceRechoose)
{
  string tintPath=milledPath()+"/tint2pts.dat";
  
  if(file_exists(tintPath) and not forceRechoose)
    {
      //cout<<"Loading stored 2pts intervals"<<endl;
      raw_file_t(tintPath,"r").bin_read(tint2pts);
    }
  else
    {
      //cout<<"Choosing tint from scratch"<<endl;
      
      const string rangeDir=mesPlotsPath+"/range";
      mkdir(rangeDir);
      
      grace_file_t fitTxt(rangeDir+"/range.txt");
      
      for(size_t iMesKin=0;iMesKin<ens.indMesKin.max();iMesKin++)
	{
	  const vector<size_t> c=ens.indMesKin(iMesKin);
	  const size_t iMom1=c[0];
	  const size_t iMom2=c[1];
	  
	  const string rangePath=combine("%s/%s.xmg",rangeDir.c_str(),kinTag(iMom1,iMom2).c_str());
	  //cout<<rangePath<<endl;
	  
	  const djvec_t temp=smooth(forward_derivative(effective_mass(corrPP[iMesKin])),3);
	  
	  CompatibilityRangeFinder<> comp(temp,rangePath);
	  
	  //! Initial value of enlarging factor
	  double nSigEnl=2.0;
	  
	  const auto recipe=[&comp,&nSigEnl](const bool verb=false)
			     {
			       double nSigSel=1;
			       comp
				 .setVerbose(verb)
				 .setRangeToConsiderByFraction(1.0/12)
				 .selectClosestCompatiblePointWithinNsigma(0,nSigSel)
				 .plotSelected()
				 .selectEnlargingError(2,nSigEnl)
				 .plotSelected()
				 .mergeSelectionByDistanceSize()
				 .plotSelected()
				 .selectLargestRange()
				 .plotSelected()
				 .setRangeToConsiderByFraction(0.0)
				 .extendRightWithCompatiblePoints(nSigEnl);
			     };
	  
	  recipe();
	  
	  //! Minimal length
	  const size_t lengthMin=ens.T/6;
	  
	  // //! Trim some points if allowed
	  // const size_t nDesel=2;
	  // if(comp.getSelectionRange(0).size()>=lengthMin+nDesel)
	  //   comp
	  //   .deSelectLeftermostPoints(nDesel)
	  //   .plotSelected();
	  
	  //! Function to check if redo
	  const auto checkReDo=[&comp,&lengthMin](){return comp.getSelectionRange(0).size()<lengthMin;};
	  
	  //! If needed to redo, verbose
	  if(checkReDo())
	    recipe(false);
	  
	  // Extend
	  int nReDo=4;
	  
	  while(checkReDo() and nReDo>0)
	    {
	      nSigEnl+=1.0;
	      
	      //cout<<"WARNING, compatibility range too short, "<<comp.getSelectionRange(0)<<" when at least "<<lengthMin<<" needed, enlarging right of "<<nSigEnl<<endl;
	      
	      comp
		.setRangeToConsiderByFraction(0.0)
		.extendRightWithCompatiblePoints(nSigEnl)
		.plotSelected()
		.selectLargestRange()
		.plotSelected();
	      
	      nReDo--;
	    }
	  if(nReDo<=0)
	    cout<<"WARNING, unable to satisfy the check!"<<endl;
	  
	  const Range& fitRange=tint2pts[iMesKin]=comp.getSelectionRange(0);
	  
	  //cout<<iMesKin<<" "<<ens.indMesKin.descr(iMesKin)<<": range "<<fitRange<<endl;
	  
	  const size_t &b=fitRange.begin,&e=fitRange.end;
	  fitTxt.write_ave_err(iMesKin,{(b+e)/2.0,(e-b)/2.0});
	}
      
      //cout<<"Storing chosen 2pts intervals"<<endl;
      raw_file_t(tintPath,"w").bin_write(tint2pts);
    }
  
  return *this;
}

template <typename TV>
void permes_combo_t<TV>::computeAxialPseudoCouplings(const string& mesPlotsPath)
{
  string apePath=milledPath()+"/APE.dat";
  
  const bool load=file_exists(apePath);
  
  if(not load)
    {
      const size_t& T=ens.T;
      
      grace_file_t ZPPlot(combine("%s/ZP.xmg",mesPlotsPath.c_str()));
      grace_file_t ZAPlot(combine("%s/ZA.xmg",mesPlotsPath.c_str()));
      
      for(auto dir : {"CPP_fit","CA0P_fit","CA3P_fit","CPP_A0P_fit","range"})
	{
	  const std::string path=combine("%s/%s/",mesPlotsPath.c_str(),dir);
	  if(not dir_exists(path))
	    mkdir(path);
	}
      
      for(size_t iMesKin=0;iMesKin<ens.indMesKin.max();iMesKin++)
	{
	  const vector<size_t> c=ens.indMesKin(iMesKin);
	  const size_t iMom1=c[0];
	  const size_t iMom2=c[1];
	  
	  const Range& tint=tint2pts[iMesKin];
	  const size_t tMin=tint.begin;
	  const size_t tMax=tint.end+1;
	  
	  if(iMom1!=iMom2)
	    {
	      const double P=ens.pMes[iMesKin];
	      
	      //perform a preliminary fit
	      djack_t fit_E,ZPP,ZA0P,ZA3P=1.0;
	      two_pts_fit(ZPP, fit_E,corrPP[iMesKin], T/2,tMin,tMax,combine("%s/CPP_fit/prel_%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),"",+1);
	      two_pts_fit(ZA0P,fit_E,corrA0P[iMesKin],T/2,tMin,tMax-1,combine("%s/CA0P_fit/prel_%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),"",-1);
	      two_pts_fit(ZA3P,fit_E,corrA3P[iMesKin],T/2,tMin,tMax,combine("%s/CA3P_fit/prel_%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),"",+1);
	      const djack_t fit_ZP=sqrt(ZPP);
	      const djack_t M=sqrt(sqr(fit_E)-sqr(P));
	      const djack_t ZA0=ZA0P/fit_ZP,ZA3=ZA3P/fit_ZP;
	      const djack_t fit_ZA=(ZA0*fit_E-ZA3*P)/M;
	      
	      //parameters to fit
	      minimizer_pars_t pars;
	      pars.add("E",fit_E.ave(),fit_E.err());
	      pars.add("ZA",fit_ZA.ave(),fit_ZA.err());
	      pars.add("ZP",fit_ZP.ave(),fit_ZP.err());
	      
	      //! fit for real
	      size_t iel=0;
	      auto x=vector_up_to<double>(corrPP[iMesKin].size());
	      multi_ch2_t<djvec_t> two_pts_fit_obj({x,x,x},{tMin,tMin,tMin},{tMax,tMax-1,tMax},{corrPP[iMesKin],corrA0P[iMesKin],corrA3P[iMesKin]},
						   {two_pts_corr_fun_t(T/2,+1),two_pts_corr_fun_t(T/2,-1),two_pts_corr_fun_t(T/2,+1)},
						   [&P](const vector<double> &p,size_t icontr)
						   {
						     switch(icontr)
						       {
						       case 0:return vector<double>({p[2]*p[2],p[0]});break;
						       case 1:return vector<double>({p[1]*p[2]*p[0],p[0]});break;
						       case 2:return vector<double>({p[1]*p[2]*P,p[0]});break;
						       default: CRASH("Unknown contr %zu",icontr);return p;
						       }
						   },iel);
	      
	      minimizer_t minimizer(two_pts_fit_obj,pars);
	      
	      for(iel=0;iel<=njacks;iel++)
		{
		  //minimize and print the result
		  vector<double> par_min=minimizer.minimize();
		  E[iMesKin][iel]=par_min[0];
		  ZA[iMesKin][iel]=par_min[1];
		  ZP[iMesKin][iel]=par_min[2];
		}
	      
	      write_constant_fit_plot(combine("%s/CPP_fit/%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),tMin,tMax,E[iMesKin],effective_mass(corrPP[iMesKin],T/2,+1));
	      write_constant_fit_plot(combine("%s/CA0P_fit/%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),tMin,tMax,E[iMesKin],effective_mass(corrA0P[iMesKin],T/2,-1));
	      write_constant_fit_plot(combine("%s/CA3P_fit/%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),tMin,tMax,E[iMesKin],effective_mass(corrA3P[iMesKin],T/2,+1));
	    }
	  else
	    {
	      two_pts_SL_fit(ZP[iMesKin],ZA[iMesKin],E[iMesKin],corrA0P[iMesKin],corrPP[iMesKin],T/2,tMin,tMax-1,combine("%s/CPP_A0P_fit/%s.xmg",mesPlotsPath.c_str(),kinTag(iMom1,iMom2).c_str()),-1,+1);
	      ZA[iMesKin]/=E[iMesKin];
	    }
	  
	  ZAPlot.write_ave_err(E[iMesKin].ave(),ZA[iMesKin].ave_err());
	  ZPPlot.write_ave_err(E[iMesKin].ave(),ZP[iMesKin].ave_err());
	}
    }
  
  
  //cout<<(load?"Loading":"Storing")<<" Axial and Pseudoscalar coupling and Energies"<<endl;
  
  raw_file_t file(apePath,(load?"r":"w"));
  
  for(auto& ape: {&ZA,&ZP,&E})
    if(load)
      file.bin_read(*ape);
    else
      file.bin_write(*ape);
  
  for(size_t iMesKin=0;iMesKin<ens.indMesKin.max();iMesKin++)
    {
      fP[iMesKin]=ZP[iMesKin]*(ens.mass[iMs]+ens.mass[iMt])/sqr(E[iMesKin]);
      fPbare[iMesKin]=ZA[iMesKin];
    }
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::fit2pts(const char* fitTag,const bool forceRechoose)
{
  const string mesPlotsPath=combine("%s/plots/%s/2pts_fit_%s",ens.dirPath.c_str(),mesTag.c_str(),fitTag);
  mkdir(mesPlotsPath);
  
  choose2ptsTint(mesPlotsPath,forceRechoose);
  
  computeAxialPseudoCouplings(mesPlotsPath);
  plotDispRel(mesPlotsPath);
  
  return *this;
}

#endif
