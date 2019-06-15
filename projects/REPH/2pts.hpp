#ifndef _2PTS_HPP
#define _2PTS_HPP

#include <REPH/permes_combo.hpp>

std::string kinTag(const size_t iMom1,const size_t iMom2)
{
  return combine("imom1_%zu_imom2_%zu",iMom1,iMom2);
}

djvec_t permes_combo_t::load2ptsPP(const size_t iMoms,const size_t iMomt)
{
  const index_t ind({{"iks",ens.nMass},
		     {"ikt",ens.nMass},
		     {"moms",ens.nMoms},
		     {"momt",ens.nMoms},
		     {"gamma",1},
		     {"reim",2}});
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,0,0});
  
  return read_djvec(combine("%s/jacks/oPPo-ss",ens.dirPath.c_str()),ens.T,i).symmetrized()/ens.spatVol;
}

djvec_t permes_combo_t::load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma)
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
  
  return sign*read_djvec(combine("%s/jacks/oAmuPo-ss",ens.dirPath.c_str()),ens.T,i).symmetrized(par)/ens.spatVol;
}

void permes_combo_t::load2pts()
{
  cout<<"Loading 2pts correlators"<<endl;
  
  const string plotsPath=combine("%s/plots/%s/2pts_corr/",ens.dirPath.c_str(),mesComboTag.c_str());
  
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
      
      corrPP[iMesKin]=load2ptsPP(iMom1,iMom2);
      corrPP[iMesKin].ave_err().write(combine("%s/PP/%s.xmg",plotsPath.c_str(),kinTag.c_str()));
      eEff[iMesKin]=effective_mass(corrPP[iMesKin]);
      
      corrA0P[iMesKin]=load2ptsAP(iMom1,iMom2,0);
      corrA0P[iMesKin].ave_err().write(combine("%s/A0P/%s.xmg",plotsPath.c_str(),kinTag.c_str()));
      
      corrA3P[iMesKin]=load2ptsAP(iMom1,iMom2,3);
      corrA3P[iMesKin].ave_err().write(combine("%s/AP3/%s.xmg",plotsPath.c_str(),kinTag.c_str()));
    }
}

void permes_combo_t::plotDispRel() const
{
  grace_file_t dispRel(combine("%s/dispRel.xmg",mesPlotsPath.c_str()));
  dispRel.write_vec_ave_err(ens.pMes,E.ave_err());
  
  const djack_t &M=E[0];
  
  dispRel.set_legend("Lattice");
  dispRel.write_polygon([M](double x){return latt_en_1D(M,x);},0,ens.pMesMax,grace::GREEN);
  dispRel.set_legend("Continuum");
  dispRel.write_polygon([M](double x){return cont_en_1D(M,x);},0,ens.pMesMax,grace::VIOLET);
}

void permes_combo_t::computeAxialPseudoCouplings()
{
  const size_t& T=ens.T;
  
  const size_t tMin=12,tMax=T/2+1;
  
  grace_file_t ZPPlot(combine("%s/ZP.xmg",mesPlotsPath.c_str()));
  grace_file_t ZAPlot(combine("%s/ZA.xmg",mesPlotsPath.c_str()));
  
  for(auto dir : {"CPP_fit","CA0P_fit","CA3P_fit","CPP_A0P_fit"})
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
      
      fP[iMesKin]=ZP[iMesKin]*(ens.mass[iMs]+ens.mass[iMt])/sqr(E[iMesKin]);
      fPbare[iMesKin]=ZA[iMesKin];
    }
}

permes_combo_t& permes_combo_t::fit2pts(const char* fitTag)
{
  mesPlotsPath=combine("%s/plots/%s/%s",ens.dirPath.c_str(),mesComboTag.c_str(),fitTag);
  mkdir(mesPlotsPath);
  
  computeAxialPseudoCouplings();
  plotDispRel();
  
  return *this;
}

#endif
