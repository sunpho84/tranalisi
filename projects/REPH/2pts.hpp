#ifndef _2PTS_HPP
#define _2PTS_HPP

#include <REPH/permes_combo.hpp>

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

void permes_combo_t::plotDispRel() const
{
  grace_file_t dispRel(combine("%s/dispRel.xmg",mesPlotsPath.c_str()));
  dispRel.write_vec_ave_err(ens.pMes,E.ave_err());
  
  const djack_t &M=E[0];
  
  dispRel.write_polygon([M](double x){return latt_en_1D(M,x);},0,ens.pMesMax,grace::GREEN);
  dispRel.write_polygon([M](double x){return cont_en_1D(M,x);},0,ens.pMesMax,grace::VIOLET);
}

void permes_combo_t::computeAxialPseudoCouplings()
{
  const size_t& T=ens.T;
  
  const size_t tMin=12,tMax=T/2+1;
  
  grace_file_t ZPPlot(combine("%s/ZP.xmg",mesPlotsPath.c_str()));
  grace_file_t ZAPlot(combine("%s/ZA.xmg",mesPlotsPath.c_str()));

  for(auto dir : {"2pts_PP_corr","2pts_AP0_corr","2pts_AP3_corr","CPP","CA0P","CA3P","CPP_A0P"})
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
      
      const std::string kinTag=combine("im1_%zu_im2_%zu",iMom1,iMom2);
      
      const djvec_t corrPP=load2ptsPP(iMom1,iMom2);
      corrPP.ave_err().write(combine("%s/2pts_PP_corr/%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()));
      eEff[iMesKin]=effective_mass(corrPP);
      const djvec_t corrA0P=load2ptsAP(iMom1,iMom2,0);
      corrA0P.ave_err().write(combine("%s/2pts_AP0_corr/%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()));
      
      const djvec_t corrA3P=load2ptsAP(iMom1,iMom2,3);
      corrA3P.ave_err().write(combine("%s/2pts_AP3_corr/%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()));
      
      if(iMom1!=iMom2)
	{
	  const double P=ens.pMes[iMesKin];
	  
	  //perform a preliminary fit
	  djack_t fit_E,ZPP,ZA0P,ZA3P=1.0;
	  two_pts_fit(ZPP, fit_E,corrPP, T/2,tMin,tMax,combine("%s/CPP/%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),"",+1);
	  two_pts_fit(ZA0P,fit_E,corrA0P,T/2,tMin,tMax-1,combine("%s/CA0P/%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),"",-1);
	  two_pts_fit(ZA3P,fit_E,corrA3P,T/2,tMin,tMax,combine("%s/CA3P/%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),"",+1);
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
	  auto x=vector_up_to<double>(corrPP.size());
	  multi_ch2_t<djvec_t> two_pts_fit_obj({x,x,x},{tMin,tMin,tMin},{tMax,tMax-1,tMax},{corrPP,corrA0P,corrA3P},
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
	  
	  write_constant_fit_plot(combine("%s/CPP/fit_%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),tMin,tMax,E[iMesKin],effective_mass(corrPP,T/2,+1));
	  write_constant_fit_plot(combine("%s/CA0P/fit_%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),tMin,tMax,E[iMesKin],effective_mass(corrA0P,T/2,-1));
	  write_constant_fit_plot(combine("%s/CA3P/fit_%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),tMin,tMax,E[iMesKin],effective_mass(corrA3P,T/2,+1));
	}
      else
	{
	  two_pts_SL_fit(ZP[iMesKin],ZA[iMesKin],E[iMesKin],corrA0P,corrPP,T/2,tMin,tMax-1,combine("%s/CPP_A0P/fit_%s.xmg",mesPlotsPath.c_str(),kinTag.c_str()),-1,+1);
	  ZA[iMesKin]/=E[iMesKin];
	}
      
      ZAPlot.write_ave_err(E[iMesKin].ave(),ZA[iMesKin].ave_err());
      ZPPlot.write_ave_err(E[iMesKin].ave(),ZP[iMesKin].ave_err());
      
      fP[iMesKin]=ZP[iMesKin]*(ens.mass[iMs]+ens.mass[iMt])/sqr(E[iMesKin]);
      fPbare[iMesKin]=ZA[iMesKin];
    }
}

#endif
