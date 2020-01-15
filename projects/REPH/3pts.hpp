#ifndef _3PTS_HPP
#define _3PTS_HPP

#include <REPH/permes_combo.hpp>

template <typename TV>
djvec_t permes_combo_t<TV>::load3pts(const string& plotDirPath,const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0,const char* extraTag)
{
  const size_t& T=ens.T;
  const size_t& L=ens.L;
  const size_t& nMass=ens.nMass;
  const size_t& nMoms=ens.nMoms;
  
  //! Sign of the contribution
  const double s[2][2][2]={{{-1,+1},{-1,-1}},
			   {{-1,-1},{+1,-1}}};
  
  //! Output correlator
  djvec_t corr(T);
  corr=0.0;
  
  const int iDecKin=ens.indDecKin({iMoms,iMomt,iMom0});
  const double& Eg=ens.Eg[iDecKin];
  
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
  
  //! Real part for Axial, Imaginary for vector
  const size_t iReIm=(iVA==0)?1:0;
  
  //! Negative parity for axial, positive for vector
  const int par=(iVA==0)?-1:+1;
  
  for(size_t iPol=0;iPol<2;iPol++)
    for(size_t iGamma=1;iGamma<=2;iGamma++)
      {
	const string path=combine("%s/jacks/o%smuGPo-gs",ens.dirPath.c_str(),VA_tag[iVA]);
	const size_t i=ind({iMs,iMt,iMoms,iMomt,iMom0,iPol,iGamma,iReIm});
	
	const djvec_t contr=-ens.loadCorr(path,i,"3")/((double)L*L*L);
	
	corr+=contr*s[iVA][iPol][iGamma-1];
	
	grace_file_t plot(combine("%s/o%smuGPo-gs_%s_pol_%zu_gamma_%zu_%s.xmg",plotDirPath.c_str(),VA_tag[iVA],ens.decKinTag(iMoms,iMomt,iMom0).c_str(),iPol,iGamma,extraTag));
	plot.set_title(combine("Ind %zu",i));
	plot.write_vec_ave_err(contr.ave_err());
      }
  corr/=4;
  
  corr.symmetrize(par);
  const double EgT=ens.EgT(iDecKin);
  for(int t=0;t<=ens.T/2;t++)
    {
      corr[t]/=exp(-(ens.T/2-t)*Eg);
      
      if(oldNormalization==1 or Eg>1e-5)
	corr[t]*=2*EgT;
    }
  
  corr.ave_err().write(combine("%s/o%smuGPo-gs_%s_%s.xmg",plotDirPath.c_str(),VA_tag[iVA],ens.decKinTag(iMoms,iMomt,iMom0).c_str(),extraTag));
  
  return corr;
}

template <typename TV>
void permes_combo_t<TV>::load3pts(const bool forceLoad)
{
  cout<<"Loading 3pts correlators"<<endl;
  
  //! Coefficient to combine the two insertions
  const double coeff[2]={+1.0,-1.0};
  
  //! Path where to store data
  const string dataPath=milledPath()+"/data3pts.dat";
  
  //! Load or not the compacted data
  const bool loadCompact=file_exists(dataPath) and not forceLoad;
  
  if(not loadCompact)
    {
      //! Main path
      const string plotDirPath=ens.dirPath+"/plots/"+mesTag+"/3pts_corr/";
      if(not dir_exists(plotDirPath))
	mkdir(plotDirPath);
      
      cout<<"Loading 3pts correlators from scratch"<<endl;
      
      for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	{
	  //! Components of the decay
	  const vector<size_t> c=ens.indDecKin(iDecKin);
	  
	  // Load fixing signs and charges
	  for(int iVA=0;iVA<2;iVA++)
	    {
	      corrPX[iVA][iDecKin]=
		load3pts(plotDirPath,iVA,iMs,iMt,c[0],c[1],c[2],"insOnT")*eT;
	      corrPX[iVA][iDecKin]+=coeff[iVA]*
		load3pts(plotDirPath,iVA,iMt,iMs,c[0],c[1],c[2],"insOnS")*eS;
	    }
	}
      
      // Plot ratio with symmetric
      // for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
      // 	if(ens.considerDec[iDecKin] and ens.hasSymmDec[iDecKin])
      // 	  for(int iVA=0;iVA<2;iVA++)
      // 	    {
      // 	      const djvec_t& a=corrPX[iVA][iDecKin];
      // 	      const djvec_t& b=corrPX[iVA][ens.symmOfDec[iDecKin]];
      // 	      const djvec_t& r=a/b;
	      
      // 	      const vector<size_t> c=ens.indDecKin(iDecKin);
	      
      // 	      const string pathR=combine("%s/o%smuGPo-gs_%s_rat_with_symm.xmg",plotDirPath.c_str(),VA_tag[iVA],ens.decKinTag(c[0],c[1],c[2]).c_str());
      // 	      r.ave_err().write(pathR);
      // 	    }
    }
  
  cout<<(loadCompact?"Loading":"Storing")<<" compacted 3pts"<<endl;
  
  //! File where to load or store
  raw_file_t file(dataPath,(loadCompact?"r":"w"));
  
  if(loadCompact)
    file.bin_read(corrPX);
  else
    file.bin_write(corrPX);
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::prepareKinematics(const bool useAnalytic)
{
  cout<<"Preparing 3pts kinematics"<<endl;
  
  const djack_t& mass=E[0];
  
  for(size_t iDecKin=0;iDecKin<ens.nDecKin;iDecKin++)
    if(ens.considerDec[iDecKin])
      {
	const size_t iMesKin=ens.iMesKinOfDecKin[iDecKin];
	const djack_t& EnNum=E[iMesKin];
	const djack_t EnAna=latt_en_1D(mass,ens.pMes[iMesKin]);;
	const djack_t& Energy=useAnalytic?EnAna:EnNum;
	const double& Eg=ens.Eg[iDecKin];
	
	dEdec[iDecKin]=Energy-Eg;
	PKdec[iDecKin]=Energy*Eg-ens.pMes[iMesKin]*ens.kHatDec[iDecKin];
	X[iDecKin]=2*PKdec[iDecKin]/sqr(Energy);
      }
  
  return *this;
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::choose3ptsTint(const string& mesPlotsPath,const char* fitTag,const bool forceRechoose)
{
  //! Path of the 3pts time interval file
  string tintPath=milledPath()+"/tint3pts_"+fitTag+".dat";
  
  //! Directory path of the range search
  const string rangeDir=ens.dirPath+"/plots/"+mesTag+"/3pts_range_"+fitTag+"/";
  mkdir(rangeDir);
  
  if(file_exists(tintPath) and not forceRechoose)
    {
      cout<<"Loading stored 3pts intervals"<<endl;
      raw_file_t(tintPath,"r").bin_read(tint3pts);
    }
  else
    {
      cout<<"Choosing tint from scratch"<<endl;
      
      for(size_t iVA=0;iVA<2;iVA++)
	{
	  //! Path of the fit range file
	  grace_file_t fitTxt(rangeDir+"/range_"+VA_tag[iVA]+".txt");
	  
	  for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	    {
	      //! Three points ratio derivative
	      const djvec_t r=forward_derivative(getCorrRat(iVA,iDecKin));
	      
	      Range fitRange;
	      if(r[0].err()==0)
		fitRange={(size_t)((int)ens.T/4-1),(size_t)((int)(ens.T/4)+1)};
	      else
		{
		  //! Path where to store the range search
		  const string outRangePath=combine("%s/c%s_%s.xmg",rangeDir.c_str(),VA_tag[iVA],ens.decKinTag(iDecKin).c_str());
		  cout<<outRangePath<<endl;
		  
		  //! Range finder
		  CompatibilityRangeFinder<> comp(r,outRangePath);
		  
		  //! Initial value of enlarging factor
		  double nSigEnl=2.0;
		  double enlFact=1.8;
		  
		  //! Minimal length
		  const size_t lengthMin=4;
		  
		  const auto recipe=[T=ens.T,&comp,&enlFact,&nSigEnl](const bool verb=false)
				    {
				      double nSigSel=1;
				      comp
					.setVerbose(verb)
					.setRangeToConsiderByFraction(1.0/8)
					.selectClosestCompatiblePointWithinNsigma(T/4,nSigSel)
					.plotSelected()
					.selectEnlargingError(enlFact,nSigEnl)
					.plotSelected()
					.mergeSelectionByDistanceSize()
					.plotSelected()
					.selectLargestRange()
					.plotSelected();
				      
				      //! Trim some points if allowed
				      const size_t nDesel=2;
				      if(comp.getSelectionRange(0).size()>=lengthMin+nDesel)
					comp
					  .deSelectLeftermostPoints(nDesel)
					  .deSelectRightermostPoints(nDesel)
					  .plotSelected();
				    };
		  
		  recipe();
		  
		  // Extend
		  while(comp.getSelectionRange(0).size()<lengthMin)
		    {
		      nSigEnl+=1.0;
		      enlFact*=1.2;
		      
		      cout<<"WARNING, compatibility range too short, "<<comp.getSelectionRange(0)<<" when at least "<<lengthMin<<" needed, enlarging of "<<nSigEnl<<endl;
		      
		      recipe(true);
		    }
		  
		  //! Store the obtained fit range
		  fitRange=comp.getSelectionRange(0);
		}
	      
	      tint3pts[iVA][iDecKin]=fitRange;
	      
	      cout<<iDecKin<<" "<<ens.indDecKin.descr(iDecKin)<<": range "<<fitRange<<endl;
	      
	      const size_t &b=fitRange.begin,&e=fitRange.end;
	      fitTxt.write_ave_err(iDecKin,{(b+e)/2.0,(e-b)/2.0});
	    }
	  
	  cout<<"Storing chosen 3pts intervals"<<endl;
	  raw_file_t(tintPath,"w").bin_write(tint3pts);
	}
    }
  
  const int tmin=2;
  const int tmax=ens.T/2-2;
  
  for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
    {
      const vector<size_t> c=ens.indDecKin(iDecKin);
      if(c[1]!=c[2])
	{
	  const djvec_t r=getCorrRat(1,iDecKin)/X[iDecKin];
	  
	  const size_t nFitPars=5;
	  djvec_t pFit(nFitPars);
	  jack_fit_t fit;
	  fit.add_fit_par(pFit[0],"p[0]",1.26e-1,0.1);
	  fit.add_fit_par_limits(pFit[1],"p[1]",7e-1,0.1, 0.15,4);
	  fit.add_fit_par(pFit[2],"p[2]",-9e-1,0.1);
	  fit.add_fit_par_limits(pFit[3],"p[3]",3e-1,0.1, 0.15,4);
	  fit.add_fit_par(pFit[4],"p[4]",r[ens.T/4].ave(),0.1);
	  
	  // for(size_t i=0;i<nFitPars;i++)
	  //   {
	      
	  //     if(i==1 or i==3) 
	  //     else             fit.add_fit_par(pFit[i],combine("p[%zu]",i),0,0.1);
	      
	  //     //if(i==0 or i==1) fit.fix_par(i);
	  //   }
	  
	  for(int t=tmin;t<tmax;t++)
	    fit.add_point(r[t],[=](const vector<double>& p,int ijack)->double
			       {
				 return
				   p[0]*exp(-p[1]*t)+
				   p[2]*exp(-p[3]*(ens.T/2-t))+
				   p[4];
			       });
	  
	  fit.fit();
	  
	  const string outRangePath=combine("%s/d%s_%s.xmg",rangeDir.c_str(),VA_tag[1],ens.decKinTag(iDecKin).c_str());
	  grace_file_t plot(outRangePath);
	  
	  plot.write_vec_ave_err(r.ave_err());
	  
	  plot.write_polygon([=](const double x) -> djack_t
			     {
			       return
				 pFit[0]*exp(-pFit[1]*x)+
				 pFit[2]*exp(-pFit[3]*(ens.T/2-x))+
				 pFit[4];
			     },tmin,tmax);
	  
	  plot.write_constant_band(tmin,tmax,pFit[4]);
	  
	  {
	    const string outRangePath=combine("%s/e%s_%s.xmg",rangeDir.c_str(),VA_tag[1],ens.decKinTag(iDecKin).c_str());
	    
	    grace_file_t plot(outRangePath);
	    plot.write_vec_ave_err(r.ave_err());
	    
	    plot.write_polygon([=](const double x) -> djack_t
			     {
			       return
				 1/(pFit[0]*exp(-pFit[1]*x)+
				    pFit[2]*exp(-pFit[3]*(ens.T/2-x)));
			     },tmin,tmax);
	  
	    plot.write_polygon([=](const double x) -> djack_t
			     {
			       return
				 1/(abs(pFit[0]*exp(-pFit[1]*x))+
				    abs(pFit[2]*exp(-pFit[3]*(ens.T/2-x))));
			     },tmin,tmax);
	    
	  }
	  
	  for(size_t i=0;i<5;i++) cout<<"par"<<i<<" "<<X[iDecKin]<<" "<<pFit[i].ave_err()<<endl;
	}
    }
  
  CRASH("");
  
  return *this;
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::choose3ptsTintCommon(const size_t& length)
{
  for(size_t iVA=0;iVA<2;iVA++)
    {
      vector<size_t> count(ens.T/2+1);
      
      for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	{
	  const Range& r=tint3pts[iVA][iDecKin];
	  for(size_t t=r.begin;t<r.end;t++)
	    count[t]++;
	}
      
      size_t tMinBest=0;
      auto f=[&](const size_t tStart){return accumulate(count.begin()+tStart,count.begin()+tStart+length,0);};
      for(size_t t=1;t<count.size()-length;t++)
	{
	  if(f(t)>f(tMinBest))
	    tMinBest=t;
	}
      commonTint3pts[iVA]={tMinBest,tMinBest+length};
      
      cout<<"Common fit range for "<<VA_tag[iVA]<<": "<<commonTint3pts[iVA]<<endl;
    }
  
  return *this;
}

template <typename TV>
TV permes_combo_t<TV>::getCorrRat(const int iVA,const size_t iDecKin)
{
  const vector<size_t> c=ens.indDecKin(iDecKin);
  vector<size_t> c0=c;
  c0[1]=c0[2];
  const auto iDecSmallestKin=ens.indDecKin(c0);
  
  TV r=(corrPX[iVA][iDecKin]-corrPX[iVA][iDecSmallestKin])/corrPX[1][iDecSmallestKin];
  
  return r;
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::fit3pts(const bool& useCommonRange,const char* fitTag,const bool forceRechoose)
{
  const std::string mesPlotsPath=ens.dirPath+"/plots/"+mesTag+"/3pts_"+fitTag;
  if(not dir_exists(mesPlotsPath))
    mkdir(mesPlotsPath);
  
  choose3ptsTint(mesPlotsPath,fitTag,forceRechoose);
  choose3ptsTintCommon();
  
  string ffPath=milledPath()+"/ff_"+fitTag+".dat";
  if(file_exists(ffPath) and not forceRechoose)
    {
      cout<<"Loading stored form factors and mel"<<endl;
      raw_file_t file(ffPath,"r");
      file.bin_read(ff);
      file.bin_read(this->quality);
    }
  else
    {
      cout<<"Fitting 3pts correlators"<<endl;
      
      for(int iVA=0;iVA<2;iVA++)
	for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	  {
	    const size_t iMesKin=ens.iMesKinOfDecKin[iDecKin];
	    
	    if(oldNormalization)
	      CRASH("NOT POSSIBLE ANY MORE");
	    
	    this->quality[iVA][iDecKin]=true;
	    
	    const Range& tint=
	      useCommonRange?
	      commonTint3pts[iVA]:
	      tint3pts[iVA][iDecKin];
	    
	    const string plotFile=combine("%s/melRat_c%s_%s.xmg",mesPlotsPath.c_str(),VA_tag[iVA],ens.decKinTag(iDecKin).c_str());
	    
	    const TV r=getCorrRat(iVA,iDecKin);
	    
	    const djack_t m=constant_fit(r,tint.begin,tint.end,plotFile);
	    
	    if(iVA==0)
	      ff[iVA][iDecKin]=m*E[iMesKin]/(ens.Eg[iDecKin]*ens.pMes[iMesKin]-E[iMesKin]*ens.kHatDec[iDecKin]);
	    else
	      ff[iVA][iDecKin]=m/X[iDecKin];
	    
	    // //! Put back fpi
	    // const djack_t fpi=fP[0]*(eT-eS);
	    // ff[iVA][iDecKin]*=fpi;
	  }
      
      cout<<"Storing matrix element and ff"<<endl;
      raw_file_t file(ffPath,"w");
      file.bin_write(ff);
      file.bin_write(this->quality);
    }
  
  return *this;
}

template <typename TV>
void permes_t<TV>::plotFf(const string& tag) const
{
  const string path=ens.dirPath+"/plots/"+mesTag;
  mkdir(path);
  
  cout<<"Plotting ff"<<endl;
  
  for(int iVA=0;iVA<2;iVA++)
    {
      using namespace grace;
      grace_file_t ffPlot(combine("%s/ff_%s" "%s" ".xmg",path.c_str(),VA_tag[iVA],((tag=="")?tag:("_"+tag)).c_str()));
      
      vector<decltype(make_tuple(&ffPlot,&ff))> list={make_tuple(&ffPlot,&ff)};
      
      for(auto& fy : list)
	{
	  grace_file_t& f=*get<0>(fy);
	  auto& y=*get<1>(fy);
	  
	  f.set_line_type(grace::line_type_t::STRAIGHT_LINE);
	  f.set_color_scheme({RED,BLUE,GREEN4,VIOLET,ORANGE});
	  
	  //for(size_t toC=0;toC<2;toC++)
	    for(size_t iMom1=0;iMom1<ens.nMoms;iMom1++)
	      for(size_t iMom2=0;iMom2<ens.nMoms;iMom2++)
		{
		  // if(toC==0)
		  //   f.set_transparency(0.4);
		  
		  for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
		    {
		      const vector<size_t> c=ens.indDecKin(iDecKin);
		      
		      if(ens.considerDec[iDecKin] and c[0]==iMom1 and c[2]==iMom2 // and toC==quality[iVA][iDecKin]
			 )
			{
			  f<<"# "<<iDecKin<<" "<<ens.indDecKin.descr(iDecKin)<<"\n";
			  f.write_ave_err(X[iDecKin].ave(),y[iVA][iDecKin].ave_err());
			}
		    }
		  f.set_comment(combine("iMoms=%zu, iMom0=%zu",iMom1,iMom2));
		  f.new_data_set();
		}
	}
    }
}

#endif
