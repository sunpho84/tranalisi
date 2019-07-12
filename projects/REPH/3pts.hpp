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
	if(ens.considerDec[iDecKin])
	  {
	    const vector<size_t> c=ens.indDecKin(iDecKin);
	    
	    for(int iVA=0;iVA<2;iVA++)
	      {
		corrPX[iVA][insOnT][iDecKin]=
		  load3pts(plotDirPath,iVA,iMs,iMt,c[0],c[1],c[2],"insOnT")*eT;
		corrPX[iVA][insOnS][iDecKin]=coeff[iVA]*
		  load3pts(plotDirPath,iVA,iMt,iMs,c[0],c[1],c[2],"insOnS")*eS;
	      }
	  }
      
      for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	if(ens.considerDec[iDecKin] and ens.hasSymmDec[iDecKin])
	  for(int iVA=0;iVA<2;iVA++)
	    for(int insOnST=0;insOnST<2;insOnST++)
	      {
		const djvec_t& a=corrPX[iVA][insOnST][iDecKin];
		const djvec_t& b=corrPX[iVA][insOnST][ens.symmOfDec[iDecKin]];
		const djvec_t& r=a/b;
		
		const vector<size_t> c=ens.indDecKin(iDecKin);
		
		const string pathR=combine("%s/o%smuGPo-gs_%s_insOn%s_rat_with_symm.xmg",plotDirPath.c_str(),VA_tag[iVA],ens.decKinTag(c[0],c[1],c[2]).c_str(),insOnST?"T":"S");
		r.ave_err().write(pathR);
	      }
	  }
  
  cout<<(loadCompact?"Loading":"Storing")<<" compacted 3pts"<<endl;
  
  raw_file_t file(dataPath,(loadCompact?"r":"w"));
  
  if(loadCompact)
    file.bin_read(corrPX);
  else
    file.bin_write(corrPX);
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::prepare3ptsNormalization(const bool useAnalytic,const bool& timeDependentEnergy)
{
  //! Main path
  const string plotDirPath=ens.dirPath+"/plots/"+mesTag+"/3pts_normalization/";
  if(not dir_exists(plotDirPath))
    mkdir(plotDirPath);
  
  cout<<"Preparing 3pts normalization"<<endl;
  
  const djack_t& mass=E[0];
  
  for(size_t iDecKin=0;iDecKin<ens.nDecKin;iDecKin++)
    if(ens.considerDec[iDecKin])
      {
	const size_t iMesKin=ens.iMesKinOfDecKin[iDecKin];
	const djack_t& EnNum=E[iMesKin];
	const djack_t EnAna=latt_en_1D(mass,ens.pMes[iMesKin]);;
	const djack_t& Energy=useAnalytic?EnAna:EnNum;
	const double& Eg=ens.Eg[iDecKin];
	const double EgT=ens.EgT(iDecKin);
	
	dEdec[iDecKin]=Energy-Eg;
	PKdec[iDecKin]=Energy*Eg-ens.pMes[iMesKin]*ens.kHatDec[iDecKin];
	X[iDecKin]=2*PKdec[iDecKin]/sqr(Energy);
	
	for(int t=0;t<=(int)ens.T/2;t++)
	  {
	    const djack_t& A=timeDependentEnergy?
	      eEff[iMesKin][std::min(std::min(t,ens.T-t),ens.T/2-1)]:
	      Energy;
	    
	    //const djack_t& A=E[iMesKin];
	    //const djack_t &B=Eeff[iMes][std::min(std::min(t,T-t),T/2-1)];
	    const djack_t& W=A;
	    const djack_t& Z=useAnalytic?ZP[0]:ZP[iMesKin];
	    //const djack_t& Z=ZP[iMesKin];
	    //cout<<t<<" "<<smart_print(A.ave_err())<<" "<<smart_print(B.ave_err())<<endl;
	    normaliz[iDecKin][t]=4*Energy*EgT/(Z*exp(-t*W-(ens.T/2-t)*Eg));
	    // cout<<"  "<<ens.indDecKin.descr(iDecKin)<<"W: "<<W.ave_err()<<", Eg: "<<Eg<<", Zp: "<<Z.ave_err()<<", n["<<t<<"]: "<<normaliz[iDecKin][t].ave_err()<<endl;
	  }
	
	const string normPath=combine("%s/%s.xmg",plotDirPath.c_str(),ens.decKinTag(iDecKin).c_str());
	grace_file_t plot(normPath);
	plot.set_title(combine("Ind %zu",iDecKin));
	plot.write_vec_ave_err(normaliz[iDecKin].ave_err());
      }
  
  return *this;
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::choose3ptsTint(const string& mesPlotsPath,const char* fitTag,const bool forceRechoose)
{
  string tintPath=milledPath()+"/tint3pts_"+fitTag+".dat";
  
  if(file_exists(tintPath) and not forceRechoose)
    {
      cout<<"Loading stored 3pts intervals"<<endl;
      raw_file_t(tintPath,"r").bin_read(tint3pts);
    }
  else
    {
      cout<<"Choosing tint from scratch"<<endl;
      
      const string rangeDir=ens.dirPath+"/plots/"+mesTag+"/3pts_range_"+fitTag+"/";
      mkdir(rangeDir);
      
      for(size_t iVA=0;iVA<2;iVA++)
	for(size_t iST=0;iST<2;iST++)
	  {
	    grace_file_t fitTxt(rangeDir+"/range_"+VA_tag[iVA]+"_"+ST_tag[iST]+".txt");
	    
	    for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	      if(ens.considerDec[iDecKin])
		{
		  // grace_file_t _outrange(combine("%s/_c%s_%s.xmg",path.c_str(),VA_tag[iVA],ens.decKinTag(iDecKin).c_str()));
		  
		  //! Three points effective mass
		  const djvec_t& c=corrPX[iVA][iST][iDecKin];
		  djvec_t y=c;
		  // if(0)
		    y*=normaliz[iDecKin];
		  // else
		  //   for(size_t t=0;t<y.size();t++)
		  //     y[t]/=exp(-dEdec[iDecKin]*t);
		  y=forward_derivative(y);
		  
		  // const djvec_t eff=effective_mass(c,T/2,0);
		  
		  const string outRangePath=combine("%s/c%s_insOn%c_%s.xmg",rangeDir.c_str(),VA_tag[iVA],ST_tag[iST],ens.decKinTag(iDecKin).c_str());
		  cout<<outRangePath<<endl;
		  
		  CompatibilityRangeFinder<> comp(y,outRangePath);
		  
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
					 .setRangeToConsiderByFraction(1.0/5)
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
		  
		  const Range& fitRange=tint3pts[iVA][iST][iDecKin]=comp.getSelectionRange(0);
		  
		  cout<<iDecKin<<" "<<ens.indDecKin.descr(iDecKin)<<": range "<<fitRange<<endl;
		  
		  const size_t &b=fitRange.begin,&e=fitRange.end;
		  fitTxt.write_ave_err(iDecKin,{(b+e)/2.0,(e-b)/2.0});
		}
	    
	    cout<<"Storing chosen 3pts intervals"<<endl;
	    raw_file_t(tintPath,"w").bin_write(tint3pts);
	    
		  /////////////////////////////////////////////////////////////////
		  
		  // vector<double> a(y.size()),e(y.size());
		  // for(size_t i=0;i<y.size();i++)
		  //   {
		  // 	a[i]=y[i].ave();
		  // 	e[i]=y[i].err();
		  //   }
		  
		  // const size_t firstCompDef=T/12;
		  // const size_t lastCompDef=compW0.size()-T/12;
		  // size_t firstComp=firstCompDef;
		  // size_t lastComp=lastCompDef;
		  
		  // //! loop until a decent compatibility range is obtained
		  // double nSig=1.0;
		  // do
		  //   {
		  // 	do firstComp++;
		  // 	while(compW0[firstComp]>nSig and firstComp<compW0.size());
		  
		  // 	do lastComp--;
		  // 	while(compW0[lastComp]>nSig and firstComp>0);
		  
		  // 	nSig*=1.1;
		  //   }
		  // while(lastComp==lastCompDef or firstComp==firstCompDef);
		  // cout<<"nSig: "<<nSig<<endl;
		  
		  // const size_t D=std::min(5,(int)(lastComp-firstComp));
		  // int tMinBest=0,tMaxBest=0;
		  // double pMax=0;
		  
		  // for(size_t tMin=firstComp;tMin<lastComp-D;tMin++)
		  //   for(size_t tMax=tMin+D;tMax<lastComp;tMax++)
		  // 	{
		  // 	  const size_t d=tMax-tMin+1;
		  
		  // 	  double Ch2Unc=0;
		  
		  // 	  for(size_t i=tMin;i<=tMax;i++)
		  // 	    Ch2Unc+=sqr(compW0[i]);
		  
		  // 	  Matrix<double,Dynamic,Dynamic> c(d,d);
		  // 	  for(size_t i=0;i<d;i++)
		  // 	    for(size_t j=0;j<d;j++)
		  // 	      c(i,j)=cov(y[i+tMin],y[j+tMin]);
		  
		  // 	  Matrix<double,Dynamic,Dynamic> cInv=c.inverse();
		  // 	  double Ch2Corr=0;
		  // 	  for(size_t i=0;i<d;i++)
		  // 	    for(size_t j=0;j<d;j++)
		  // 	      Ch2Corr+=a[i+tMin]*a[j+tMin]*cInv(i,j);
		  
		  // 	  //Ch2Corr/=tMax-tMin+1;
		  
		  // 	  if(Ch2Corr<1000 and Ch2Corr>0)
		  // 	    {
		  // 	      const double p=ch2Distr(Ch2Corr,d);
		  // 	      cout<<tMin<<" "<<tMax<<" "<<Ch2Corr<<" "<<Ch2Unc<<" "<<p<<endl;
		  // 	      if(p>pMax)
		  // 		{
		  // 		  pMax=p;
		  // 		  tMinBest=tMin;
		  // 		  tMaxBest=tMax;
		  // 		}
		  // 	    }
		  // 	}
		  
		  // if(tMinBest==0 and tMaxBest==0)
		  //   {
		  // 	tMinBest=T/4-D/2;
		  // 	tMaxBest=T/4+D/2;
		  //   }
		  
		  // cout<<"Best: "<<ens.indDecKin.descr(iDecKin)<<" ["<<tMinBest<<":"<<tMaxBest<<"] "<<pMax<<" within ["<<firstComp<<":"<<lastComp<<"]"<<endl;
		  // outrange.write_line([](double){return 0;},tMinBest-0.5,tMaxBest+0.5);
		  
		  // tMax=tMaxBest;
		  // tMin=tMinBest;
		  
		  /////////////////////////////////////////////////////////////////
	  }
    }
  
  return *this;
}

template <typename TV>
permes_combo_t<TV>& permes_combo_t<TV>::choose3ptsTintCommon(const size_t& length)
{
  for(size_t iVA=0;iVA<2;iVA++)
    for(size_t iST=0;iST<2;iST++)
      {
	vector<size_t> count(ens.T/2+1);
	
	for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	  if(ens.considerDec[iDecKin])
	    {
	      const Range& r=tint3pts[iVA][iST][iDecKin];
	      for(size_t t=r.begin;t<=r.end;t++)
		count[t]++;
	    }
	
	size_t tMinBest=0;
	auto f=[&](const size_t tStart){return accumulate(count.begin()+tStart,count.begin()+tStart+length,0);};
	for(size_t t=1;t<count.size()-length;t++)
	  {
	    if(f(t)>f(tMinBest))
	      tMinBest=t;
	  }
	commonTint3pts[iVA][iST]={tMinBest,tMinBest+length};
	
	cout<<"Common fit range for "<<VA_tag[iVA]<<" "<<ST_tag[iST]<<": "<<commonTint3pts[iVA][iST]<<endl;
      }
  
  
  return *this;
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
      cout<<"Loading stored form factors"<<endl;
      raw_file_t file(ffPath,"r");
      file.bin_read(ff);
      file.bin_read(this->quality);
    }
  else
    {
      cout<<"Fitting 3pts correlators"<<endl;
      
      for(int iVA=0;iVA<2;iVA++)
	for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	  if(ens.considerDec[iDecKin])
	    {
	      const int iMesKin=ens.iMesKinOfDecKin[iDecKin];
	      
	      djack_t H=0.0;
	      djvec_t yTot(ens.T/2+1);
	      yTot=0.0;
	      
	      this->quality[iVA][iDecKin]=true;
	      
	      for(size_t iST=0;iST<2;iST++)
		{
		  const Range& tint=
		    useCommonRange?
		    commonTint3pts[iVA][iST]:
		    tint3pts[iVA][iST][iDecKin];
		  const djvec_t y=corrPX[iVA][iST][iDecKin]*normaliz[iDecKin];
		  yTot+=y;
		  
		  // corrPX[iVA][iST][iDecKin].ave_err().write(combine("%s/c%s_insOn%c_%s_f1.xmg",mesPlotsPath.c_str(),VA_tag[iVA],ST_tag[iST],ens.decKinTag(iDecKin).c_str()));
		  // normaliz[iDecKin].ave_err().write(combine("%s/c%s_insOn%c_%s_f2.xmg",mesPlotsPath.c_str(),VA_tag[iVA],ST_tag[iST],ens.decKinTag(iDecKin).c_str()));
		  
		  const string plotFile=combine("%s/c%s_insOn%c_%s.xmg",mesPlotsPath.c_str(),VA_tag[iVA],ST_tag[iST],ens.decKinTag(iDecKin).c_str());
		  const djack_t Hcontr=constant_fit(y,tint.begin,tint.end,plotFile);
		  
		  H+=Hcontr;
		  
		  djack_t ch2=0;
		  for(size_t t=tint.begin;t<=tint.end;t++)
		    {
		      const djack_t contr=(y[t]-Hcontr)/y[t].err();
		      ch2+=sqr(contr);
		    }
		  
		  const size_t nDof=tint.size()-1;
		  cout<<"Chi2 "<<plotFile<<": "<<smart_print(ch2)<<" / "<<nDof<<endl;
		  const double q=ch2.ave()/nDof;
		  
		  if(q>2.0) this->quality[iVA][iDecKin]=false;
		}
	      
	      yTot.ave_err().write(combine("%s/c%s_%s.xmg",mesPlotsPath.c_str(),VA_tag[iVA],ens.decKinTag(iDecKin).c_str()));
	      
	      const double eps=1/sqrt(2.0);
	      
	      if(iVA==1)
		{
		  const djack_t& M=E[0];
		  const djack_t& PK=PKdec[iDecKin];
		  const djack_t div=fPbare[0]*(eT-eS);
		  
		  ff[iVA][iDecKin]=M*(H/eps-div)/PK;
		}
	      else
		ff[iVA][iDecKin]=H/eps*E[iMesKin]/(ens.Eg[iDecKin]*ens.pMes[iMesKin]-E[iMesKin]*ens.kHatDec[iDecKin]);
	    }
      
      cout<<"Storing form factors"<<endl;
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
      
      ffPlot.set_line_type(grace::line_type_t::STRAIGHT_LINE);
      ffPlot.set_color_scheme({RED,BLUE,GREEN4,VIOLET,ORANGE});

      for(size_t toC=0;toC<2;toC++)
	for(size_t iMom1=0;iMom1<ens.nMoms;iMom1++)
	  for(size_t iMom2=0;iMom2<ens.nMoms;iMom2++)
	    {
	      if(toC==0)
		ffPlot.set_transparency(0.4);
	      
	      for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
		{
		  const vector<size_t> c=ens.indDecKin(iDecKin);
		  
		  if(ens.considerDec[iDecKin] and c[1]==iMom1 and c[2]==iMom2 and toC==quality[iVA][iDecKin])
		    {
		      ffPlot<<"# "<<iDecKin<<" "<<ens.indDecKin.descr(iDecKin)<<"\n";
		      ffPlot.write_ave_err(X[iDecKin].ave(),ff[iVA][iDecKin].ave_err());
		    }
		}
	      ffPlot.new_data_set();
	    }
    }
}

#endif
