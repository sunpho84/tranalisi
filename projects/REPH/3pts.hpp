#ifndef _3PTS_HPP
#define _3PTS_HPP

#include <REPH/permes_combo.hpp>

djvec_t permes_combo_t::load3pts(const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0)
{
  const string path=ens.dirPath+"/plots/"+mesComboTag+"/3pts_corr/";
  mkdir(path);
  
  const std::string milledPath=ens.dirPath+"/milledData/";
  mkdir(milledPath);
  
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
	
	contr.ave_err().write(combine("%s/o%smuGPo-gs_%s_pol_%zu_gamma_%zu.xmg",path.c_str(),VA_tag[iVA],ens.decKinTag(iMoms,iMomt,iMom0).c_str(),iPol,iGamma));
      }
  corr/=2*sqrt(2);
  
  corr/=-((double)L*L*L);
  
  corr.symmetrize(par);
  corr.ave_err().write(combine("%s/o%smuGPo-gs_%s.xmg",path.c_str(),VA_tag[iVA],ens.decKinTag(iMoms,iMomt,iMom0).c_str()));
  
  corr.bin_write(combine("%s/c%s_mS%zu_mT%zu_momS%zu_momT%zu_momO%zu.dat",milledPath.c_str(),VA_tag[iVA],iMs,iMt,iMoms,iMomt,iMom0));
  
  return corr;
}

void permes_combo_t::load3pts()
{
  cout<<"Loading 3pts correlators"<<endl;
  
  //! Coefficient to combine the two insertions
  const double coeff[2]={+1.0,-1.0};
  
  for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
    if(ens.considerDec[iDecKin])
      {
	const vector<size_t> c=ens.indDecKin(iDecKin);
	
	for(int iVA=0;iVA<2;iVA++)
	  {
	    corrPX[iVA][insOnT][iDecKin]=
	      load3pts(iVA,iMs,iMt,c[0],c[1],c[2])*eT;
	    corrPX[iVA][insOnS][iDecKin]=coeff[iVA]*
	      load3pts(iVA,iMt,iMs,c[0],c[1],c[2])*eS;
	    // cout<<eT<<" "<<eS<<" "<<corrPX[iVA][iDecKin][10].ave_err()<<endl;
	  }
      }
}

permes_combo_t& permes_combo_t::prepare3ptsNormalization()
{
  cout<<"Preparing 3pts normalization"<<endl;
  
  for(size_t iDecKin=0;iDecKin<ens.nDecKin;iDecKin++)
    {
      const size_t iMesKin=ens.iMesKinOfDecKin[iDecKin];
      const double& Eg=ens.Eg[iDecKin];
      const double EgT=ens.EgT(iDecKin);
      
      dEdec[iDecKin]=E[iMesKin]-Eg;
      PKdec[iDecKin]=E[iMesKin]*Eg-ens.pMes[iMesKin]*ens.kHatDec[iDecKin];
      X[iDecKin]=2*PKdec[iDecKin]/sqr(E[iMesKin]);
      
      for(int t=0;t<=(int)ens.T/2;t++)
	{
	  const djack_t &A=E[iMesKin];
	  //const djack_t &B=Eeff[iMes][std::min(std::min(t,T-t),T/2-1)];
	  const djack_t &W=A;
	  //cout<<t<<" "<<smart_print(A.ave_err())<<" "<<smart_print(B.ave_err())<<endl;
	  normaliz[iDecKin][t]=4*W*EgT/(ZP[iMesKin]*exp(-t*W-(ens.T/2-t)*Eg));
	}
    }
  
  return *this;
}

permes_combo_t& permes_combo_t::chooseTint()
{
  const string path=ens.dirPath+"/plots/"+mesComboTag+"/3pts_range/";
  mkdir(path);
  
  cout<<"Choosing the 3pts fit interval"<<endl;
  
  //! Proxy T
  const size_t& T=ens.T;
  
  for(size_t iVA=0;iVA<2;iVA++)
    for(size_t iST=0;iST<2;iST++)
      {
	for(size_t i3ptsKin=0;i3ptsKin<ens.indDecKin.max();i3ptsKin++)
	  if(ens.considerDec[i3ptsKin])
	    {
	      // grace_file_t _outrange(combine("%s/_c%s_%s.xmg",path.c_str(),VA_tag[iVA],ens.decKinTag(i3ptsKin).c_str()));
	      
	      //! Three points effective mass
	      const djvec_t& c=corrPX[iVA][iST][i3ptsKin];
	      const djvec_t eff=effective_mass(c,T/2,0);
	      
	      //! Three points effective mass subtracted of the expected behaviour
	      const djvec_t y=eff-dEdec[i3ptsKin];
	      
	      const string outRangePath=combine("%s/c%s_insOn%c_%s.xmg",path.c_str(),VA_tag[iVA],ST_tag[iST],ens.decKinTag(i3ptsKin).c_str());
	      cout<<outRangePath<<endl;
	      
	      CompatibilityRangeFinder<> comp(y,outRangePath);
	      
	      //! Initial value of enlarging factor
	      double nSigEnl=2.0;
	      
	      const auto receipt=[T=ens.T,&comp,&nSigEnl](const bool verb=false)
				 {
				   double nSigSel=1;
				   comp
				     .setVerbose(verb)
				     .setRangeToConsiderByFraction(1.0/12)
				     .selectClosestCompatiblePointWithinNsigma(T/4,nSigSel)
				     .plotSelected()
				     .selectEnlargingError(2,nSigEnl)
				     .plotSelected()
				     .mergeSelectionByDistanceSize()
				     .plotSelected()
				     .selectLargestRange()
				     .plotSelected();
				 };
	      
	      receipt();
	      
	      //! Minimal length
	      const size_t lengthMin=3;
	      
	      // Extend
	      while(comp.getSelectionRange(0).size()<lengthMin)
		{
		  nSigEnl+=1.0;
		  
		  cout<<"WARNING, compatibility range too short, "<<comp.getSelectionRange(0)<<" when at least "<<lengthMin<<" needed, enlarging of "<<nSigEnl<<endl;
		  
		  receipt(true);
		}
	      
	      const Range& fitRange=tint3pts[iVA][iST][i3ptsKin]=comp.getSelectionRange(0);
	      
	      cout<<i3ptsKin<<" "<<ens.indDecKin.descr(i3ptsKin)<<": range "<<fitRange<<endl;
	      
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
	      
	      // cout<<"Best: "<<ens.indDecKin.descr(i3ptsKin)<<" ["<<tMinBest<<":"<<tMaxBest<<"] "<<pMax<<" within ["<<firstComp<<":"<<lastComp<<"]"<<endl;
	      // outrange.write_line([](double){return 0;},tMinBest-0.5,tMaxBest+0.5);
	      
	      // tMax=tMaxBest;
	      // tMin=tMinBest;
	      
	      /////////////////////////////////////////////////////////////////
	    }
      }
  
  return *this;
}

permes_combo_t& permes_combo_t::fit3pts(const char* fitTag)
{
  cout<<"Fitting 3pts correlators"<<endl;
  
  const std::string path=ens.dirPath+"/plots/"+mesComboTag+"/3pts_fit";
  if(not dir_exists(path))
    mkdir(path);
  
  for(int iVA=0;iVA<2;iVA++)
    for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
      if(ens.considerDec[iDecKin] // and hasSymm[iDecKin]
	 )
	{
	  const int iMesKin=ens.iMesKinOfDecKin[iDecKin];
	  
	  djack_t H=0.0;
	  for(size_t iST=0;iST<2;iST++)
	    {
	      const Range& tint=tint3pts[iVA][iST][iDecKin];
	      const djvec_t y=corrPX[iVA][iST][iDecKin]*normaliz[iDecKin];
	      H+=constant_fit(y,tint.begin,tint.end
			      ,combine("%s/c%s_insOn%c_%s.xmg",path.c_str(),VA_tag[iVA],ST_tag[iST],ens.decKinTag(iDecKin).c_str()));
	    }
	  
	  if(iVA==1)
	    ff[iVA][iDecKin]=(H-fPbare[iMesKin]*(eT-eS))/PKdec[iDecKin];
	  else
	    ff[iVA][iDecKin]=H*E[iMesKin]/(ens.Eg[iDecKin]*ens.pMes[iMesKin]-E[iMesKin]*ens.kHatDec[iDecKin]);
	}
  
  return *this;
}

permes_combo_t& permes_combo_t::plotFf()
{
  const string path=ens.dirPath+"/plots/"+mesComboTag;
  mkdir(path);
  
  cout<<"Plotting ff"<<endl;
  
  for(int iVA=0;iVA<2;iVA++)
    {
      using namespace grace;
      grace_file_t ffPlot(combine("%s/ff_%s.xmg",path.c_str(),VA_tag[iVA]));
      ffPlot.set_line_type(grace::line_type_t::STRAIGHT_LINE);
      ffPlot.set_color_scheme({RED,BLUE,GREEN4,VIOLET,ORANGE});
      
      for(size_t iMom1=0;iMom1<ens.nMoms;iMom1++)
	for(size_t iMom2=0;iMom2<ens.nMoms;iMom2++)
	  {
	    for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	      {
		const vector<size_t> c=ens.indDecKin(iDecKin);
		
		if(ens.considerDec[iDecKin] and c[1]==iMom1 && c[2]==iMom2)
		  {
		    ffPlot<<"# "<<iDecKin<<" "<<ens.indDecKin.descr(iDecKin)<<"\n";
		    ffPlot.write_ave_err(X[iDecKin].ave(),ff[iVA][iDecKin].ave_err());
		  }
	      }
	    ffPlot.new_data_set();
	  }
    }
  
  return *this;
}

#endif
