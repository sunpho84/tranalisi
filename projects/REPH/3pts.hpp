#ifndef _3PTS_HPP
#define _3PTS_HPP

#include <REPH/permes_combo.hpp>

std::string decKinTag(const size_t iMoms,const size_t iMomt,const size_t iMom0)
{
  return combine("iMomS%zu_iMomT%zu_iMomO%zu",iMoms,iMomt,iMom0);
}

djvec_t permes_combo_t::load3pts(const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0)
{
  const std::string plot3ptsPath=combine("%s/plots/%s/3pts_corr/",ens.dirPath.c_str(),mesComboTag.c_str());
  if(not dir_exists(plot3ptsPath))
    mkdir(plot3ptsPath);
  
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
	
	contr.ave_err().write(combine("%s/o%smuGPo-gs_%s_pol_%zu_gamma_%zu.xmg",plot3ptsPath.c_str(),VA_tag[iVA],decKinTag(iMoms,iMomt,iMom0).c_str(),iPol,iGamma));
      }
  corr/=2*sqrt(2);
  
  corr/=-((double)L*L*L);
  
  corr.symmetrize(par);
  corr.ave_err().write(combine("%s/o%smuGPo-gs_%s.xmg",plot3ptsPath.c_str(),VA_tag[iVA],decKinTag(iMoms,iMomt,iMom0).c_str()));
  
  return corr;
}

void permes_combo_t::load3pts()
{
  cout<<"Loading 3pts correlators"<<endl;
  
  //! Coefficient to combine the two insertions
  const double coeff[2]={+1.0,-1.0};
  
  for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
    {
      const vector<size_t> c=ens.indDecKin(iDecKin);
      
      for(int iVA=0;iVA<2;iVA++)
	{
	  corrPX[iVA][iDecKin]=
	    load3pts(iVA,iMs,iMt,c[0],c[1],c[2])*eT+
	    coeff[iVA]*
	    load3pts(iVA,iMt,iMs,c[0],c[1],c[2])*eS;
	  // cout<<eT<<" "<<eS<<" "<<corrPX[iVA][iDecKin][10].ave_err()<<endl;
	}
    }
}

//! Prepare the 3pts normalization
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

// void fuffa()
  // {
    
    // /////////////////////////////////////////////////////////////////
    
    // vector<size_t> tMin(ind3ptsKin.max(),0);
    // vector<size_t> tMax(ind3ptsKin.max(),0);
    
    // for(int iVA=0;iVA<2;iVA++)
    //   {
    // 	for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    // 	  if(consider[i3ptsKin])
    // 	    {
    // 	      grace_file_t outplot(combine("plots/3pts_%s_range_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
    // 	      const djvec_t eff=effective_mass(corr[iVA][i3ptsKin],T/2,0);
    // 	      forward_derivative(eff).ave_err().write(combine("plots/3pts_%s_eff_eff_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
    // 	      const djvec_t y=eff-dE[i3ptsKin];
    // 	      outplot.write_vec_ave_err(y.ave_err());
	      
    // 	      //compute the compatibility level
    // 	      vector<double> err(T/2),comp(T/2);
    // 	      grace_file_t outcomp(combine("plots/3pts_%s_comp_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
    // 	      for(size_t t=T/12;t<T/2-T/12;t++)
    // 		{
    // 		  const djack_t e=eff[t]-dE[i3ptsKin];
    // 		  if(not isnan(e.err()))
    // 		    {
    // 		      err[t]=e.err();
    // 		      comp[t]=fabs(e.ave()/e.err());
    // 		      outcomp.write_xy(t,comp[t]);
    // 		    }
    // 		}
	      
    // 	      if(hasSymm[i3ptsKin])
    // 		{
    // 		  const int i3ptsKinSymm=symmOf[i3ptsKin];
    // 		  cout<<dE[i3ptsKin].ave_err()<<" "<<dE[i3ptsKinSymm].ave_err()<<endl;
    // 		  grace_file_t out(combine("plots/3pts_%s_symm_rat_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
    // 		  const djvec_t s=-corr[iVA][i3ptsKin]/corr[iVA][symmOf[i3ptsKin]]-1.0;
    // 		  out.write_vec_ave_err(s.ave_err());
    // 		}
	      
    // 	      //plot the compatibility level
    // 	      outcomp.new_data_set();
    // 	      outcomp.set_all_colors(grace::BLUE);
    // 	      const double &err_min=*min_element(err.begin()+T/12,err.end()-T/12);
    // 	      const double &comp_min=*min_element(comp.begin()+T/12,comp.end()-T/12);
    // 	      cout<<"comp min: "<<comp_min<<endl;
    // 	      cout<<"err min: "<<err_min<<endl;
    // 	      for(size_t t=0;t<T/2;t++)
    // 		outcomp.write_xy(t,err[t]/err_min);
	      
    // 	      //! Number of stddev for compatibility
    // 	      const double nDev=std::max(comp_min*1.1,1.0);
    // 	      cout<<"nDev: "<<nDev<<endl;
	      
    // 	      vector<size_t> cl;
    // 	      bool incl=false;
    // 	      for(size_t t=T/12;t<T/2-T/12;t++)
    // 		{
    // 		  const djack_t e=eff[t]-dE[i3ptsKin];
    // 		  const double c=fabs(e.ave()/e.err());
    // 		  const bool newIncl=(c<nDev);
		  
    // 		  if(incl and not newIncl) cl.push_back(t-1);
    // 		  if((not incl) and newIncl) cl.push_back(t);
		  
    // 		  incl=newIncl;
    // 		}
	      
    // 	      if(incl)
    // 		cl.push_back(T/2-T/12);
	      
    // 	      if(cl.size()==0)
    // 		CRASH("No compatible point for kinematic %s",ind3ptsKin.descr(i3ptsKin).c_str());
    // 	      else
    // 		{
    // 		  cout<<"NRanges: "<<cl.size()/2<<"  ";
    // 		  for(size_t iCl=0;iCl<cl.size()/2;iCl++)
    // 		    cout<<"["<<cl[iCl*2]<<";"<<cl[iCl*2+1]<<"] ";
    // 		  cout<<endl;
    // 		}
	      
    // 	      if(cl.size()%2)
    // 		CRASH("Size %d of %d should be even",(int)i3ptsKin,(int)cl.size());
	      
    // 	      // for(size_t iCl=0;iCl<cl.size()/2;iCl++)
    // 	      //   outplot.write_line(cl[iCl*2]-0.5,cl[iCl*2+1]+0.5,0.0);
	      
    // 	      cout<<"Merging "<<i3ptsKin<<endl;
    // 	      size_t iCl=0;
    // 	      int nMerged=0;
    // 	      do
    // 		{
    // 		  while(cl.size()/2>1 and iCl<cl.size()/2)
    // 		    {
    // 		      auto begThis=cl.begin()+iCl;
    // 		      auto endThis=cl.begin()+iCl+1;
    // 		      auto begNext=cl.begin()+iCl+2;
    // 		      auto endNext=cl.begin()+iCl+3;
		      
    // 		      const size_t gap=((endNext-begNext)+(endThis-begThis))/2;
		      
    // 		      if(*(endThis+1)<=(*endThis)+gap)
    // 			{
    // 			  cout<<"Merged ["<<*begThis<<";"<<*endThis<<"] with ["<<*begNext<<";"<<*endNext<<"], result";
    // 			  cl.erase(endThis,begNext+1);
    // 			  cout<<" ["<<*begThis<<";"<<*endThis<<"]"<<endl;
			  
    // 			  nMerged++;
    // 			}
    // 		      else
    // 			{
    // 			  cout<<"Not merged ["<<*begThis<<";"<<*endThis<<"] with ["<<*begNext<<";"<<*endNext<<"]"<<endl;
    // 			  iCl+=2;
    // 			}
    // 		    }
    // 		  cout<<"NMerged: "<<nMerged<<endl;
    // 		  cout<<endl;
    // 		}
    // 	      while(nMerged!=0);
	      
    // 	      // for(size_t iCl=0;iCl<cl.size()/2;iCl++)
    // 	      //   outplot.write_constant_band(cl[iCl*2]-0.5,cl[iCl*2+1]+0.5,djack_t(dE[i3ptsKin]-dE[i3ptsKin].err()*3));
	      
    // 	      outplot.set_title(to_string(i3ptsKin));
	      
    // 	      //select tmin/max taking largest interval
    // 	      for(size_t iCl=0;iCl<cl.size()/2;iCl++)
    // 		if(tMax[i3ptsKin]-tMin[i3ptsKin]<cl[iCl*2+1]-cl[iCl*2])
    // 		  {
    // 		    tMax[i3ptsKin]=cl[iCl*2+1];
    // 		    tMin[i3ptsKin]=cl[iCl*2];
    // 		  }
	      
    // 	      cout<<i3ptsKin<<" "<<ind3ptsKin.descr(i3ptsKin)<<": range ["<<tMin[i3ptsKin]<<";"<<tMax[i3ptsKin]<<"]"<<endl;
	      
    // 	      /////////////////////////////////////////////////////////////////
	      
    // 	      vector<double> a(y.size()),e(y.size());
    // 	      for(size_t i=0;i<y.size();i++)
    // 		{
    // 		  a[i]=y[i].ave();
    // 		  e[i]=y[i].err();
    // 		}
	      
    // 	      const size_t firstCompDef=T/12;
    // 	      const size_t lastCompDef=comp.size()-T/12;
    // 	      size_t firstComp=firstCompDef;
    // 	      size_t lastComp=lastCompDef;
	      
    // 	      //! loop until a decent compatibility range is obtained
    // 	      double nSig=1.0;
    // 	      do
    // 		{
    // 		  do firstComp++;
    // 		  while(comp[firstComp]>nSig and firstComp<comp.size());
		  
    // 		  do lastComp--;
    // 		  while(comp[lastComp]>nSig and firstComp>0);
		  
    // 		  nSig*=1.1;
    // 		}
    // 	      while(lastComp==lastCompDef or firstComp==firstCompDef);
	      
    // 	      const size_t D=std::min(5,(int)(lastComp-firstComp));
    // 	      int tMinBest=0,tMaxBest=0;
    // 	      double pMax=0;
	      
    // 	      for(size_t tMin=firstComp;tMin<lastComp-D;tMin++)
    // 		for(size_t tMax=tMin+D;tMax<lastComp;tMax++)
    // 		  {
    // 		    const size_t d=tMax-tMin+1;
		    
    // 		    double Ch2Unc=0;
		    
    // 		    for(size_t i=tMin;i<=tMax;i++)
    // 		      Ch2Unc+=sqr(comp[i]);
		    
    // 		    Matrix<double,Dynamic,Dynamic> c(d,d);
    // 		    for(size_t i=0;i<d;i++)
    // 		      for(size_t j=0;j<d;j++)
    // 			c(i,j)=cov(y[i+tMin],y[j+tMin]);
		    
    // 		    Matrix<double,Dynamic,Dynamic> cInv=c.inverse();
    // 		    double Ch2Corr=0;
    // 		    for(size_t i=0;i<d;i++)
    // 		      for(size_t j=0;j<d;j++)
    // 			Ch2Corr+=a[i+tMin]*a[j+tMin]*cInv(i,j);
		    
    // 		    //Ch2Corr/=tMax-tMin+1;
		    
    // 		    if(Ch2Corr<1000 and Ch2Corr>0)
    // 		      {
    // 			const double p=ch2Distr(Ch2Corr,d);
    // 			cout<<tMin<<" "<<tMax<<" "<<Ch2Corr<<" "<<Ch2Unc<<" "<<p<<endl;
    // 			if(p>pMax)
    // 			  {
    // 			    pMax=p;
    // 			    tMinBest=tMin;
    // 			    tMaxBest=tMax;
    // 			  }
    // 		      }
    // 		  }
	      
    // 	      if(tMinBest==0 and tMaxBest==0)
    // 		{
    // 		  tMinBest=T/4-D/2;
    // 		  tMaxBest=T/4+D/2;
    // 		}
	      
    // 	      cout<<"Best: "<<ind3ptsKin.descr(i3ptsKin)<<" ["<<tMinBest<<":"<<tMaxBest<<"] "<<pMax<<" within ["<<firstComp<<":"<<lastComp<<"]"<<endl;
    // 	      outplot.write_line([](double){return 0;},tMinBest-0.5,tMaxBest+0.5);
	      
    // 	      tMax[i3ptsKin]=tMaxBest;
    // 	      tMin[i3ptsKin]=tMinBest;
	      
    // 	      /////////////////////////////////////////////////////////////////
    // 	    }
	
permes_combo_t& permes_combo_t::fit3pts(const char* fitTag)
{
  cout<<"Fitting 3pts correlators"<<endl;
  
  decPlotsPath=combine("%s/plots/%s/3pts_fit/%s",ens.dirPath.c_str(),mesComboTag.c_str(),fitTag);
  mkdir(decPlotsPath);
  
  for(int iVA=0;iVA<2;iVA++)
    for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
      if(ens.considerDec[iDecKin] // and hasSymm[iDecKin]
	 )
	{
	  const int iMesKin=ens.iMesKinOfDecKin[iDecKin];
	  
	  const vector<size_t> cDec=ens.indDecKin(iDecKin);
	  
	  const djvec_t y=corrPX[iVA][iDecKin]*normaliz[iDecKin];
	  const djack_t H=constant_fit(y,ens.T/4-1// +0*tMin[iDecKin]
				       ,ens.T/4+1// +0*tMax[iDecKin]
				       ,combine("%s/c%s_%s.xmg",decPlotsPath.c_str(),VA_tag[iVA],decKinTag(cDec[0],cDec[1],cDec[2]).c_str()));
	
  	if(iVA==1)
  	  ff[iVA][iDecKin]=(H-fPbare[iMesKin]*(eT-eS))/PKdec[iDecKin];
  	else
  	  ff[iVA][iDecKin]=H*E[iMesKin]/(ens.Eg[iDecKin]*ens.pMes[iMesKin]-E[iMesKin]*ens.kHatDec[iDecKin]);
	}
  
  return *this;
}

permes_combo_t& permes_combo_t::plotFf()
{
  cout<<"Plotting ff"<<endl;
  
  for(int iVA=0;iVA<2;iVA++)
    {
      grace_file_t ffPlot(combine("%s/plots/%s/ff_%s.xmg",ens.dirPath.c_str(),mesComboTag.c_str(),VA_tag[iVA]));
      ffPlot.set_line_style(grace::NO_LINE);
      
      for(size_t iDecKin=0;iDecKin<ens.indDecKin.max();iDecKin++)
	if(ens.considerDec[iDecKin] // and hasSymm[iDecKin]
	   )
	  {
	    ffPlot<<"# "<<iDecKin<<" "<<ens.indDecKin.descr(iDecKin)<<"\n";
	    ffPlot.write_ave_err(X[iDecKin].ave(),ff[iVA][iDecKin].ave_err());
	  }
    }
  
  return *this;
}

#endif
