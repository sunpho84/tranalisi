#ifndef _PERMES_COMBO_HPP
#define _PERMES_COMBO_HPP

#include <REPH/perens.hpp>

//! Holds all info for a combination of quarks
struct permes_combo_t
{
  //! Reference ensemble
  const perens_t& ens;
  
  //! Index of the spectator quark
  const size_t iMs;
  
  //! Index of the forward line quark
  const size_t iMt;
  
  //! Combination relative to the meson combination
  const std::string mesComboTag;
  
  //! Path where to store mesons plots
  const std::string mesPlotsPath;
  
  //! Charge of the spectator quark
  const double eS;
  
  //! Charge of the forward line quark
  const size_t eT;
  
  //! Pseudoscalar coupling
  djvec_t ZP;
  
  //! Axial coupling
  djvec_t ZA;
  
  //! Decay constant taken from Pseudoscalar current
  djvec_t fP;
  
  //! Decay constant taken from Axial current
  djvec_t fPbare;
  
  //! Energy
  djvec_t E;
  
  //! Effective mass
  vector<djvec_t> eEff;
  
  //! Load the PP correlation function
  djvec_t load2ptsPP(const size_t iMoms,const size_t iMomt);
  
  //! Load the AP correlation function
  djvec_t load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma);
  
  //! Load the three points correlation functions
  djvec_t load3pts(const size_t iVA,const size_t iMoms,const size_t iMomt,const size_t iMom0);
  
  //! Computes the axial couplings
  void computeAxialPseudoCouplings();
  
  //! Plot the dispersion relation
  void plotDispRel() const;
  
  // void fuffa()
  // {
  //   //set proxies
  //   const vector<std::array<double,3>>& moms=ens.moms;
  //   const size_t& nMoms=ens.nMoms;
  //   const size_t& T=ens.T;
  //   const size_t& L=ens.L;
    
    /////////////////////////////////////////////////////////////////
    
    // const index_t ind3ptsKin({{"iMoms",nMoms},{"iMomt",nMoms},{"iMom0",nMoms}});
    
    // vector<vector<djvec_t>> corr(2,vector<djvec_t>(ind3ptsKin.max(),djvec_t{(size_t)T/2+1}));
    // vector<djvec_t> normaliz(ind3ptsKin.max(),djvec_t{(size_t)T/2+1});
    // djvec_t dE(ind3ptsKin.max());
    // djvec_t PK(ind3ptsKin.max());
    
    // vector<bool> consider(ind3ptsKin.max());
    // vector<bool> hasSymm(ind3ptsKin.max());
    // vector<int> symmOf(ind3ptsKin.max(),-1);
    // for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    //   {
    // 	const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
    // 	const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
    // 	consider[i3ptsKin]=(iMomt!=iMom0);
    // 	hasSymm[i3ptsKin]=(iMoms==iMom0);
    // 	if(hasSymm[i3ptsKin])
    // 	  symmOf[i3ptsKin]=ind3ptsKin({iMomt,iMoms,iMomt});
    //   }
    
    // vector<int> iMesOf3ptsKin(ind3ptsKin.max());
    // vector<double> P(ind3ptsKin.max()),Eg(ind3ptsKin.max()),k(ind3ptsKin.max()),khat(ind3ptsKin.max());
    // for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    //   {
    // 	const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
    // 	const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
    // 	const size_t iMes=iMesOf3ptsKin[i3ptsKin]=ens.indMesKin({iMoms,iMom0});
    // 	P[i3ptsKin]=2*M_PI*(moms[iMom0][2]-moms[iMoms][2])/L;
    // 	k[i3ptsKin]=2*M_PI*(moms[iMom0][2]-moms[iMomt][2])/L;
    // 	khat[i3ptsKin]=2*sin(k[i3ptsKin]/2);
	
    // 	// const djack_t E=latt_en_1D(M,P);
    // 	Eg[i3ptsKin]=2*asinh(fabs(khat[i3ptsKin])/2);
    // 	const double EgT=sinh(Eg[i3ptsKin])*(1-exp(-T*Eg[i3ptsKin]));
    // 	dE[i3ptsKin]=E[iMes]-Eg[i3ptsKin];
    // 	PK[i3ptsKin]=E[iMes]*Eg[i3ptsKin]-P[i3ptsKin]*khat[i3ptsKin];
	
    // 	cout<<iMoms<<" "<<iMomt<<" "<<iMom0<<endl;
    // 	cout<<" P: "<<P[i3ptsKin]<<endl;
    // 	cout<<" Pg: "<<k[i3ptsKin]<<endl;
    // 	cout<<" Pghat: "<<khat[i3ptsKin]<<endl;
    // 	cout<<" E: "<<E[iMes].ave_err()<<endl;
    // 	cout<<" Eg: "<<Eg[i3ptsKin]<<endl;
    // 	cout<<" EgT: "<<EgT<<endl;
    // 	cout<<" dE: "<<dE[i3ptsKin].ave_err()<<endl;
	
    // 	for(int iVA=0;iVA<2;iVA++)
    // 	  corr[iVA][i3ptsKin]=
    // 	    load3pts(iVA,iMoms,iMomt,iMom0)*eT+
    // 	    ((iVA==0)?+1.0:-1.0)*
    // 	    load3pts(iVA,iMoms,iMomt,iMom0)*eS;
    // 	// load3pts(iVA,iMs,iMt,iMoms,iMomt,iMom0);
	
    // 	for(size_t t=0;t<=T/2;t++)
    // 	  {
    // 	    const djack_t &A=E[iMes];
    // 	    //const djack_t &B=Eeff[iMes][std::min(std::min(t,T-t),T/2-1)];
    // 	    const djack_t &W=A;
    // 	    //cout<<t<<" "<<smart_print(A.ave_err())<<" "<<smart_print(B.ave_err())<<endl;
    // 	    normaliz[i3ptsKin][t]=4*W*EgT/(ZP[iMes]*exp(-t*W-(T/2-t)*Eg[i3ptsKin]));
    // 	  }
    //   }
    
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
	
    // 	grace_file_t H_plot(combine("plots/H_%s.xmg",VA_tag[iVA]));
    // 	grace_file_t ff_plot(combine("plots/ff_%s.xmg",VA_tag[iVA]));
    // 	ff_plot.set_line_style(grace::NO_LINE);
	
    // 	for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    // 	  if(consider[i3ptsKin] // and hasSymm[i3ptsKin]
    // 	     )
    // 	    {
    // 	      const int iMes=iMesOf3ptsKin[i3ptsKin];
	      
    // 	      const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
	      
    // 	      const djvec_t y=corr[iVA][i3ptsKin]*normaliz[i3ptsKin];
    // 	      const djack_t H=constant_fit(y,// T/4-1+0*
    // 					   tMin[i3ptsKin],// T/4+1+0*
    // 					   tMax[i3ptsKin],combine("plots/3pts_%s_fit_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
    // 	      const djack_t x=2*PK[i3ptsKin]/sqr(E[0]);
    // 	      djack_t f;
	      
    // 	      if(iVA==1)
    // 		f=(H-fPbare[iMes]*(eT-eS))/PK[i3ptsKin];
    // 	      else
    // 		f=H*E[0]/(Eg[i3ptsKin]*P[i3ptsKin]-E[iMes]*khat[i3ptsKin]);
	      
    // 	      H_plot.write_ave_err(1/PK[i3ptsKin].ave(),H.ave_err());
    // 	      ff_plot<<"# "<<i3ptsKin<<" "<<ind3ptsKin.descr(i3ptsKin)<<"\n";
    // 	      ff_plot.write_ave_err(x.ave(),f.ave_err());
    // 	    }
    //   }
  // }
  
  //! Constructor
  permes_combo_t(const perens_t& ens,const size_t& iMs,const size_t& iMt,const double& eS,const double& eT,const std::string tag) :
    ens(ens),
    iMs(iMs),
    iMt(iMt),
    mesComboTag(combine("Mes_iMs%zu" "_" "iMt%zu" "_%s",iMs,iMt,tag.c_str())),
    mesPlotsPath(combine("%s/plots/%s",ens.dirPath.c_str(),mesComboTag.c_str())),
    eS(eS),
    eT(eT)
  {
    mkdir(mesPlotsPath);
    
    ZA.resize(ens.nMesKin);
    ZP.resize(ens.nMesKin);
    fP.resize(ens.nMesKin);
    fPbare.resize(ens.nMesKin);
    E.resize(ens.nMesKin);
    eEff.resize(ens.nMesKin);
    
    computeAxialPseudoCouplings();
    plotDispRel();
  }
};

#endif
