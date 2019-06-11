#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

//! Incapsulate all info and operations
class perens_t
{
  int L,T,spatVol;
  size_t nMass;
  vector<double> mass;
  size_t nMoms;
  vector<array<double,3>> moms;
  
  djvec_t load2ptsPP(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt);
  
  djvec_t load2ptsAP(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iGamma);
  
  djvec_t load3pts(const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0)
  {
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
	  
	  contr.ave_err().write(combine("plots/o%smuGPo-gs_%s.xmg",VA_tag[iVA],ind.descr(i).c_str()));
	}
    corr/=2*sqrt(2);
    
    corr/=-(L*L*L);
    
    corr.symmetrize(par);
    const size_t iave=ind_ave({iMs,iMt,iMoms,iMomt,iMom0});
    corr.ave_err().write(combine("plots/o%smuGPo-gs_%s.xmg",VA_tag[iVA],ind_ave.descr(iave).c_str()));
    
    return corr;
  }
  
  void getAxialPseudoCouplings(djack_t& ZP,djack_t& ZA,djack_t& E,
			       const djvec_t& C_PP,const djvec_t& C_A0P,const djvec_t& C_A3P,const double P,const size_t tMin,const size_t tMax,const size_t iMom1,const size_t iMom2)
  {
    if(iMom1!=iMom2)
      {
	//perform a preliminary fit
	djack_t ZPP,ZA0P,ZA3P=1.0;
	two_pts_fit(ZPP,E,C_PP,T/2,tMin,tMax,combine("plots/CPP_%zu_%zu.xmg",iMom1,iMom2),"",+1);
	two_pts_fit(ZA0P,E,C_A0P,T/2,tMin,tMax-1,combine("plots/CA0P_%zu_%zu.xmg",iMom1,iMom2),"",-1);
	two_pts_fit(ZA3P,E,C_A3P,T/2,tMin,tMax,combine("plots/CA3P_%zu_%zu.xmg",iMom1,iMom2),"",+1);
	ZP=sqrt(ZPP);
	const djack_t M=sqrt(E*E-P*P);
	const djack_t ZA0=ZA0P/ZP,ZA3=ZA3P/ZP;
	ZA=(ZA0*E-ZA3*P)/M;
	
	//parameters to fit
	minimizer_pars_t pars;
	pars.add("E",E.ave(),E.err());
	pars.add("ZA",ZA.ave(),ZA.err());
	pars.add("ZP",ZP.ave(),ZP.err());
	
	//! fit for real
	size_t iel=0;
	auto x=vector_up_to<double>(C_PP.size());
	multi_ch2_t<djvec_t> two_pts_fit_obj({x,x,x},{tMin,tMin,tMin},{tMax,tMax-1,tMax},{C_PP,C_A0P,C_A3P},
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
	    E[iel]=par_min[0];
	    ZA[iel]=par_min[1];
	    ZP[iel]=par_min[2];
	  }
	
	write_constant_fit_plot(combine("plots/CPP_fit_%zu_%zu.xmg",iMom1,iMom2),tMin,tMax,E,effective_mass(C_PP,T/2,+1));
	write_constant_fit_plot(combine("plots/CA0P_fit_%zu_%zu.xmg",iMom1,iMom2),tMin,tMax,E,effective_mass(C_A0P,T/2,-1));
	write_constant_fit_plot(combine("plots/CA3P_fit_%zu_%zu.xmg",iMom1,iMom2),tMin,tMax,E,effective_mass(C_A3P,T/2,+1));
      }
    else
      {
	two_pts_SL_fit(ZP,ZA,E,C_A0P,C_PP,T/2,tMin,tMax-1,combine("plots/CPP_A0P_fit_%zu_%zu.xmg",iMom1,iMom2),-1,+1);
	ZA/=E;
      }
  }
  
public:
  
  perens_t()
  {
    raw_file_t fin("jacks/input.txt","r");
    
    L=fin.read<size_t>("L");
    spatVol=L*L*L;
    T=fin.read<size_t>("T");
    
    nMass=fin.read<size_t>("NMass");
    mass.resize(nMass);
    for(size_t iMass=0;iMass<nMass;iMass++)
      mass[iMass]=fin.read<double>();
    
    nMoms=fin.read<size_t>("NMoms");
    moms.resize(nMoms);
    for(size_t iMom=0;iMom<nMoms;iMom++)
      for(size_t mu=0;mu<3;mu++)
	moms[iMom][mu]=fin.read<double>();
    
    /////////////////////////////////////////////////////////////////
    
    // grace_file_t fit2ptsPlot("plots/fit2pts.xmg");
    grace_file_t ZPPlot("plots/ZP.xmg");
    grace_file_t ZAPlot("plots/ZA.xmg");
    
    const size_t iMs=1;
    const size_t iMt=3;
    
    const double es=-1.0/3;
    const double et=+2.0/3;
    
    // const size_t iMs=3;
    // const size_t iMt=1;
    
    // const double es=+2.0/3;
    // const double et=-1.0/3;
    
    const index_t indMesKin({{"mom1",nMoms},{"mom2",nMoms}});
    const size_t nMesKin=indMesKin.max();
    vector<double> Pmes(nMesKin);
    djvec_t E(nMesKin),ZP(nMesKin),ZA(nMesKin);
    vector<djvec_t> Eeff(nMesKin,djvec_t(T));
    
    double pMax=0;
    for(size_t iMom1=0;iMom1<nMoms;iMom1++)
      for(size_t iMom2=iMom1;iMom2<nMoms;iMom2++)
	{
	  const size_t i1=iMom1*nMoms+iMom2;
	  const size_t i2=iMom2*nMoms+iMom1;
	  
	  const int tMin=12,tMax=T/2+1;
	  
	  const djvec_t corrPP=load2ptsPP(iMs,iMt,iMom1,iMom2);
	  corrPP.ave_err().write(combine("plots/2pts_PP_corr_%d.xmg",i1));
	  Eeff[i1]=Eeff[i2]=effective_mass(corrPP);
	  const djvec_t corrA0P=load2ptsAP(iMs,iMt,iMom1,iMom2,0);
	  corrA0P.ave_err().write(combine("plots/2pts_AP0_corr_%d.xmg",i1));
	  
	  const djvec_t corrA3P=load2ptsAP(iMs,iMt,iMom1,iMom2,3);
	  corrA3P.ave_err().write(combine("plots/2pts_AP3_corr_%d.xmg",i1));
	  
	  // djack_t Z2,EFit;
	  // two_pts_fit(Z2,EFit,corrPP,T/2,tMin,tMax);
	  // Z[i1]=Z[i2]=sqrt(Z2);
	  // E[i1]=E[i2]=EFit;
	  
	  // fit2ptsPlot.write_vec_ave_err(effective_mass(corrPP).ave_err());
	  // fit2ptsPlot.write_constant_band(tMin,tMax,EFit);
	  
	  double P=2*M_PI*(moms[iMom2][2]-moms[iMom1][2])/L;
	  Pmes[i1]=+P;
	  Pmes[i2]=-P;
	  
	  pMax=std::max(pMax,fabs(P));
	  
	  //ZPPlot.write_ave_err(EFit.ave(),Z[i1].ave_err());
	  
	  getAxialPseudoCouplings(ZP[i1],ZA[i1],E[i1],corrPP,corrA0P,corrA3P,Pmes[i1],tMin,tMax,iMom1,iMom2);
	  ZA[i2]=ZA[i1];
	  ZP[i2]=ZP[i1];
	  E[i2]=E[i1];
	  
	  ZAPlot.write_ave_err(E[i1].ave(),ZA[i1].ave_err());
	  ZPPlot.write_ave_err(E[i1].ave(),ZP[i1].ave_err());
	}
    
    const djack_t fP=ZP[0]*(mass[iMs]+mass[iMt])/sqr(E[0]);
    const djack_t fP_bis=ZA[0];
    cout<<"fP: "<<fP.ave_err()<<endl;
    cout<<"fP_bis: "<<fP_bis.ave_err()<<endl;
    const djack_t Zv=fP/fP_bis;
    cout<<"Zv: "<<Zv.ave_err()<<endl;
    
    grace_file_t dispRel("plots/dispRel.xmg");
    dispRel.write_vec_ave_err(Pmes,E.ave_err());
    
    const djack_t &M=E[0];
    
    dispRel.write_polygon([M](double x){return latt_en_1D(M,x);},0,pMax,grace::GREEN);
    dispRel.write_polygon([M](double x){return cont_en_1D(M,x);},0,pMax,grace::VIOLET);
    
    /////////////////////////////////////////////////////////////////
    
    const index_t ind3ptsKin({{"iMoms",nMoms},{"iMomt",nMoms},{"iMom0",nMoms}});
    
    vector<vector<djvec_t>> corr(2,vector<djvec_t>(ind3ptsKin.max(),djvec_t{(size_t)T/2+1}));
    vector<djvec_t> normaliz(ind3ptsKin.max(),djvec_t{(size_t)T/2+1});
    djvec_t dE(ind3ptsKin.max());
    djvec_t PK(ind3ptsKin.max());
    
    vector<bool> consider(ind3ptsKin.max());
    vector<bool> hasSymm(ind3ptsKin.max());
    vector<int> symmOf(ind3ptsKin.max(),-1);
    for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
      {
	const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
	const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
	consider[i3ptsKin]=(iMomt!=iMom0);
	hasSymm[i3ptsKin]=(iMoms==iMom0);
	if(hasSymm[i3ptsKin])
	  symmOf[i3ptsKin]=ind3ptsKin({iMomt,iMoms,iMomt});
      }
    
    vector<int> iMesOf3ptsKin(ind3ptsKin.max());
    vector<double> P(ind3ptsKin.max()),Eg(ind3ptsKin.max()),k(ind3ptsKin.max()),khat(ind3ptsKin.max());
    for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
      {
	const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
	const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
	const size_t iMes=iMesOf3ptsKin[i3ptsKin]=indMesKin({iMoms,iMom0});
	P[i3ptsKin]=2*M_PI*(moms[iMom0][2]-moms[iMoms][2])/L;
	k[i3ptsKin]=2*M_PI*(moms[iMom0][2]-moms[iMomt][2])/L;
	khat[i3ptsKin]=2*sin(k[i3ptsKin]/2);
	
	// const djack_t E=latt_en_1D(M,P);
	Eg[i3ptsKin]=2*asinh(fabs(khat[i3ptsKin])/2);
	const double EgT=sinh(Eg[i3ptsKin])*(1-exp(-T*Eg[i3ptsKin]));
	dE[i3ptsKin]=E[iMes]-Eg[i3ptsKin];
	PK[i3ptsKin]=E[iMes]*Eg[i3ptsKin]-P[i3ptsKin]*khat[i3ptsKin];
	
	cout<<iMoms<<" "<<iMomt<<" "<<iMom0<<endl;
	cout<<" P: "<<P[i3ptsKin]<<endl;
	cout<<" Pg: "<<k[i3ptsKin]<<endl;
	cout<<" Pghat: "<<khat[i3ptsKin]<<endl;
	cout<<" E: "<<E[iMes].ave_err()<<endl;
	cout<<" Eg: "<<Eg[i3ptsKin]<<endl;
	cout<<" EgT: "<<EgT<<endl;
	cout<<" dE: "<<dE[i3ptsKin].ave_err()<<endl;
	
	for(int iVA=0;iVA<2;iVA++)
	  corr[iVA][i3ptsKin]=
	    load3pts(iVA,iMs,iMt,iMoms,iMomt,iMom0)*et+
	    ((iVA==0)?+1.0:-1.0)*
	    load3pts(iVA,iMt,iMs,iMoms,iMomt,iMom0)*es;
	// load3pts(iVA,iMs,iMt,iMoms,iMomt,iMom0);
	
	for(int t=0;t<=T/2;t++)
	  {
	    const djack_t &A=E[iMes];
	    //const djack_t &B=Eeff[iMes][std::min(std::min(t,T-t),T/2-1)];
	    const djack_t &W=A;
	    //cout<<t<<" "<<smart_print(A.ave_err())<<" "<<smart_print(B.ave_err())<<endl;
	    normaliz[i3ptsKin][t]=4*W*EgT/(ZP[iMes]*exp(-t*W-(T/2-t)*Eg[i3ptsKin]));
	  }
      }
    
    /////////////////////////////////////////////////////////////////
    
    vector<size_t> tMin(ind3ptsKin.max(),0);
    vector<size_t> tMax(ind3ptsKin.max(),0);
    
    for(int iVA=0;iVA<2;iVA++)
      {
	for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
	  if(consider[i3ptsKin])
	    {
	      grace_file_t outplot(combine("plots/3pts_%s_range_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
	      const djvec_t eff=effective_mass(corr[iVA][i3ptsKin],T/2,0);
	      forward_derivative(eff).ave_err().write(combine("plots/3pts_%s_eff_eff_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
	      const djvec_t y=eff-dE[i3ptsKin];
	      outplot.write_vec_ave_err(y.ave_err());
	      
	      //compute the compatibility level
	      vector<double> err(T/2),comp(T/2);
	      grace_file_t outcomp(combine("plots/3pts_%s_comp_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
	      for(int t=T/12;t<T/2-T/12;t++)
		{
		  const djack_t e=eff[t]-dE[i3ptsKin];
		  if(not isnan(e.err()))
		    {
		      err[t]=e.err();
		      comp[t]=fabs(e.ave()/e.err());
		      outcomp.write_xy(t,comp[t]);
		    }
		}
	      
	      if(hasSymm[i3ptsKin])
		{
		  const int i3ptsKinSymm=symmOf[i3ptsKin];
		  cout<<dE[i3ptsKin].ave_err()<<" "<<dE[i3ptsKinSymm].ave_err()<<endl;
		  grace_file_t out(combine("plots/3pts_%s_symm_rat_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
		  const djvec_t s=-corr[iVA][i3ptsKin]/corr[iVA][symmOf[i3ptsKin]]-1.0;
		  out.write_vec_ave_err(s.ave_err());
		}
	      
	      //plot the compatibility level
	      outcomp.new_data_set();
	      outcomp.set_all_colors(grace::BLUE);
	      const double &err_min=*min_element(err.begin()+T/12,err.end()-T/12);
	      const double &comp_min=*min_element(comp.begin()+T/12,comp.end()-T/12);
	      cout<<"comp min: "<<comp_min<<endl;
	      cout<<"err min: "<<err_min<<endl;
	      for(int t=0;t<T/2;t++)
		outcomp.write_xy(t,err[t]/err_min);
	      
	      //! Number of stddev for compatibility
	      const double nDev=std::max(comp_min*1.1,1.0);
	      cout<<"nDev: "<<nDev<<endl;
	      
	      vector<size_t> cl;
	      bool incl=false;
	      for(int t=T/12;t<T/2-T/12;t++)
		{
		  const djack_t e=eff[t]-dE[i3ptsKin];
		  const double c=fabs(e.ave()/e.err());
		  const bool newIncl=(c<nDev);
		  
		  if(incl and not newIncl) cl.push_back(t-1);
		  if((not incl) and newIncl) cl.push_back(t);
		  
		  incl=newIncl;
		}
	      
	      if(incl)
		cl.push_back(T/2-T/12);
	      
	      if(cl.size()==0)
		CRASH("No compatible point for kinematic %s",ind3ptsKin.descr(i3ptsKin).c_str());
	      else
		{
		  cout<<"NRanges: "<<cl.size()/2<<"  ";
		  for(size_t iCl=0;iCl<cl.size()/2;iCl++)
		    cout<<"["<<cl[iCl*2]<<";"<<cl[iCl*2+1]<<"] ";
		  cout<<endl;
		}
	      
	      if(cl.size()%2)
		CRASH("Size %d of %d should be even",(int)i3ptsKin,(int)cl.size());
	      
	      // for(size_t iCl=0;iCl<cl.size()/2;iCl++)
	      //   outplot.write_line(cl[iCl*2]-0.5,cl[iCl*2+1]+0.5,0.0);
	      
	      cout<<"Merging "<<i3ptsKin<<endl;
	      size_t iCl=0;
	      int nMerged=0;
	      do
		{
		  while(cl.size()/2>1 and iCl<cl.size()/2)
		    {
		      auto begThis=cl.begin()+iCl;
		      auto endThis=cl.begin()+iCl+1;
		      auto begNext=cl.begin()+iCl+2;
		      auto endNext=cl.begin()+iCl+3;
		      
		      const size_t gap=((endNext-begNext)+(endThis-begThis))/2;
		      
		      if(*(endThis+1)<=(*endThis)+gap)
			{
			  cout<<"Merged ["<<*begThis<<";"<<*endThis<<"] with ["<<*begNext<<";"<<*endNext<<"], result";
			  cl.erase(endThis,begNext+1);
			  cout<<" ["<<*begThis<<";"<<*endThis<<"]"<<endl;
			  
			  nMerged++;
			}
		      else
			{
			  cout<<"Not merged ["<<*begThis<<";"<<*endThis<<"] with ["<<*begNext<<";"<<*endNext<<"]"<<endl;
			  iCl+=2;
			}
		    }
		  cout<<"NMerged: "<<nMerged<<endl;
		  cout<<endl;
		}
	      while(nMerged!=0);
	      
	      // for(size_t iCl=0;iCl<cl.size()/2;iCl++)
	      //   outplot.write_constant_band(cl[iCl*2]-0.5,cl[iCl*2+1]+0.5,djack_t(dE[i3ptsKin]-dE[i3ptsKin].err()*3));
	      
	      outplot.set_title(to_string(i3ptsKin));
	      
	      //select tmin/max taking largest interval
	      for(size_t iCl=0;iCl<cl.size()/2;iCl++)
		if(tMax[i3ptsKin]-tMin[i3ptsKin]<cl[iCl*2+1]-cl[iCl*2])
		  {
		    tMax[i3ptsKin]=cl[iCl*2+1];
		    tMin[i3ptsKin]=cl[iCl*2];
		  }
	      
	      cout<<i3ptsKin<<" "<<ind3ptsKin.descr(i3ptsKin)<<": range ["<<tMin[i3ptsKin]<<";"<<tMax[i3ptsKin]<<"]"<<endl;
	      
	      /////////////////////////////////////////////////////////////////
	      
	      vector<double> a(y.size()),e(y.size());
	      for(size_t i=0;i<y.size();i++)
		{
		  a[i]=y[i].ave();
		  e[i]=y[i].err();
		}
	      
	      const size_t firstCompDef=T/12;
	      const size_t lastCompDef=comp.size()-T/12;
	      size_t firstComp=firstCompDef;
	      size_t lastComp=lastCompDef;
	      
	      //! loop until a decent compatibility range is obtained
	      double nSig=1.0;
	      do
		{
		  do firstComp++;
		  while(comp[firstComp]>nSig and firstComp<comp.size());
		  
		  do lastComp--;
		  while(comp[lastComp]>nSig and firstComp>0);
		  
		  nSig*=1.1;
		}
	      while(lastComp==lastCompDef or firstComp==firstCompDef);
	      
	      const size_t D=std::min(5,(int)(lastComp-firstComp));
	      int tMinBest=0,tMaxBest=0;
	      double pMax=0;
	      
	      for(size_t tMin=firstComp;tMin<lastComp-D;tMin++)
		for(size_t tMax=tMin+D;tMax<lastComp;tMax++)
		  {
		    const size_t d=tMax-tMin+1;
		    
		    double Ch2Unc=0;
		    
		    for(size_t i=tMin;i<=tMax;i++)
		      Ch2Unc+=sqr(comp[i]);
		    
		    Matrix<double,Dynamic,Dynamic> c(d,d);
		    for(size_t i=0;i<d;i++)
		      for(size_t j=0;j<d;j++)
			c(i,j)=cov(y[i+tMin],y[j+tMin]);
		    
		    Matrix<double,Dynamic,Dynamic> cInv=c.inverse();
		    double Ch2Corr=0;
		    for(size_t i=0;i<d;i++)
		      for(size_t j=0;j<d;j++)
			Ch2Corr+=a[i+tMin]*a[j+tMin]*cInv(i,j);
		    
		    //Ch2Corr/=tMax-tMin+1;
		    
		    if(Ch2Corr<1000 and Ch2Corr>0)
		      {
			const double p=ch2Distr(Ch2Corr,d);
			cout<<tMin<<" "<<tMax<<" "<<Ch2Corr<<" "<<Ch2Unc<<" "<<p<<endl;
			if(p>pMax)
			  {
			    pMax=p;
			    tMinBest=tMin;
			    tMaxBest=tMax;
			  }
		      }
		  }
	      
	      if(tMinBest==0 and tMaxBest==0)
		{
		  tMinBest=T/4-D/2;
		  tMaxBest=T/4+D/2;
		}
	      
	      cout<<"Best: "<<ind3ptsKin.descr(i3ptsKin)<<" ["<<tMinBest<<":"<<tMaxBest<<"] "<<pMax<<" within ["<<firstComp<<":"<<lastComp<<"]"<<endl;
	      outplot.write_line([](double){return 0;},tMinBest-0.5,tMaxBest+0.5);
	      
	      tMax[i3ptsKin]=tMaxBest;
	      tMin[i3ptsKin]=tMinBest;
	      
	      /////////////////////////////////////////////////////////////////
	    }
	
	grace_file_t H_plot(combine("plots/H_%s.xmg",VA_tag[iVA]));
	grace_file_t ff_plot(combine("plots/ff_%s.xmg",VA_tag[iVA]));
	ff_plot.set_line_style(grace::NO_LINE);
	
	for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
	  if(consider[i3ptsKin] // and hasSymm[i3ptsKin]
	     )
	    {
	      const int iMes=iMesOf3ptsKin[i3ptsKin];
	      
	      const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
	      
	      const djvec_t y=corr[iVA][i3ptsKin]*normaliz[i3ptsKin];
	      const djack_t H=constant_fit(y,// T/4-1+0*
					   tMin[i3ptsKin],// T/4+1+0*
					   tMax[i3ptsKin],combine("plots/3pts_%s_fit_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
	      const djack_t x=2*PK[i3ptsKin]/sqr(E[0]);
	      djack_t f;
	      
	      if(iVA==1)
		f=(H-fP_bis*(et-es))/PK[i3ptsKin];
	      else
		f=H*E[0]/(Eg[i3ptsKin]*P[i3ptsKin]-E[iMes]*khat[i3ptsKin]);
	      
	      H_plot.write_ave_err(1/PK[i3ptsKin].ave(),H.ave_err());
	      ff_plot<<"# "<<i3ptsKin<<" "<<ind3ptsKin.descr(i3ptsKin)<<"\n";
	      ff_plot.write_ave_err(x.ave(),f.ave_err());
	    }
      }
  }
};

#endif
