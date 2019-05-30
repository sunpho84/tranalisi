#include <tranalisi.hpp>

int L,T,spatVol;
size_t nMass;
vector<double> mass;
size_t nMoms;
vector<array<double,3>> moms;


  //               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
  //               }                           |
  //     --- -> Q0 --- X --- -> Qt ---         |
  //    /                             \        |
  //   /                               \       |
  // P ---------- <- QS -------------V-A (t)
  

const char VA_tag[2][2]={"V","A"};

djvec_t load2ptsPP(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt)
{
  const index_t ind({{"iks",nMass},
		     {"ikt",nMass},
		     {"moms",nMoms},
		     {"momt",nMoms},
		     {"gamma",1},
		     {"reim",2}});
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,0,0});
  
  return read_djvec("jacks/oPPo-ss",T,i).symmetrized()/spatVol;
}

djvec_t load2ptsAP(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iGamma)
{
  const index_t ind({{"iks",nMass},
		     {"ikt",nMass},
		     {"moms",nMoms},
		     {"momt",nMoms},
		     {"gamma",4},
		     {"reim",2}});
  
  const size_t iReIm=(iGamma==0)?RE:IM;
  const int par=(iGamma==0)?-1:+1;
  const int sign=(iGamma==0)?-1:+1;
  
  const size_t i=ind({iMs,iMt,iMoms,iMomt,iGamma,iReIm});
  
  return sign*read_djvec("jacks/oAmuPo-ss",T,i).symmetrized(par)/spatVol;
}

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
  
  corr/=(L*L*L);

  corr.symmetrize(par);
  const size_t iave=ind_ave({iMs,iMt,iMoms,iMomt,iMom0});
  corr.ave_err().write(combine("plots/o%smuGPo-gs_%s.xmg",VA_tag[iVA],ind_ave.descr(iave).c_str()));
  
  return corr;
}

vector<double> compatibilityWith(const djvec_t& corr,const djack_t& c)
{
  const size_t n=corr.size();
  vector<double> out(n*n);
  
  vector<double> y(n);
  for(size_t i=0;i<n;i++)
    {
      const djack_t temp=corr[i]-c;
      const ave_err_t ae=temp.ave_err();
      if(not isnan(ae.ave()))
	y[i]=sqr(ae.ave()/ae.err());
      else
	y[i]+=1.0;
    }
  
  for(size_t iMin=0;iMin<n;iMin++)
    {
      double Ch2=0;
      for(size_t iMax=iMin+1;iMax<n-1;iMax++)
	{
	  Ch2+=y[iMax];
	  out[iMax+n*iMin]=Ch2;///(iMax-iMin+1);
	}
    }
  
  return out;
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

// void fix_by_markov()
//   {
//   vector<vector<double>> Ch2(ind3ptsKin.max());
//   vector<array<double,3>> p3pts(ind3ptsKin.max());
//   for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
//     {
//       const vector<size_t> c=ind3ptsKin(i3ptsKin);
//       for(size_t i=0;i<3;i++)
// 	p3pts[i3ptsKin][i]=moms[c[i]][2];
//       if(consider[i3ptsKin])
//       	Ch2[i3ptsKin]=compatibilityWith(effective_mass(corr[i3ptsKin],T/2,0),dE[i3ptsKin]);
//     }
  
//   vector<array<int,6>> neigh(ind3ptsKin.max());
//   vector<int> tMin(ind3ptsKin.max(),T/4-1);
//   vector<int> tMax(ind3ptsKin.max(),T/4+1);
//   for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
//     {
//       vector<size_t> c=ind3ptsKin(i3ptsKin);
//       for(size_t i=0;i<3;i++)
// 	{
// 	  if(c[i]!=0)
// 	    {
// 	      c[i]--;
// 	      neigh[i3ptsKin][i]=ind3ptsKin(c);
// 	      c[i]++;
// 	    }
// 	  else
// 	    neigh[i3ptsKin][i]=-1;
	  
// 	  if(c[i]!=nMoms-1)
// 	    {
// 	      c[i]++;
// 	      neigh[i3ptsKin][i+3]=ind3ptsKin(c);
// 	      c[i]--;
// 	    }
// 	  else
// 	    neigh[i3ptsKin][i+3]=-1;
// 	}
//     }
  
//   auto totCh2Fun=[&]()
// 		 {
// 		   double C=0;
// 		   int n=0;
// 		   for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
// 		     if(consider[i3ptsKin])
// 		       {
// 			 C+=Ch2[i3ptsKin][tMax[i3ptsKin]+T/2*tMin[i3ptsKin]];
// 			 n+=tMax[i3ptsKin]-tMin[i3ptsKin]+1;
			 
// 			 //cout<<" "<<i3ptsKin<<" ["<<tMax[i3ptsKin]<<";"<<tMin[i3ptsKin]<<"], "<<Ch2[i3ptsKin][tMax[i3ptsKin]+T/2*tMin[i3ptsKin]]<<endl;
// 		       }
		   
// 		   return C-n;
// 		 };
  
//   const int nR=000000;
//   mt19937 gen(1213);
//   double Ch2Old=totCh2Fun();
//   for(int iR=0;iR<nR;iR++)
//     {
//       const int i3ptsKin=uniform_int_distribution<int>(0,ind3ptsKin.max()-1)(gen);
//       if(consider[i3ptsKin])
// 	for(auto& tptr : {&tMin,&tMax})
// 	  {
// 	    int& t=(*tptr)[i3ptsKin];
	    
// 	    const int shift=uniform_int_distribution<int>(0,1)(gen)*2-1;
// 	    //cout<<"kin "<<i3ptsKin<<", ["<<tMin[i3ptsKin]<<";"<<tMax[i3ptsKin]<<"] shifting "<<t<<" of "<<shift<<endl;
	    
// 	    vector<size_t> c=ind3ptsKin(i3ptsKin);
// 	    int oldT=t;
// 	    if(t+shift<(int)T/2 and t>0)
// 	      {
// 		t+=shift;
// 		//cout<<" test "<<i3ptsKin<<", ["<<tMin[i3ptsKin]<<";"<<tMax[i3ptsKin]<<"]"<<endl;
		
// 		for(int d=0;d<6;d++)
// 		  {
// 		    const int n=neigh[i3ptsKin][d];
		    
// 		    if(n!=-1)
// 		      if(abs((*tptr)[n]-(*tptr)[i3ptsKin])>1)
// 			{
// 			  // cout<<(*tptr)[n]<<" too far for neighbor "<<n<<" in dir "<<d<<", "<<(*tptr)[i3ptsKin]<<endl;
// 			  t=oldT;
// 			}
// 		  }
		
// 		if(tMax[i3ptsKin]-tMin[i3ptsKin]<3)
// 		  {
// 		    // cout<<"Unordered tMin tMax "<<tMin[i3ptsKin]<<" "<<tMax[i3ptsKin]<<endl;
// 		    t=oldT;
// 		  }
// 		else
// 		  {
// 		    double Ch2New=totCh2Fun();
// 		    double delta=Ch2New-Ch2Old;
// 		    double pAcc=exp(-delta);
// 		    if(uniform_real_distribution<>(0,1)(gen)<pAcc)
// 		      {
// 			cout<<i3ptsKin<<" changed to ["<<tMin[i3ptsKin]<<";"<<tMax[i3ptsKin]<<"], new Ch2: "<<Ch2New<<endl;
// 			Ch2Old=Ch2New;
// 		      }
// 		    else
// 		      {
// 			// cout<<"Ch2New>Ch2Old: "<<Ch2New<<" "<<Ch2Old<<endl;
// 			t=oldT;
// 		      }
// 		  }
// 	      }
// 	    // else
// 	    //    cout<<"Cannot move "<<t<<"of "<<shift<<endl;
// 	  }
//     }
  
//   for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
//     if(consider[i3ptsKin])
//       {
// 	const vector<size_t> c=ind3ptsKin(i3ptsKin);
// 	cout<<c[0]<<" "<<c[1]<<" "<<c[2]<<"   ["<<tMin[i3ptsKin]<<":"<<tMax[i3ptsKin]<<"]"<<endl;
//       }
  
//   }

  // for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
  //   {
  //     const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
  //     const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
      
  //     if(consider[i3ptsKin])
  // 	{
  // 	  const djvec_t eff=effective_mass(corr[i3ptsKin],T/2,0);
  // 	  const djvec_t y=eff-dE[i3ptsKin];
	  
  // 	  const size_t s=y.size();
  // 	  //cout<<y.ave_err()<<endl;
  // 	  vector<double> a(s),e(s);
  // 	  for(size_t i=0;i<s;i++)
  // 	    {
  // 	      a[i]=y[i].ave();
  // 	      e[i]=y[i].err();
  // 	    }
	  
  // 	  size_t tMinBestCorr=0,tMaxBestCorr=s-1;
  // 	  size_t tMinBestUnc=0,tMaxBestUnc=s-1;
  // 	  double Ch2BestCorr=1e300;
  // 	  double Ch2BestUnc=1e300;
  // 	  const size_t D=5;
  // 	  for(size_t tMin=D;tMin<T/2-D;tMin++)
  // 	    for(size_t tMax=tMin+D;tMax<T/2-D;tMax++)
  // 	      {
  // 		const size_t d=tMax-tMin+1;
  // 		Matrix<double,Dynamic,Dynamic> c(d,d);
  // 		for(size_t i=0;i<d;i++)
  // 		  for(size_t j=0;j<d;j++)
  // 		    c(i,j)=cov(y[i+tMin],y[j+tMin]);
		
  // 		Matrix<double,Dynamic,Dynamic> cInv=c.inverse();
		
  // 		double Ch2Unc=0;
  // 		for(size_t i=tMin;i<=tMax;i++)
  // 		  Ch2Unc+=sqr(a[i]/e[i]);
  // 		Ch2Unc/=tMax-tMin+1;
		
  // 		double Ch2Corr=0;
  // 		for(size_t i=0;i<d;i++)
  // 		  for(size_t j=0;j<d;j++)
  // 		    Ch2Corr+=a[i+tMin]*a[j+tMin]*cInv(i,j);
		
  // 		Ch2Corr/=tMax-tMin+1;
		
  // 		//cout<<tMin<<" "<<tMax<<" "<<Ch2Corr<<" "<<Ch2Unc<<endl;
		
  // 		if(Ch2Corr>0 and Ch2Corr<100 and Ch2Corr<Ch2BestCorr)
  // 		  {
  // 		    Ch2BestCorr=Ch2Corr;
  // 		    tMinBestCorr=tMin;
  // 		    tMaxBestCorr=tMax;
  // 		  }
		
  // 		if(Ch2Unc<Ch2BestUnc)
  // 		  {
  // 		    Ch2BestUnc=Ch2Unc;
  // 		    tMinBestUnc=tMin;
  // 		    tMaxBestUnc=tMax;
  // 		  }
  // 	      }
	  
  // 	  cout<<"Best: "<<iMoms<<" "<<iMomt<<" "<<iMom0<<" "<<endl;
  // 	  cout<<"   tBestCorr: "<<tMinBestCorr<<" "<<tMaxBestCorr<<" Ch2: "<<Ch2BestCorr<<endl;
  // 	  cout<<"   tBestUnc:  "<<tMinBestUnc<<" "<<tMaxBestUnc<<" Ch2: "<<Ch2BestUnc<<endl;
	  
  // 	  grace_file_t outplot(combine("plots/3pts_range_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
  // 	  // outplot.write_vec_ave_err(eff.ave_err());
  // 	  // outplot.write_constant_band(tMinBestCorr,tMaxBestCorr,dE[i3ptsKin],grace::GREEN);
  // 	  // outplot.write_constant_band(tMinBestUnc,tMaxBestUnc,dE[i3ptsKin],grace::VIOLET);
	  
  // 	  vector<double> Dev(s);
  // 	  for(size_t i=0;i<s;i++)
  // 	    outplot.write_xy(i,fabs(y[i].ave())/y[i].err());

  // 	}
  //   }
  
  // for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
  //   {
  //     grace_file_t outplot(combine("plots/3pts_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
  //     outplot.write_vec_ave_err(effective_mass(corr[i3ptsKin],T/2,0).ave_err());
  //     outplot.write_constant_band(0,T/2,dE[i3ptsKin]);
      
  //     // djvec_t y=corr[i3ptsKin]*normaliz[i3ptsKin];
      
  //     // outplot.write_vec_ave_err(y.ave_err());
  //   }
  
	// cout<<"Best range: "<<tMinBest<<" "<<tMaxBest<<" "<<TotCh2Best<<endl;
  
  // for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
  //   if(consider[i3ptsKin])
  //     {
  // 	grace_file_t outplot(combine("plots/3pts_range_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
  // 	outplot.write_vec_ave_err(effective_mass(corr[i3ptsKin],T/2,0).ave_err());
  // 	outplot.write_constant_band(tMinBest,tMaxBest,dE[i3ptsKin]);
  //     }
  
  // const size_t s=T/2;
  // size_t tMinBestCorr=0,tMaxBestCorr=s-1;
  // size_t tMinBestUnc=0,tMaxBestUnc=s-1;
  // double Ch2BestCorr=1e300;
  // double Ch2BestUnc=1e300;
  // const size_t D=5;
  // for(size_t tMin=D;tMin<T/2-D;tMin++)
  //   for(size_t tMax=tMin+D;tMax<T/2-D;tMax++)
  //     {
  // 	const size_t d=tMax-tMin+1;
		
  // 	double Ch2Unc=0;
	
  // 	for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
  // 	  {
  // 	    const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
  // 	    const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
	    
  // 	    if(iMomt!=iMom0)
  // 	      {
  // 		const djvec_t eff=effective_mass(corr[i3ptsKin],T/2,0);
  // 		const djvec_t y=eff-dE[i3ptsKin];
		
  // 		for(size_t i=tMin;i<=tMax;i++)
  // 		  Ch2Unc+=sqr(a[i]/e[i]);
  // 		Ch2Unc/=tMax-tMin+1;
		
  // 		double Ch2Corr=0;
  // 		for(size_t i=0;i<d;i++)
  // 		  for(size_t j=0;j<d;j++)
  // 		    Ch2Corr+=a[i+tMin]*a[j+tMin]*cInv(i,j);
		
  // 		Ch2Corr/=tMax-tMin+1;
		
  // 		//cout<<tMin<<" "<<tMax<<" "<<Ch2Corr<<" "<<Ch2Unc<<endl;
		
  // 		if(Ch2Corr>0 and Ch2Corr<100 and Ch2Corr<Ch2BestCorr)
  // 		  {
  // 		    Ch2BestCorr=Ch2Corr;
  // 		    tMinBestCorr=tMin;
  // 		    tMaxBestCorr=tMax;
  // 		  }
		
  // 		if(Ch2Unc<Ch2BestUnc)
  // 		  {
  // 		    Ch2BestUnc=Ch2Unc;
  // 		    tMinBestUnc=tMin;
  // 		    tMaxBestUnc=tMax;
  // 		  }
  // 	      }
	  
  // 	  cout<<"Best: "<<iMoms<<" "<<iMomt<<" "<<iMom0<<" "<<endl;
  // 	  cout<<"   tBestCorr: "<<tMinBestCorr<<" "<<tMaxBestCorr<<" Ch2: "<<Ch2BestCorr<<endl;
  // 	  cout<<"   tBestUnc:  "<<tMinBestUnc<<" "<<tMaxBestUnc<<" Ch2: "<<Ch2BestUnc<<endl;
	  
  // 	  grace_file_t outplot(combine("plots/3pts_range_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
  // 	  // outplot.write_vec_ave_err(eff.ave_err());
  // 	  // outplot.write_constant_band(tMinBestCorr,tMaxBestCorr,dE[i3ptsKin],grace::GREEN);
  // 	  // outplot.write_constant_band(tMinBestUnc,tMaxBestUnc,dE[i3ptsKin],grace::VIOLET);
	  
  // 	  vector<double> Dev(s);
  // 	  for(size_t i=0;i<s;i++)
  // 	    outplot.write_xy(i,fabs(y[i].ave())/y[i].err());

  // 	}
  //   }
  

  // for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
  //   if(consider[i3ptsKin])
  //     {
  // 	const djack_t N=corr[i3ptsKin][T/4]/exp(-dE[i3ptsKin]*T/4);
	
  // 	//parameters to fit
  // 	minimizer_pars_t pars;
  // 	// pars.add("Z0",N.ave(),N.err());
  // 	pars.add("dE10",2,0.1);
  // 	pars.add("Z1",N.ave(),N.err());
  // 	pars.setlimits("dE10",1.0,10);
  // 	// pars.add("dE20",2,0.1);
  // 	// pars.add("Z2",N.ave(),N.err());
  // 	// pars.setlimits("dE20",1.0,10);
	
  // 	//! fit for real
  // 	size_t iel=0;
  // 	auto x=vector_up_to<double>(corr[i3ptsKin].size());
  // 	auto fw=[&](const auto& N,const auto& dE,const auto& dE10,const auto& Z1// ,const auto& dE20,const auto& Z2
  // 		    ,double x)
  // 		{
  // 		  return N*exp(-dE*x)+Z1*exp(-dE*dE10*x)// +Z2*exp(-dE*dE20*(T/2-x))
  // 		    ;
  // 		};
  // 	auto f=[&](const vector<double> &p,const double &x)
  // 	       {
  // 		 return fw(N[iel],dE[i3ptsKin][iel],p[0],p[1],// p[2],p[3],
  // 			   x);
  // 	       };
  // 	simple_ch2_t<djvec_t> fit_obj(vector_up_to<double>(corr[i3ptsKin].size()),tMin[i3ptsKin],tMax[i3ptsKin],corr[i3ptsKin],f,iel);
	
  // 	minimizer_t minimizer(fit_obj,pars);
	
  // 	djack_t Z0,dE10,Z1// ,dE20,Z2
  // 	  ;
  // 	for(iel=0;iel<njacks;iel++)
  // 	  {
  // 	    //minimize and print the result
  // 	    vector<double> par_min=minimizer.minimize();
  // 	    // Z0[iel]=par_min[0];
  // 	    dE10[iel]=par_min[0];
  // 	    Z1[iel]=par_min[1];
  // 	    // dE20[iel]=par_min[2];
  // 	    // Z2[iel]=par_min[3];
  // 	  }
	
  // 	cout<<"dE10 "<<dE10.ave_err()<<endl;
  // 	cout<<"Z1/N2: "<<djack_t(Z1/N).ave_err()<<endl;
  // 	// cout<<"dE20 "<<dE20.ave_err()<<endl;
  // 	// cout<<"Z2/N2: "<<djack_t(Z2/N).ave_err()<<endl;
	
  // 	auto der=[&](double x) -> djack_t
  // 		 {
  // 		   auto f1=[&](double x) -> djack_t
  // 			   {
  // 			     return N*exp(-dE[i3ptsKin]*x)+Z1*exp(-(dE[i3ptsKin]*dE10)*x)// +Z2*exp(-dE[i3ptsKin]*dE20*(T/2-x))
  // 			       ;
  // 			   };
		   
  // 		   const double dx=1e-3;
  // 		   return -log((f1(x+dx)/f1(x)))/dx;
  // 		 };
	
  // 	const djvec_t eff=effective_mass(corr[i3ptsKin],T/2,0);
  // 	grace_file_t multifit(combine("plots/3pts_bis_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
  // 	write_fit_plot(multifit,tMin[i3ptsKin],tMax[i3ptsKin],der,vector_up_to<double>(eff.size()),eff);
	
  // 	const double target=0.005;
  // 	const djack_t tMinEst=-log(target*N/Z1)/(dE[i3ptsKin]*(dE10-1));
  // 	const djack_t tMaxEst=tMinEst+4;//T/2+log(target*abs(N/Z2))/(dE[i3ptsKin]*(dE20-1));
  // 	cout<<tMinEst.ave_err()<<endl;
  // 	cout<<tMaxEst.ave_err()<<endl;
  // 	multifit.set_title(smart_print(tMinEst));
  // 	multifit.write_constant_band(tMinEst.ave(),tMaxEst.ave(),dE[i3ptsKin]);
  //     }


int main(int narg,char **arg)
{
  set_njacks(15);
  
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
  
  size_t iMs=2;
  size_t iMt=0;
  
  const index_t indMesKin({{"mom1",nMoms},{"mom2",nMoms}});
  const size_t nMesKin=indMesKin.max();
  vector<double> Pmes(nMesKin);
  djvec_t E(nMesKin),ZP(nMesKin),ZA(nMesKin);
  
  double pMax=0;
  for(size_t iMom1=0;iMom1<nMoms;iMom1++)
    for(size_t iMom2=iMom1;iMom2<nMoms;iMom2++)
      {
	const size_t i1=iMom1*nMoms+iMom2;
	const size_t i2=iMom2*nMoms+iMom1;
	
	const int tMin=12,tMax=T/2+1;
	
	const djvec_t corrPP=load2ptsPP(iMs,iMt,iMom1,iMom2);
	corrPP.ave_err().write(combine("plots/2pts_PP_corr_%d.xmg",i1));
	
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
  for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    {
      const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
      const size_t iMomt=c3pts[1],iMom0=c3pts[2];
      consider[i3ptsKin]=iMomt!=iMom0;
    }
  
  vector<double> k(ind3ptsKin.max()),khat(ind3ptsKin.max());
  for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    {
      const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
      const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
      const size_t iMes=indMesKin({iMoms,iMom0});
      const double P=2*M_PI*(moms[iMom0][2]-moms[iMoms][2])/L;
      k[i3ptsKin]=2*M_PI*(moms[iMom0][2]-moms[iMomt][2])/L;
      khat[i3ptsKin]=2*sin(k[i3ptsKin]/2);
      
      // const djack_t E=latt_en_1D(M,P);
      const double Eg=2*asinh(fabs(khat[i3ptsKin])/2);
      const double EgT=sinh(Eg)*(1-exp(-T*Eg));
      dE[i3ptsKin]=E[iMes]-Eg;
      PK[i3ptsKin]=E[iMes]*Eg-P*k[i3ptsKin];
      
      
      cout<<iMoms<<" "<<iMomt<<" "<<iMom0<<endl;
      cout<<" P: "<<P<<endl;
      cout<<" Pg: "<<k[i3ptsKin]<<endl;
      cout<<" Pghat: "<<khat[i3ptsKin]<<endl;
      cout<<" E: "<<E[iMes].ave_err()<<endl;
      cout<<" Eg: "<<Eg<<endl;
      cout<<" EgT: "<<EgT<<endl;
      cout<<" dE: "<<dE[i3ptsKin].ave_err()<<endl;
      
      for(int iVA=0;iVA<2;iVA++)
	corr[iVA][i3ptsKin]=load3pts(iVA,iMs,iMt,iMoms,iMomt,iMom0);
      
      for(int t=0;t<=T/2;t++)
	normaliz[i3ptsKin][t]=4*E[iMes]*EgT/(ZP[iMes]*exp(-t*E[iMes]-(T/2-t)*Eg));
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
	    outplot.write_vec_ave_err(eff.ave_err());
	    
	    vector<size_t> cl;
	    bool incl=false;
	    for(int t=0;t<T/2-1;t++)
	      {
		const djack_t e=eff[t]-dE[i3ptsKin];
		const double c=fabs(e.ave()/e.err());
		const int nDev=2;
		const bool newIncl=(c<nDev);
		
		if(incl and not newIncl) cl.push_back(t-1);
		if((not incl) and newIncl) cl.push_back(t);
		
		incl=newIncl;
	      }
	    
	    if(incl)
	      cl.push_back(T/2);
	    
	    if(cl.size()%2)
	      CRASH("Size %d of %d should be even",(int)i3ptsKin,(int)cl.size());
	    
	    for(size_t iCl=0;iCl<cl.size()/2;iCl++)
	      outplot.write_constant_band(cl[iCl*2]-0.5,cl[iCl*2+1]+0.5,dE[i3ptsKin]);
	    
	    cout<<"Merging "<<i3ptsKin<<endl;
	    size_t iCl=0;
	    const size_t gap=3;
	    int nMerged=0;
	    while(cl.size()/2>1 and iCl<cl.size()/2)
	      {
		auto begThis=cl.begin()+iCl;
		auto endThis=cl.begin()+iCl+1;
		auto begNext=cl.begin()+iCl+2;
		auto endNext=cl.begin()+iCl+3;
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
	    
	    for(size_t iCl=0;iCl<cl.size()/2;iCl++)
	      outplot.write_constant_band(cl[iCl*2]-0.5,cl[iCl*2+1]+0.5,djack_t(dE[i3ptsKin]-dE[i3ptsKin].err()*3));
	    
	    outplot.set_title(to_string(i3ptsKin));
	    
	    for(size_t iCl=0;iCl<cl.size()/2;iCl++)
	      if(tMax[i3ptsKin]-tMin[i3ptsKin]<cl[iCl*2+1]-cl[iCl*2])
		{
		  tMax[i3ptsKin]=cl[iCl*2+1];
		  tMin[i3ptsKin]=cl[iCl*2];
		}
	  }
      
      grace_file_t H(combine("plots/H_%s.xmg",VA_tag[iVA]));
      grace_file_t ff(combine("plots/ff_%s.xmg",VA_tag[iVA]));
      for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
	if(consider[i3ptsKin])
	  {
	    const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
	    
	    const djvec_t y=corr[iVA][i3ptsKin]*normaliz[i3ptsKin]/PK[i3ptsKin];
	    const djack_t C=constant_fit(y,tMin[i3ptsKin],tMax[i3ptsKin],combine("plots/3pts_%s_fit_%s.xmg",VA_tag[iVA],ind3ptsKin.descr(i3ptsKin).c_str()));
	    const djack_t x=2*PK[i3ptsKin]/sqr(E[0]);
	    const djack_t f=C+fP_bis/PK[i3ptsKin];
	    
	    H.write_ave_err(1/PK[i3ptsKin].ave(),C.ave_err());
	    ff.write_ave_err(x.ave(),f.ave_err());
	  }
    }
  
  return 0;
}
