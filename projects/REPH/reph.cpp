#include <tranalisi.hpp>

size_t L,T,spatVol;
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
  


djvec_t load2pts(const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt)
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

djvec_t load3pts(const string VA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0)
{
  const double s[2][2]={-1,-1,+1,-1};
  
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
  
  const size_t iReIm=0;
  for(size_t iPol=0;iPol<2;iPol++)
    for(size_t iGamma=1;iGamma<=2;iGamma++)
      {
	const size_t i=ind({iMs,iMt,iMoms,iMomt,iMom0,iPol,iGamma,iReIm});
	
	corr+=read_djvec("jacks/o" + VA + "muGPo-gs",T,i)*s[iPol][iGamma-1];
      }
  corr/=2*sqrt(2);
  
  return corr.symmetrized()/(L*L*L);
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
  
  grace_file_t fit2ptsPlot("plots/fit2pts.xmg");
  grace_file_t Z2Plot("plots/Z2.xmg");
  
  size_t iMs=0;
  size_t iMt=0;
  
  const index_t indMesKin({{"mom1",nMoms},{"mom2",nMoms}});
  const size_t nMesKin=indMesKin.max();
  vector<double> Pmes(nMesKin);
  djvec_t E(nMesKin),Z(nMesKin);
  
  double pMax=0;
  for(size_t iMom1=0;iMom1<nMoms;iMom1++)
    for(size_t iMom2=iMom1;iMom2<nMoms;iMom2++)
      {
	const size_t i1=iMom1*nMoms+iMom2;
	const size_t i2=iMom2*nMoms+iMom1;
	
	const int tmin=12,tmax=T/2+1;
	
	const djvec_t corr=load2pts(iMs,iMt,iMom1,iMom2);
	
	djack_t Z2,EFit;
	two_pts_fit(Z2,EFit,corr,T/2,tmin,tmax);
	Z[i1]=Z[i2]=sqrt(Z2);
	E[i1]=E[i2]=EFit;
	
	fit2ptsPlot.write_vec_ave_err(effective_mass(corr).ave_err());
	fit2ptsPlot.write_constant_band(tmin,tmax,EFit);
	
	double P=2*3.14159*(moms[iMom2][2]-moms[iMom1][2])/L;
	Pmes[i1]=+P;
	Pmes[i2]=-P;
	
	pMax=std::max(pMax,fabs(P));
	
	Z2Plot.write_ave_err(EFit.ave(),Z[i1].ave_err());
      }
  
  grace_file_t dispRel("plots/dispRel.xmg");
  dispRel.write_vec_ave_err(Pmes,E.ave_err());
  
  const djack_t &M=E[0];
  
  dispRel.write_polygon([M](double x){return latt_en_1D(M,x);},0,pMax,grace::GREEN);
  dispRel.write_polygon([M](double x){return cont_en_1D(M,x);},0,pMax,grace::VIOLET);
  
  /////////////////////////////////////////////////////////////////
  
  const index_t ind3ptsKin({{"iMoms",nMoms},{"iMomt",nMoms},{"iMom0",nMoms}});
  
  vector<djvec_t> corr(ind3ptsKin.max(),djvec_t{T/2+1});
  vector<djvec_t> normaliz(ind3ptsKin.max(),djvec_t{T/2+1});
  djvec_t dE(ind3ptsKin.max());
  
  vector<bool> consider(ind3ptsKin.max());
  for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    {
      const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
      const size_t iMomt=c3pts[1],iMom0=c3pts[2];
      consider[i3ptsKin]=iMomt!=iMom0;
    }

  for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    {
      const vector<size_t> c3pts=ind3ptsKin(i3ptsKin);
      const size_t iMoms=c3pts[0],iMomt=c3pts[1],iMom0=c3pts[2];
      const size_t iMes=indMesKin({iMoms,iMom0});
      const double P=Pmes[iMes];
      const double Pg=2*M_PI*(moms[iMomt][2]-moms[iMom0][2])/L;
      
      const djack_t E=latt_en_1D(M,P);
      const double Eg=fabs(Pg);
      dE[i3ptsKin]=E-Eg;
      
      cout<<iMoms<<" "<<iMomt<<" "<<iMom0<<endl;
      cout<<" P: "<<P<<endl;
      cout<<" Pg: "<<Pg<<endl;
      cout<<" E: "<<E.ave_err()<<endl;
      cout<<" Eg: "<<Eg<<endl;
      cout<<" dE: "<<dE[i3ptsKin].ave_err()<<endl;
      
      corr[i3ptsKin]=load3pts("A",iMs,iMt,iMoms,iMomt,iMom0);
      
      for(size_t t=0;t<=T/2;t++)
	normaliz[i3ptsKin][t]=4*E*Eg/(Z[0]*exp(-dE[i3ptsKin]*t)*exp(-(T*Eg/2)));
    }
  
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
  
  vector<vector<double>> Ch2(ind3ptsKin.max());
  vector<array<double,3>> p3pts(ind3ptsKin.max());
  for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    {
      const vector<size_t> c=ind3ptsKin(i3ptsKin);
      for(size_t i=0;i<3;i++)
	p3pts[i3ptsKin][i]=moms[c[i]][2];
      if(consider[i3ptsKin])
      	Ch2[i3ptsKin]=compatibilityWith(effective_mass(corr[i3ptsKin],T/2,0),dE[i3ptsKin]);
    }
  
  size_t tMinBest=0,tMaxBest=T/2-1;
  double TotCh2Best=1e300;
  int i=0,j=0;
  for(size_t tMin=4;tMin<T/2-4;tMin++)
    for(size_t tMax=tMin+4;tMax<T/2-4;tMax++)
      for(size_t tMin1=4;tMin1<T/2-4;tMin1++)
	for(size_t tMax1=tMin1+4;tMax1<T/2-4;tMax1++)
	  for(size_t tMin2=4;tMin2<T/2-4;tMin2++)
	    for(size_t tMax2=tMin2+4;tMax2<T/2-4;tMax2++)
	      for(size_t tMin3=4;tMin3<T/2-4;tMin3++)
		for(size_t tMax3=tMin3+4;tMax3<T/2-4;tMax3++)
		  {
		    double a=tMin;
		    double b=tMax;
		    
		    double a1=((double)tMin1-tMin)/moms.back()[2]+1e-10;
		    double a2=((double)tMin2-tMin)/moms.back()[2]+1e-10;
		    double a3=((double)tMin3-tMin)/moms.back()[2]+1e-10;
		    double b1=((double)tMax1-tMax)/moms.back()[2]+1e-10;
		    double b2=((double)tMax2-tMax)/moms.back()[2]+1e-10;
		    double b3=((double)tMax3-tMax)/moms.back()[2]+1e-10;
		    
		    // cout<<"tMin: "<<tMin<<" "<<tMin1<<" "<<tMin2<<" "<<tMin3<<endl;
		    // cout<<"tMax: "<<tMax<<" "<<tMax1<<" "<<tMax2<<" "<<tMax3<<endl;
		    
		    // cout<<" a: "<<a1<<" "<<a2<<" "<<a3<<endl;
		    // cout<<" b: "<<b1<<" "<<b2<<" "<<b3<<endl;
		    
		    i++;
		    double TotCh2=0;
		    int n=0;
		    for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
		      if(consider[i3ptsKin])
			{
			  const array<double,3> &p=p3pts[i3ptsKin];
			  int tMi=round(a+a1*p[0]+a2*p[1]+a3*p[2]);
			  int tMa=round(b+b1*p[0]+b2*p[1]+b3*p[2]);

			  if(tMa<tMi) CRASH("%d %d  %lg %lg %lg",tMi,tMa,p[0],p[1],p[2]);
			  
			  //cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<tMi<<" "<<tMa<<endl;
			  
			  TotCh2+=Ch2[i3ptsKin][tMa+T/2*tMi];//*(tMa-tMi+1);
			  n+=tMa-tMi+1;
			}
		    TotCh2/=n;
		    
		    if(TotCh2<TotCh2Best)
		      {
			TotCh2Best=TotCh2;
			tMinBest=tMin;
			tMaxBest=tMax;
			
			cout<<"New best: "<<tMinBest<<" "<<tMaxBest<<" "<<tMin1<<" "<<tMax1<<" "<<tMin2<<" "<<tMax2<<" "<<tMin3<<" "<<tMax3<<" "<<TotCh2<<" "<<n<<endl;
		      }
		    
		    if(i>j*19448100)
		      {
			j++;
			cout<<j<<endl;
		      }
		  }
  
  cout<<"Best range: "<<tMinBest<<" "<<tMaxBest<<" "<<TotCh2Best<<endl;
  
  for(size_t i3ptsKin=0;i3ptsKin<ind3ptsKin.max();i3ptsKin++)
    if(consider[i3ptsKin])
      {
	grace_file_t outplot(combine("plots/3pts_range_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
	outplot.write_vec_ave_err(effective_mass(corr[i3ptsKin],T/2,0).ave_err());
	outplot.write_constant_band(tMinBest,tMaxBest,dE[i3ptsKin]);
      }
  
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
  //   {
  //     grace_file_t outplot(combine("plots/3pts_%s.xmg",ind3ptsKin.descr(i3ptsKin).c_str()));
  //     outplot.write_vec_ave_err(effective_mass(corr[i3ptsKin],T/2,0).ave_err());
  //     outplot.write_constant_band(0,T/2,dE[i3ptsKin]);
      
  //     // djvec_t y=corr[i3ptsKin]*normaliz[i3ptsKin];
      
  //     // outplot.write_vec_ave_err(y.ave_err());
  //   }

  
  return 0;
}
