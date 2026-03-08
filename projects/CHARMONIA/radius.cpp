#include <tranalisi.hpp>

const size_t T{128};
const size_t tw{35};
const size_t L{64};
const size_t nIns{2};
const size_t nConfs{400};
const size_t nCopies{2};
const size_t nEpsIn{5};
const std::string insName[2]{"2pt","3pt"};

const double eps[2]{0.05,0.2};

const size_t tMin2pts=30;
const size_t tMax2pts=40;

const size_t tMin2ptsIns=35;
const size_t tMax2ptsIns=45;

const size_t tMin3pts=22;
const size_t tMax3pts=37;

int main()
{
  const index_t rawId{{{"iEps",nEpsIn},{"ins",nIns},{"conf",nConfs},{"copy",nCopies},{"T",T}}};
  
  std::vector<double> rawData(rawId.max(),0.0);
  
  raw_file_t fin("sanfDat.txt","r");
  for(double& d : rawData)
    fin.read(d);
  // if(not fin.feof())
  //   CRASH("not at the end of file, expecting %zu lines",rawId.max());
  
  useEnlargedJacks=true;
  useFastAveErr=false;
  set_njacks(50);
  
  const index_t id{{{"iEpsIn",nEpsIn},{"ins",nIns}}};
  
  std::vector<djvec_t> cData(id.max(),djvec_t(T));
  for(size_t iEpsIn=0;iEpsIn<nEpsIn;iEpsIn++)
      for(size_t iIns=0;iIns<nIns;iIns++)
	{
	  djvec_t& dout=cData[id({iEpsIn,iIns})];
	  dout=0.0;
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& w)
			 {
			   for(size_t t=0;t<T;t++)
			     {
			       double s=0.0;
			       for(size_t iCopy=0;iCopy<nCopies;iCopy++)
				 {
				   const size_t iLine=rawId({iEpsIn,iIns,iConf,iCopy,t});
				   if(iLine>=rawId.max())
				     CRASH("unable to read iEpsIn %lu, iIns %zu, iConf %zu, iCopy %zu, t %zu",iEpsIn,iIns,iConf,iCopy,t);
				   s+=w*rawData[iLine];
				 }
			       dout[t][iClust]+=s;
			     }
			 });
	  dout.clusterize((double)nConfs/njacks);
	  dout/=nCopies;
	  
	  dout.ave_err().write("plots/"+insName[iIns]+std::to_string(iEpsIn)+".xmg");
	  //effective_mass(dout.subset(0,T/2),T/2,(iIns==0)?1:0).ave_err().write("plots/"+insName[iIns]+std::to_string(eps[iEps])+".xmg");
	  
	}
  
  const auto c2ptFun=
    [&cData,
     &id](const size_t iEps) -> djvec_t
    {
      return cData[id({iEps,0})].symmetrized();
    };
  
  [gg=ofstream("/tmp/gg"),&rawData,&rawId]() mutable
  {
    const size_t t=12;
    const size_t binSize=nConfs/njacks;
    vector<double> j(njacks+1,0);
    for(size_t iConf=0;iConf<nConfs;iConf++)
      {
	const size_t iJack=iConf/binSize;
	for(size_t iCopy=0;iCopy<2;iCopy++)
	  {
	    j[iJack]+=rawData[rawId({1,0,iConf,iCopy,t})]/4;
	    j[iJack]+=rawData[rawId({1,0,iConf,iCopy,T-t})]/4;
	  }
      }
    
    for(size_t iJack=0;iJack<njacks;iJack++)
      j[njacks]+=j[iJack];
    for(size_t iJack=0;iJack<njacks;iJack++)
      j[iJack]=(j[njacks]-j[iJack])/(nConfs-binSize);
    j[njacks]/=nConfs;
    
    gg.precision(16);
    for(size_t iJack=0;iJack<=njacks;iJack++)
      gg<<j[iJack]<<endl;
  }();
  
  for(size_t iEps=0;iEps<nEpsIn;iEps++)
    c2ptFun(iEps).ave_err().write("plots/c2pts"+std::to_string(iEps)+".xmg");
  
  const auto c3ptFun=
    [&cData,
     &id](const size_t& iEps)
    {
      return cData[id({iEps,1})].subset(tw,T/2);
    };
  
  auto dumpJack=
    [](const string& path,
       const djvec_t& v)
  {
    ofstream f(path);
    
    f.precision(16);
    for(size_t t=0;t<v.size();t++)
      for(size_t ijack=0;ijack<=njacks;ijack++)
	f<<t<<" "<<ijack<<" "<<v[t][ijack]<<endl;
  };
  
  dumpJack("/tmp/hh",c2ptFun(1));
  
  for(size_t iEps=0;iEps<nEpsIn;iEps++)
    c3ptFun(iEps).ave_err().write("plots/c3pts"+std::to_string(iEps)+".xmg");
  
  /////////////////////////////////////////////////////////////////
  
  const djvec_t c2pt=c2ptFun(0);
  const djvec_t eff_mass=effective_mass(c2pt,T/2);
  dumpJack("/tmp/m.txt",eff_mass);
  const djvec_t eff_sq_coupling=effective_squared_coupling(c2pt,eff_mass,T/2);
  dumpJack("/tmp/z2.txt",eff_sq_coupling);
  
  const djack_t myM=constant_fit(eff_mass,tMin2pts,tMax2pts,"plots/c2ptEffmass.xmg");
  
  /////////////////////////////////////////////////////////////////
  
  /// Finite eps analysis
  cout<<"---- Finite eps analysis ---"<<endl;
  
  djack_t Z2[nEpsIn],E[nEpsIn];
  djack_t F[nEpsIn];
  djvec_t ss[nEpsIn];
  for(size_t iEps=0;iEps<nEpsIn;iEps++)
    {
      const djvec_t c2=c2ptFun(iEps);
      two_pts_fit(Z2[iEps],E[iEps],c2,T/2,tMin2pts,tMax2pts,"plots/m"+to_string(iEps)+".xmg");
      const djvec_t eff_massE=effective_mass(c2,T/2);
      const djvec_t eff_sq_couplingE=effective_squared_coupling(c2,eff_massE,T/2);
      
      djvec_t cF(T/2-tw);
      ss[iEps].resize(T/2-tw);
      const djvec_t c3=c3ptFun(iEps);
      for(size_t t=0;t<T/2-tw;t++)
	{
	  const djack_t m=E[0];
	  const djack_t e=E[iEps];
	  const djack_t z20=Z2[0];
	  const djack_t z2=Z2[iEps];
	  // const djack_t m=eff_mass[t+tw];
	  // const djack_t e=eff_massE[t+tw];
	  // const djack_t z20=eff_sq_coupling[t+tw];
	  // const djack_t z2=eff_sq_couplingE[t+tw];
	  
	  ss[iEps][t]=e*c3[t]/sqrt(z2)/(m+e)
	    /exp(-e*tw)/exp(-m*t);
	  
	  cF[t]=4*m*e*c3[t]/sqrt(z20*z2)/(m+e)
	    /exp(-e*tw)/exp(-m*t);
	}
      F[iEps]=constant_fit(cF,tMin3pts,tMax3pts,"plots/F"+to_string(iEps)+".xmg");
      ss[iEps].ave_err().write("plots/ss"+to_string(iEps)+".xmg");
    }
  
  djvec_t Q2(2);
  for(size_t iEst=0;iEst<2;iEst++)
    Q2[iEst]=sqr(E[0]-E[2*iEst+1])-sqr(M_PI*eps[iEst]/L);
  
  for(size_t iEst=0;iEst<2;iEst++)
    {
      const djvec_t s=(ss[2*iEst+1]+ss[2*iEst+2])/(2*ss[0])-1;
      (s/Q2[iEst]).ave_err().write("plots/s"+to_string(iEst)+".xmg");
    }
  
  const djack_t& M=E[0];
  cout<<"M: "<<smart_print(M)<<endl;
  
  const djack_t a=M/2.9839;
  cout<<"a: "<<smart_print(a)<<" GeV^-1"<<endl;
  
  [ii=ofstream("/tmp/ii"),&F]() mutable
  {
    ii.precision(16);
    ii<<" # ieps ijack f[+eps] f[-eps] f0[]"<<endl;
    for(size_t iEst=0;iEst<2;iEst++)
      for(size_t ijack=0;ijack<=njacks;ijack++)
	ii<<iEst<<" "<<ijack<<" "<<F[2*iEst+2][ijack]<<" "<<F[2*iEst+1][ijack]<<" "<<F[0][ijack]<<endl;
  }();
  
  djack_t r2[2];
  // djvec_t Q2G(2);
  for(size_t iEst=0;iEst<2;iEst++)
    {
      //cout<<" "<<smart_print(F[2*iEps+1])<<" "<<smart_print(F[2*iEps+2])<<" "<<smart_print(F[0])<<endl;
      const djack_t G=(F[2*iEst+1]+F[2*iEst+2])/(2*F[0])-1;
      //cout<<"G: "<<smart_print(G)<<endl;
      
      const djack_t EG=sqrt(sqr(E[0])+sqr(eps[iEst]*M_PI/L));
      // Q2G[iEst]=-2.0*E[0]*(EG-E[0]);
      // cout<<"Q2 unc: "<<smart_print((Q2[iEst]-Q2G[iEst])*1e12)<<endl;
      
      r2[iEst]=6*G/Q2[iEst]*sqr(a);
      cout<<"Eps "<<eps[iEst]<<": "<<smart_print(r2[iEst])<<endl;
    }
  
  cout<<"diff: "<<smart_print(r2[1]-r2[0])<<endl;
  
  for(size_t iEst=0;iEst<2;iEst++)
    {
      const djack_t d=Q2[iEst]/(sqr(M_PI/L*eps[iEst]));
      cout<<"Q2/p2: "<<d<<endl;
    }
  
  /////////////////////////////////////////////////////////////////
  
  // Slope analysis
  cout<<"--- Slope analysis ---"<<endl;
  
  const auto getCInsFun=
    [](const auto& cFun)
    {
      return
	[&cFun](const size_t& iEst)->djvec_t
	{
	  return (cFun(2*iEst+1)+cFun(2*iEst+2)-2*cFun(0))/(2*sqr(eps[iEst]));
	};
    };
  
  const auto c2ptInsFun=
    getCInsFun(c2ptFun);
  
  for(size_t iEst=0;iEst<2;iEst++)
    c2ptInsFun(iEst).ave_err().write("plots/c2ptsIns"+std::to_string(iEst)+".xmg");
  
  djvec_t DZ2_fr_Z2(2),SL(2);
  djvec_t c2ptInsRat[2];
  djvec_t eff_slope[2];
  djvec_t eff_slope_offset[2];
  
  for(size_t iEst=0;iEst<2;iEst++)
    {
      c2ptInsRat[iEst]=c2ptInsFun(iEst)/c2pt;
      eff_slope[iEst]=effective_slope(c2ptInsRat[iEst],eff_mass,T/2);
      SL[iEst]=constant_fit(eff_slope[iEst],tMin2ptsIns,tMax2ptsIns,"plots/c2ptSlope"+std::to_string(iEst)+".xmg");
      eff_slope_offset[iEst]=effective_squared_coupling_rel_corr(c2ptInsRat[iEst],eff_mass,eff_slope[iEst],T/2);
      DZ2_fr_Z2[iEst]=constant_fit(eff_slope_offset[iEst],tMin2ptsIns,tMax2ptsIns,"plots/c2ptSlopeOffset"+std::to_string(iEst)+".xmg");
    }
  cout<<"SL0: "<<smart_print(SL[0])<<endl;
  cout<<"SL1: "<<smart_print(SL[1])<<endl;
  
  const djack_t expSL=-sqr(M_PI/L)/(2*M);
  cout<<"expSL: "<<smart_print(expSL)<<endl;
  
  const auto c3ptInsFun=getCInsFun(c3ptFun);
  
  for(size_t iEst=0;iEst<2;iEst++)
    c3ptInsFun(iEst).ave_err().write("plots/c3ptsIns"+std::to_string(iEst)+".xmg");
  
  const djvec_t c3pt=c3ptFun(0);
  
  djvec_t dff[2];
  djack_t R2[2];
  for(size_t iEst=0;iEst<2;iEst++)
    {
      djvec_t c3ptIns=c3ptInsFun(iEst);
      const djvec_t r3pt=c3ptIns/c3pt;
      r3pt.ave_err().write("plots/r3pt"+std::to_string(iEst)+".xmg");
      
      
      for(size_t t=0;t<T/2-tw;t++)
	{
	  dff[iEst].resize(T/2-tw);
	  const djack_t dz2_fr_z2=DZ2_fr_Z2[iEst];
	  const djack_t sl=SL[iEst];
	  const djack_t m=M;
	  // const djack_t dz2_fr_z2=eff_slope_offset[iEst][t+tw];
	  // const djack_t sl=eff_slope[iEst][t+tw];
	  // const djack_t m=eff_mass[t+tw];
	  dff[iEst][t]=r3pt[t];
	  dff[iEst][t]-=dz2_fr_z2/2+sl*(1/(2*m)+tw);
	  dff[iEst][t]*=-6.0*sqr(L)/sqr(M_PI)*sqr(a);
	}
      R2[iEst]=constant_fit(dff[iEst],tMin3pts,tMax3pts,"plots/dff"+std::to_string(iEst)+".xmg");
      cout<<"Eps "<<eps[iEst]<<" : "<<smart_print(R2[iEst])<<endl;
    }
  
  const djvec_t diff=dff[1]-dff[0];
  diff.ave_err().write("plots/diff.xmg");
  
  cout<<"diff: "<<smart_print(R2[1]-R2[0])<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  // Check on dZ2/Z2 ans sl
  cout<<"--- Comparison of finite eps and slope ---"<<endl;
  
  djack_t dz2_fr_z2[2];
  djack_t sl[2];
  for(size_t iEps=0;iEps<2;iEps++)
    {
      dz2_fr_z2[iEps]=((Z2[2*iEps+1]+Z2[2*iEps+2])/(2*Z2[0])-1)/sqr(eps[iEps]);
      cout<<"dz/z2: "<<smart_print(dz2_fr_z2[iEps])<<" vs "<<smart_print(DZ2_fr_Z2[iEps])<<endl;
    }
  
  for(size_t iEps=0;iEps<2;iEps++)
    {
      sl[iEps]=((E[2*iEps+1]+E[2*iEps+2])/2-E[0])/sqr(eps[iEps]);
      cout<<"sl: "<<smart_print(sl[iEps])<<" vs "<<smart_print(SL[iEps])<<endl;
    }
  
  cout<<"m: "<<smart_print(M)<<" vs "<<smart_print(myM)<<endl;
  
  return 0;
}
