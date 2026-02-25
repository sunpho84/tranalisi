
#include <tranalisi.hpp>

const size_t T{128};
const size_t L{64};
const size_t nIns{2};
const size_t nConfs{400};
const size_t nCopies{2};
const size_t nEpsIn{5};
const std::string insName[2]{"2pt","3pt"};

const size_t nCombo{7};
const double eps[2]{0.05,0.2};
const std::array<double,nCombo> coef{0,0.0,0.0,0.0,0.0,sqr(eps[0]),sqr(eps[1])};
std::vector<std::vector<std::pair<size_t,double>>> coefCombo{{{{0,+1.0}},                     //0 order
							      {{1,+1.0}},
							      {{2,+1.0}},
							      {{3,+1.0}},
							      {{4,+1.0}},
  							      {{1,+0.5},{2,+0.5},{0,-1.0}},    //1st der in eps2
							      {{3,+0.5},{4,+0.5},{0,-1.0}}}};  //same

const size_t tMin2pts=30;
const size_t tMax2pts=40;

const size_t tMin3pts=20;
const size_t tMax3pts=30;

int main()
{
  const index_t rawId{{{"iEps",nEpsIn},{"ins",nIns},{"conf",nConfs},{"copy",nCopies},{"T",T}}};
  
  std::vector<double> rawData(rawId.max(),0.0);
  
  raw_file_t fin("sanfDat.txt","r");
  for(double& d : rawData)
    fin.read(d);
  // if(not fin.feof())
  //   CRASH("not at the end of file, expecting %zu lines",rawId.max());
  
  set_njacks(50);
  
  const index_t id{{{"iCombo",nCombo},{"ins",nIns}}};
  
  std::vector<djvec_t> data(id.max(),djvec_t(T));
  for(size_t iCombo=0;iCombo<nCombo;iCombo++)
    {
      cout<<"Combo "<<iCombo<<": "<<endl;
      for(const auto& [i,c] : coefCombo[iCombo])
	cout<<"("<<i<<","<<c<<")"<<endl;
      
      for(size_t iIns=0;iIns<nIns;iIns++)
	{
	  djvec_t& dout=data[id({iCombo,iIns})];
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& w)
			 {
			   for(size_t iCopy=0;iCopy<nCopies;iCopy++)
			     for(size_t t=0;t<T;t++)
			       {
				 double s{0.0};
				 for(const auto& [i,c] : coefCombo[iCombo])
				   {
				     const size_t iLine=rawId({i,iIns,iConf,iCopy,t});
				     if(iLine>=rawId.max())
				       CRASH("unable to read i %lu, iIns %zu, iConf %zu, iCopy %zu, t %zu",i,iIns,iConf,iCopy,t);
				     s+=c*rawData[iLine];
				   }
				 dout[t][iClust]+=w*s;
			       }
			 });
	  dout.clusterize((double)nConfs/njacks);
	  dout/=nCopies;
	  if(coef[iCombo])
	    dout/=coef[iCombo];
	  
	  dout.ave_err().write("plots/"+insName[iIns]+std::to_string(iCombo)+".xmg");
	  //effective_mass(dout.subset(0,T/2),T/2,(iIns==0)?1:0).ave_err().write("plots/"+insName[iIns]+std::to_string(eps[iEps])+".xmg");
	  
	}
    }
  
  // {
  //   cout.precision(16);
  //   const djack_t E1=constant_fit(effective_mass(data[id({3,0})].symmetrized()),30,35,"plots/E1.xmg");
  //   const djack_t E2=constant_fit(effective_mass(data[id({4,0})].symmetrized()),30,35,"plots/E2.xmg");
  //   const djack_t M=constant_fit(effective_mass(data[id({0,0})].symmetrized()),30,35,"plots/M.xmg");
  //   djack_t SL;
  //   {
  //     djack_t Z2,DZ2_fr_Z2,M;
  //     two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL,data[id({0,0})].symmetrized(),data[id({5,0})].symmetrized(),T/2,30,40,"plots/t0.xmg","plots/t1.xmg");
  //     cout<<"SL from 2nd der in eps2: "<<SL.ave_err()<<endl;
  //   }
  //   cout<<"M: "<<M.ave_err()<<endl;
  //   cout<<"E1: "<<E1.ave_err()<<endl;
  //   cout<<"E2: "<<E2.ave_err()<<endl;
  //   djack_t Ee1=latt_en_1D(M,M_PI/L*0.05);
  //   cout<<"Ee1: "<<Ee1.ave_err()<<endl;
  //   const djack_t t1=(E1-M)/(sqr(0.05));
  //   const djack_t t2=(E2-M)/(sqr(0.05));
  //   cout<<t1.ave_err()<<endl;
  //   cout<<t2.ave_err()<<endl;
  //   cout<<smart_print((Ee1-M)/(sqr(0.05)))<<endl;
  // }
  
  const djvec_t c2pt=data[id({0,0})].symmetrized();
  djack_t Z2,M;
  djvec_t DZ2_fr_Z2(2),SL(2);
  for(size_t iEst=0;iEst<2;iEst++)
    {
      const djvec_t c2ptIns=data[id({iEst+nEpsIn,0})].symmetrized();
      two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2[iEst],SL[iEst],c2pt,c2ptIns,T/2,tMin2pts,tMax2pts,"plots/c2pt0.xmg","plots/c2pt"+std::to_string(iEst+nEpsIn)+".xmg");
      effective_slope(c2ptIns/c2pt,effective_mass(c2ptIns),T/2).ave_err().write("plots/c2ptSlope"+std::to_string(iEst+nEpsIn)+".xmg");
    }
  cout<<"M: "<<smart_print(M)<<endl;
  cout<<"SL0: "<<smart_print(SL[0])<<endl;
  cout<<"SL1: "<<smart_print(SL[1])<<endl;
  
  const djack_t expSL=-sqr(M_PI/L)/(2*M);
  cout<<"expSL: "<<smart_print(expSL)<<endl;
  
  const djack_t a=M/2.9839;
  cout<<"a: "<<a.ave_err()<<endl;
  
  const size_t tw=35;
  djvec_t c3pt=data[id({0,1})].subset(tw,T/2);
  
  djvec_t czv(T/2+1-tw);
  for(size_t it=0;it<czv.size();it++)
    czv[it]=two_pts_corr_fun(Z2,M,T/2.0,tw+it,0)/c3pt[it];
   const djack_t zv=constant_fit(czv,tMin3pts,tMax3pts,"plots/zv.xmg");
  
  djvec_t dff[2];
  djack_t D[2];
  for(size_t iEps=0;iEps<2;iEps++)
    {
      djvec_t c3ptIns=data[id({iEps+nEpsIn,1})].subset(tw,T/2);
      const djvec_t r3pt=c3ptIns/c3pt;
      r3pt.ave_err().write("plots/r3pt"+std::to_string(iEps)+".xmg");
      
      dff[iEps]=r3pt;
      for(size_t t=0;t<dff[iEps].size();t++)
	{
	  dff[iEps][t]-=DZ2_fr_Z2[iEps]/2+SL[iEps]*(1/(2*M)+tw);
	  dff[iEps][t]*=-6.0*sqr(L)/sqr(M_PI)*sqr(a);
	}
      D[iEps]=constant_fit(dff[iEps],20,25,"plots/dff"+std::to_string(iEps)+".xmg");
      cout<<"Eps "<<eps[iEps]<<" : "<<smart_print(D[iEps])<<endl;
    }
  
  djvec_t diff=dff[1]-dff[0];
  diff.ave_err().write("plots/diff.xmg");
  
  // for(size_t t=0;t<c3pt.size();t++)
  //   c3pt[t]/=exp(-M*(tw+t))*Z2/M;
  // c3pt.ave_err().write("plots/c3pt.xmg");
  
  /////////////////////////////////////////////////////////////////
  
  {
    djack_t Z2[nEpsIn],E[nEpsIn];
    djack_t F[nEpsIn];
    for(size_t iEps=0;iEps<nEpsIn;iEps++)
      {
	two_pts_fit(Z2[iEps],E[iEps],data[id({iEps,0})].symmetrized(),T/2,tMin2pts,tMax2pts,"plots/m"+to_string(iEps)+".xmg");
	
	djvec_t cF(T/2-tw+1);
	for(size_t t=0;t<=T/2-tw;t++)
	  cF[t]=4*E[0]*E[iEps]*data[id({iEps,1})][tw+t]/sqrt(Z2[0]*Z2[iEps])/(E[0]+E[iEps])
	    /exp(-E[iEps]*tw)/exp(-E[0]*t);
	
	F[iEps]=constant_fit(cF,tMin3pts,tMax3pts,"plots/F"+to_string(iEps)+".xmg");
      }
    
    djack_t dz2_fr_z2[2];
    djack_t sl[2];
    for(size_t iEps=0;iEps<2;iEps++)
      {
	dz2_fr_z2[iEps]=((Z2[2*iEps+1]+Z2[2*iEps+2])/(2*Z2[0])-1)/sqr(eps[iEps]);
	cout<<"dz/z2: "<<smart_print(dz2_fr_z2[iEps])<<" vs "<<smart_print(DZ2_fr_Z2[iEps])<<endl;
	
	sl[iEps]=((E[2*iEps+1]+E[2*iEps+2])/2-E[0])/sqr(eps[iEps]);
	cout<<"sl: "<<smart_print(sl[iEps])<<" vs "<<smart_print(SL[iEps])<<endl;
      }
    
    djack_t r2[2];
    for(size_t iEps=0;iEps<2;iEps++)
      {
	//cout<<" "<<smart_print(F[2*iEps+1])<<" "<<smart_print(F[2*iEps+2])<<" "<<smart_print(F[0])<<endl;
	const djack_t G=(F[2*iEps+1]+F[2*iEps+2])/(2*F[0])-1;
	//cout<<"G: "<<smart_print(G)<<endl;
	
	const djack_t Q2=sqr(E[0]-E[2*iEps+1])-sqr(M_PI*eps[iEps]/L);
	//cout<<"Q2: "<<smart_print(Q2)<<endl;
	
	r2[iEps]=6*G/Q2*sqr(a);
	cout<<"r2 "<<eps[iEps]<<": "<<smart_print(r2[iEps])<<endl;
      }
    
    cout<<smart_print(r2[1]-r2[0])<<endl;
  }
  
  return 0;
}
