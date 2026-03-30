#include <tranalisi.hpp>

const string dataPath="data";
vector<string> confs;
size_t nConfs;

const size_t T=128;
const size_t L=64;
const size_t tMaxBox=26;

struct InterpDef
{
  std::string rep;
  
  std::string id;
  
  std::vector<string> rap;
};

const std::vector<InterpDef> interpDef{
  {"001","1",{"Pz","Mz"}},
  {"011","2",{"0P11","0M11"}},
  {"111","3",{"P111","M111"}},
  {"002","4",{"Pz2","Mz2"}},
  {"012","5",{"0P12","0M12"}}
};

inline void setConfs()
{
  const std::string confsListPath="confsList.dat";
  
  if(not file_exists(confsListPath))
    {
      cout<<"Searching for confs"<<endl;
      
      map<pair<int,int>,string> tmpConfs;
      for(const auto& entry : filesystem::directory_iterator(dataPath))
	if(entry.is_directory())
	  {
	    const string& conf=filesystem::path(entry.path()).filename();
	    const size_t iStream=atoi(conf.substr(6,1).c_str());
	    const size_t iConf=atoi(conf.substr(0,4).c_str());
	    const int par=(iStream%2)*2-1;
	    const int offs=iStream/2*2;
	    tmpConfs[{offs,iConf*par}]=conf;
	  }
      
      for(const auto& [dum,conf] : tmpConfs)
	confs.emplace_back(conf);
      
      raw_file_t(confsListPath,"w").bin_write(confs);
    }
  else
    {
      cout<<"Readinf conf list from "<<confsListPath<<endl;
      raw_file_t(confsListPath,"r").bin_read(confs);
    }
  
  nConfs=confs.size();
  cout<<"NConfs: "<<nConfs<<endl;
}

inline map<string,vector<vector<vector<double>>>> getRaw(const char* cachedFilePath,
							 const char* rawFileNameTemplate,
							 const std::vector<const char*>& suffixList,
							 const size_t& tMax,
							 size_t& nHits)
{
  map<string,vector<vector<vector<double>>>> rawData;
  
  if(file_exists(cachedFilePath))
    {
      cout<<"Reading cached data from file "<<cachedFilePath<<endl;
      raw_file_t(cachedFilePath,"r").bin_read(rawData);
    }
  else
    {
      cout<<"Reading all data from "<<rawFileNameTemplate<<" raw files"<<endl;
      for(const char* const& suffix : suffixList)
	for(const filesystem::path conf : confs)
	  {
	    raw_file_t file(dataPath+"/"+conf.string()+"/"+combine(rawFileNameTemplate,suffix),"r");
	    char line[1024];
	    auto readLine=[&file,
			   &line]()
	    {
	      bool r=false;
	      while((not r) and (not file.feof()))
		{
		  file.get_line(line);
		  for(char* c=line;*c!='\0' and not r;c++)
		    r|=(*c!=' ' and *c!='\0');
		};
	      
	      return r;
	    };
	    
	    map<string,vector<vector<double>>> confData;
	    
	    string baseTag;
	    while(readLine())
	      {
		char a[100]{},b[100]{},c[100];
		if(sscanf(line," # Contraction of %s ^ \\dag and %s",a,b)==2)
		  baseTag=(string)a+"__"+b;
		else
		  {
		    if(sscanf(line," # %s",c)!=1)
		      CRASH("Unable to get the corr name on line: %s",line);
		    else
		      {
			const string tag=baseTag+"__"+c;
			
			if(c==std::string("P5P5"))
			  {
			    vector<double>& data=confData[tag].emplace_back(2*tMax);
			    for(size_t t=0;t<T;t++)
			      {
				if(not readLine())
				  CRASH("Unable to read time %zu for contr %s %s",t,a,b);
				
				double r,i;
				if(sscanf(line,"%lg %lg",&r,&i)!=2)
				  CRASH("Unable to convert %s to two doubles",line);
				
				if(t<tMax)
				  {
				    data[t]=r;
				    data[t+tMax]=i;
				  }
			      }
			  }
			else
			  for(size_t t=0;t<T;t++)
			    if(not readLine())
			      CRASH("Unable to read time %zu for contr %s %s",t,a,b);
		      }
		  }
	      }
	    
	    size_t iC=0;
	    for(const auto& [tag,in] : confData)
	      {
		auto& outs=
		  rawData[tag];
		
		if(outs.empty() or in.size()==outs.front().size())
		  outs.push_back(in);
		else
		  CRASH("Please remove conf %s",confs[iC].c_str());
		
		iC++;
	      }
	  }
    }
  
  raw_file_t(cachedFilePath,"w").bin_write(rawData);
  
  if(size_t tmpNConfs=rawData.begin()->second.size();nConfs!=tmpNConfs)
    // {
    CRASH("nConfs has impossibly been changed from %zu to %zu",tmpNConfs,nConfs);
  //nConfs=tmpNConfs;
  //cout<<"Adjusted nConfs to: "<<nConfs<<endl;
  // }
  
  nHits=rawData.begin()->second.front().size();
  //cout<<"Setting nHits to: "<<nHits<<endl;
  
  return rawData;
}

void box()
{
  size_t nHits;
  
  auto getRawBox=
    [rawData=
     getRaw("rawBox.dat",
	    "mes_contr_box_src_snk%s_T25",
	    {"1","2","3","4","5"},
	    tMaxBox,
	    nHits),
     &nHits](const std::string& a,
	     const std::string& b)
    {
      const index_t idx({{"hit",nHits},{"tMax",tMaxBox},{"conf",nConfs}});
      
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("%s__%s,__P5P5",a.c_str(),b.c_str());
      // cout<<"Searching for "<<what<<endl;
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	CRASH("Unable to find %s",what.c_str());
      
      const auto& v=_v->second;
      
      for(size_t iHit=0;iHit<nHits;iHit++)
	{
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<tMaxBox;t++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		res[idx({iHit,t,iConf})]+=c;
	      }
	}
      
      return [res,
	      idx](const size_t& iHit,
		   const size_t& t,
		   const size_t& iConf)
      {
	return res[idx({iHit,t,iConf})];
      };
    };
  
  auto getBox=
    [&](const string& mso1,
	const string& mso2,
	const string& msi1,
	const string& msi2)
    {
      djvec_t res(tMaxBox);
      
      const auto d=
	getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
      for(size_t t=0;t<tMaxBox;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
		       {
			 double o=0;
			 for(size_t iHit=0;iHit<nHits;iHit++)
			   o+=d(iHit,t,iConf).real();
			 o/=nHits;
			 res[t][iClust]+=weight*o;
		       });
	}
      res.clusterize(((double)nConfs/njacks));
      
      return res;
    };
  
  for(const InterpDef& bSo : interpDef)
    for(const InterpDef& bSi : interpDef)
      {
	const std::string& repSo{bSo.rep};
	const std::string& repSi{bSi.rep};
	const std::string& so1{bSo.rap[0]};
	const std::string& so2{bSo.rap[1]};
	const std::string& si{bSi.id};
	const djvec_t A=getBox(so1,so2,"M"+si,"P"+si);
	const djvec_t B=getBox(so2,so1,"M"+si,"P"+si);
	
	// effective_mass(g,T/2).ave_err().write("plots/A_"+repSi+"_"+repSo+".xmg");
	// effective_mass(h,T/2).ave_err().write("plots/B_"+repSi+"_"+repSo+".xmg");
	A.ave_err().write("plots/A_"+repSi+"_"+repSo+".xmg");
	B.ave_err().write("plots/B_"+repSi+"_"+repSo+".xmg");
	(A-B).ave_err().write("plots/C_"+repSi+"_"+repSo+".xmg");
      }
      
  // const djvec_t mzpzm1p1=getBox("Mz","Pz","M1","P1");
  // const djvec_t pzmzm1p1=getBox("Pz","Mz","M1","P1");
  // const djvec_t mz2pz2m1p1=getBox("Mz2","Pz2","M1","P1");
  // const djvec_t pz2mz2m1p1=getBox("Pz2","Mz2","M1","P1");
  // const djvec_t mzpzm2p2=getBox("Mz","Pz","M4","P4");
  // const djvec_t pzmzm2p2=getBox("Pz","Mz","M4","P4");
  // const djvec_t mxpxm1p1=getBox("Mx","Px","M1","P1");
  // const djvec_t pxmxm1p1=getBox("Px","Mx","M1","P1");
  
  // mzpzm1p1.ave_err().write("/tmp/mzpzm1p1.xmg");
  // pzmzm1p1.ave_err().write("/tmp/pzmzm1p1.xmg");
  // mz2pz2m1p1.ave_err().write("/tmp/mz2pz2m1p1.xmg");
  // pz2mz2m1p1.ave_err().write("/tmp/pz2mz2m1p1.xmg");
  // mzpzm2p2.ave_err().write("/tmp/mzpzm2p2.xmg");
  // pzmzm2p2.ave_err().write("/tmp/pzmzm2p2.xmg");
  // mxpxm1p1.ave_err().write("/tmp/mxpxm1p1.xmg");
  // pxmxm1p1.ave_err().write("/tmp/pxmxm1p1.xmg");
  
// #define DEFINE_FUNCTOR(NAME,					
// 		       CODE...)					
//   auto NAME=							
//     [&](const std::string& a,					
// 	const std::string& b)					
//     {								
//       return							
// 	[&,							
// 	 d=getCorr(a,b)](const size_t& copy,			
// 			 const size_t& t,			
// 			 const size_t conf)			
// 	{							
// 	  return CODE;
// 	};					
//     }
  
//   DEFINE_FUNCTOR(getter,d[idx({copy,t,conf})]);
  
//   auto copyCombine=
//     [](auto c,
//        const size_t& iConf,
//        const size_t& t,
//        const auto&...obs)
//     {
//       decltype(c(obs(0,0,0)...)) res{};
      
//       for(size_t copy=0;copy<nHits;copy++)
// 	res+=c(obs(copy,t,iConf)...);
      
//       res/=nHits;
      
//       return res;
//     };
  
//   auto jackCombine=
//     [&](auto c,
// 	const auto&...obs)
//     {
//       djvec_t res(T);
      
//       jackknivesFill(nConfs,
// 		     [&](const size_t& iConf,
// 			 const size_t& iClust,
// 			 const double& weight)
// 		     {
// 		       for(size_t t=0;t<T;t++)
// 			 res[t][iClust]+=weight*c(iConf,t,obs...);
// 		     });
		     
//       res.clusterize((double)nConfs/njacks);
      
//       return res;
//     };
  
//   const auto Pi00=
//     getter("ph0_sm0_Sr0_sm0_ph0","ph0_sm0_Sr0_sm0_ph0");

//   const auto PiP1P1R0=getter("phM_smM_Sr0_smP_phP","phP_smP_Sr0_smM_phM");
//   const auto PiP1M1R0=getter("phM_smM_Sr0_smM_phM","phP_smP_Sr0_smP_phP");
//   const auto PiM1M1R0=getter("phP_smP_Sr0_smM_phM","phM_smM_Sr0_smP_phP");
  
//   const auto PiM1M1R1=getter("phP_smP_Sr1_smM_phM","phM_smM_Sr1_smP_phP");
//   const auto PiM1P1R1=getter("phP_smP_Sr1_smP_phP","phM_smM_Sr1_smM_phM");
//   const auto PiP1P1R1=getter("phM_smM_Sr1_smP_phP","phP_smP_Sr1_smM_phM");
  
//   auto real=
//     [&copyCombine](const size_t& iConf,
// 		   const size_t& t,
// 		   const auto& f)
//     {
//       return
// 	copyCombine([](const complex<double>& x)
// 	{
// 	  return x.real();
// 	},iConf,t,f);
//     };
  
//   [[maybe_unused]]
//   auto imag=
//     [&copyCombine](const size_t& iConf,
// 		   const size_t& t,
// 		   const auto& f)
//     {
//       return
// 	copyCombine([](const complex<double>& x)
// 	{
// 	  return x.imag();
// 	},iConf,t,f);
//     };
  
//   auto sub=
//     [&copyCombine](const size_t& iConf,
// 		   const size_t& t,
// 		   const auto& f1,
// 		   const auto& f2)
//     {
//       const complex<double> ab=
// 	copyCombine([](const auto& x,
// 		       const auto& y)
// 	{
// 	  return x*conj(y);
// 	},iConf,t,f1,f2);
      
//       const complex<double> a=
// 	copyCombine([](const auto& x)
// 	{
// 	  return x;
// 	},iConf,t,f1);
      
//       const complex<double> b=
// 	copyCombine([](const auto& x)
// 	{
// 	  return x;
// 	},iConf,t,f2);
      
//       return ((a*conj(b)).real()*nHits-ab.real())/(double)(nHits-1);
//     };
  
//   const djvec_t P00=jackCombine(real,Pi00);
//   const djvec_t P11=jackCombine(real,PiP1P1R0);
//   const djvec_t P11b=jackCombine(real,PiM1P1R1);
//   const djvec_t P11c=jackCombine(real,PiM1M1R0);
//   const djvec_t P11d=jackCombine(real,PiP1P1R1);
//   // const complex<djvec_t> P00ri={jackCombine(real,Pi00),jackCombine(imag,Pi00)};
//   const djvec_t P0000=jackCombine(sub,Pi00,Pi00);
//   const djvec_t PI_UD_DU_mm=jackCombine(sub,PiP1P1R0,PiM1M1R1);
//   const djvec_t PI_UD_DU_pm=jackCombine(sub,PiP1M1R0,PiM1P1R1);
//   const djvec_t PI_UD_DU_pp=jackCombine(sub,PiM1M1R0,PiP1P1R1);
  
//   // const djvec_t P11=
//   //   copyCombine(,Pi11);
  
//   effective_mass(P00.symmetrized()).ave_err().write("plots/Pi00.xmg");
//   effective_mass(P11.symmetrized()).ave_err().write("plots/Pi11.xmg");
//   effective_mass(P11b.symmetrized()).ave_err().write("plots/Pi11b.xmg");
//   effective_mass(P11c.symmetrized()).ave_err().write("plots/Pi11c.xmg");
//   effective_mass(P11d.symmetrized()).ave_err().write("plots/Pi11d.xmg");
//   effective_mass(P0000.symmetrized()).ave_err().write("plots/Pi0000.xmg");
//   effective_mass(PI_UD_DU_mm.symmetrized()).ave_err().write("plots/Pi1111a.xmg");
//   effective_mass(PI_UD_DU_pm.symmetrized()).ave_err().write("plots/Pi1111b.xmg");
//   effective_mass(PI_UD_DU_pp.symmetrized()).ave_err().write("plots/Pi1111c.xmg");
  
//   PI_UD_DU_mm.symmetrized().ave_err().write("plots/cPi1111a.xmg");
//   PI_UD_DU_pm.symmetrized().ave_err().write("plots/cPi1111b.xmg");
//   PI_UD_DU_pp.symmetrized().ave_err().write("plots/cPi1111c.xmg");
  
//   const djvec_t P=effective_mass(P00.symmetrized());
//   const size_t L=64;
//   const djvec_t onePm=sqrt(P*P+sqr(4*M_PI/L));
//   onePm.ave_err().write("plots/1Pm.xmg");
//   const djvec_t twoPm=2*effective_mass(P11.symmetrized());
//   twoPm.ave_err().write("plots/2Pm.xmg");
//   const djvec_t twoPmp=effective_mass(sqr(P11).symmetrized());
//   twoPmp.ave_err().write("plots/2Pmp.xmg");
//   const djvec_t twoP=2*P;
//   twoP.ave_err().write("plots/2P.xmg");
}


void direct()
{
  size_t nHits;
  
  auto getRawDirect=
    [rawData=
     getRaw("rawDirect.dat",
	    "mes_contr_direct_r0_P%s",
	    {"0","1","2","3","4","5"},
	    T,
	    nHits),
     &nHits](const std::string& mso,
	     const std::string& msi)
    {
      const index_t idx({{"t",T},{"conf",nConfs}});
      
      vector<double> res(idx.max());
      vector<double> par(idx.max());
      vector<double> sin(idx.max());
      
      const string what=
	combine("M0_Sr0_P0_0__%s_Sr0_P%s_0,__P5P5",msi.c_str(),mso.c_str());
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	{
	  cout<<"List of corr:"<<endl;
	  for(const auto& [key,val] : rawData)
	    cout<<key<<endl;
	  
	  CRASH("Unable to find %s",what.c_str());
	}
      
      const auto& v=_v->second;
      
      const size_t nHitsToUse=nHits;
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	for(size_t t=0;t<T;t++)
	  {
	    double sumNorm{};
	    double sumProd{};
	    complex<double> sumCompl{};
	    
	    for(size_t iHit=0;iHit<nHitsToUse;iHit++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		sumProd+=(c*c).real();
		sumNorm+=norm(c);
		sumCompl+=c;
	      }
	    
	    res[idx({t,iConf})]=L*L*L*(norm(sumCompl)-sumNorm)/(nHitsToUse*(nHitsToUse-1));
	    par[idx({t,iConf})]=L*L*L*((sumCompl*sumCompl).real()-sumProd)/(nHitsToUse*(nHitsToUse-1));
	    sin[idx({t,iConf})]=sumCompl.real()/nHitsToUse;
	  }
      
      return [data=std::array<std::vector<double>,3>{res,par,sin},
	      idx](const size_t& i,
		   const size_t& t,
		   const size_t& iConf)
      {
	return data[i][idx({t,iConf})];
      };
    };
  
  {
    const size_t iConf=1;
    cout<<confs[iConf]<<endl;
    const InterpDef& bSo = interpDef[0];
    const InterpDef& bSi = interpDef[0];
    
    const std::string& so{bSo.id};
    const std::string& si2{bSi.rap[1]};
    cout<<so<<" "<<si2<<endl;
    
    const auto _A=getRawDirect(so,si2);
    for(size_t t=0;t<T;t++)
      cout<<_A(0,t,iConf)/(L*L*L)<<endl;
  }
  
  auto getDirect=
    [&](const string& mso,
	const string& msi)
    {
      std::array<djvec_t,3> res{djvec_t(T),djvec_t(T),djvec_t(T)};
      
      const auto d=
	getRawDirect(mso,msi);
      
      for(size_t i=0;i<3;i++)
	{
	  for(size_t t=0;t<T;t++)
	    jackknivesFill(nConfs,
			   [&](const size_t& iConf,
			       const size_t& iClust,
			       const double& weight)
			   {
			     res[i][t][iClust]+=weight*d(i,t,iConf);
			   });
	  res[i].clusterize(((double)nConfs/njacks)).symmetrize();
	}
      
      return res;
    };
  
  enum{DIR,PAR,SIN};
  
  for(const InterpDef& bSo : interpDef)
    for(const InterpDef& bSi : interpDef)
      {
	const std::string& repSo{bSo.rep};
	const std::string& repSi{bSi.rep};
	const std::string& so{bSo.id};
	const std::string& si1{bSi.rap[0]};
	const std::string& si2{bSi.rap[1]};
	
	const auto _A=getDirect(so,si1);
	const auto _B=getDirect(so,si2);
	
	const djvec_t A=_A[DIR];
	const djvec_t B=_B[DIR];
	const djvec_t D=_A[PAR];
	const djvec_t E=_B[PAR];
	const djvec_t G=_A[SIN];
	const djvec_t H=_B[SIN];
      // effective_mass(g,T/2).ave_err().write("plots/A_"+repSi+"_"+repSo+".xmg");
      // effective_mass(h,T/2).ave_err().write("plots/B_"+repSi+"_"+repSo+".xmg");
	
	(effective_mass(_B[DIR])-2*effective_mass(_B[SIN])).ave_err().write("plots/dinte_"+repSi+"_"+repSo+".xmg");
	(effective_mass(_B[DIR]-_A[DIR])-2*effective_mass(_B[SIN])).ave_err().write("plots/dinte2_"+repSi+"_"+repSo+".xmg");
	
	effective_mass(A).ave_err().write("plots/dA_"+repSi+"_"+repSo+".xmg");
	effective_mass(B).ave_err().write("plots/dB_"+repSi+"_"+repSo+".xmg");
	effective_mass(D).ave_err().write("plots/dD_"+repSi+"_"+repSo+".xmg");
	effective_mass(E).ave_err().write("plots/dE_"+repSi+"_"+repSo+".xmg");
	effective_mass(G).ave_err().write("plots/dG_"+repSi+"_"+repSo+".xmg");
	effective_mass(H).ave_err().write("plots/dH_"+repSi+"_"+repSo+".xmg");
	effective_mass(A-B).ave_err().write("plots/dC_"+repSi+"_"+repSo+".xmg");
	effective_mass(H-G).ave_err().write("plots/dI_"+repSi+"_"+repSo+".xmg");
    }
}

int main()
{
  setConfs();
  njacks=50;
  
  //box();
  
  direct();
  
  return 0;
}
