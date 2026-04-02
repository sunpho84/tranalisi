#include <tranalisi.hpp>

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

inline std::vector<std::string> getConfs(const std::string& confsListPath,
					 const std::string& rawDataPath)
{
  std::vector<std::string> confs;
  
  if(not file_exists(confsListPath))
    {
      cout<<"Searching for confs in the \""<<rawDataPath<<"\" directory"<<endl;
      
      std::map<std::pair<int,int>,std::string> tmpConfs;
      for(const auto& entry : filesystem::directory_iterator(rawDataPath))
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
      cout<<"Reading conf list from "<<confsListPath<<endl;
      raw_file_t(confsListPath,"r").bin_read(confs);
    }
  
  return confs;
}

inline map<string,vector<vector<vector<double>>>> getRaw(const char* cachedFilePath,
							 const char* rawFileNameTemplate,
							 const std::vector<const char*>& suffixList,
							 const size_t& tMax,
							 const std::string& rawDataDir,
							 const std::vector<std::string> confs)
{
  const size_t nConfs=confs.size();
  
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
	    raw_file_t file(rawDataDir+"/"+conf.string()+"/"+combine(rawFileNameTemplate,suffix),"r");
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
		  {
		    cout<<"Problem in conf "<<confs[iC]<<endl;
		    cout<<"First conf: "<<outs.empty()<<endl;
		    cout<<"Size: "<<in.size()<<endl;
		    cout<<"Expected size: "<<outs.front().size()<<endl;
		    cout<<"Tag: "<<tag<<endl;
		    CRASH("Please remove conf %s",confs[iC].c_str());
		  }
		
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
  
  //cout<<"Setting nHits to: "<<nHits<<endl;
  
  return rawData;
}

std::vector<djvec_t> computeOrLoad(const index_t& idOut,
				   const std::string& path,
				   std::vector<djvec_t> compute(const index_t& idOut))
{
  std::vector<djvec_t> res;
  
  if(file_exists(path))
    {
      cout<<"Loading all data from "<<path<<endl;
      res=raw_file_t(path,"r").bin_read<std::vector<djvec_t>>();
    }
  else
    {
      res=compute(idOut);
      cout<<"Writing all data to "<<path<<endl;
      raw_file_t(path,"w").bin_write(res);
    }
  
  return res;
}

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeBox(const index_t& idOut)
{
  const std::string boxCorrPath="boxCorr";
  const std::vector<std::string> confs=getConfs("confsBoxList.dat",boxCorrPath);
  
  const size_t nConfs=confs.size();
  
  const auto rawData=
    getRaw("rawBox.dat",
	   "mes_contr_box_src_snk%s_T25",
	   {"1","2","3","4","5"},
	   tMaxBox,
	   boxCorrPath,
	   confs);
  const size_t nHits=rawData.begin()->second.front().size();
  
  const index_t idx({{"hit",nHits},{"tMax",tMaxBox},{"conf",nConfs}});
  
  auto getRawBox=
    [&rawData,
     &nHits,
     &idx,
     &nConfs](const std::string& a,
	      const std::string& b)
    {
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
      
      return res;
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
			   o+=d[idx({iHit,t,iConf})].real();
			 o/=nHits;
			 res[t][iClust]+=weight*o;
		       });
	}
      res.clusterize(((double)nConfs/njacks));
      
      return res;
    };
  
  std::vector<djvec_t> res(idOut.max());
  
  for(size_t iBSo=0;iBSo<interpDef.size();iBSo++)
    for(size_t iBSi=0;iBSi<interpDef.size();iBSi++)
      {
	const InterpDef& bSo{interpDef[iBSo]};
	const InterpDef& bSi{interpDef[iBSi]};
	
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
	const djvec_t C=A-B;
	C.ave_err().write("plots/C_"+repSi+"_"+repSo+".xmg");
	
	res[idOut({iBSo,iBSi})]=C;
      }
  
  return res;
}

auto box()
{
  index_t idOut({{"bSo",interpDef.size()},{"bSi",interpDef.size()}});
  
  return [data=computeOrLoad(idOut,"box.dat",computeBox),
	  idOut](const size_t& iBso,
		 const size_t& iBsi) -> const djvec_t&
  {
    return data[idOut({iBso,iBsi})];
  };
}

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeDirect(const index_t& idOut)
{
  const std::string directDataPath="directCorr";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",directDataPath);
  const size_t nConfs=confs.size();
  
  const auto rawData=
     getRaw("rawDirect.dat",
	    "mes_contr_direct_r0_P%s",
	    {"0","1","2","3","4","5"},
	    T,
	    directDataPath,
	    confs);
  
  const size_t nHits=rawData.begin()->second.front().size();
  cout<<"NHits: "<<nHits<<endl;
  
  auto getRawDirect=
    [&rawData,
     &nHits,
     &nConfs](const std::string& mso,
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
  
  // {
  //   const size_t iConf=1;
  //   cout<<confs[iConf]<<endl;
  //   const InterpDef& bSo = interpDef[0];
  //   const InterpDef& bSi = interpDef[0];
    
  //   const std::string& so{bSo.id};
  //   const std::string& si2{bSi.rap[1]};
  //   cout<<so<<" "<<si2<<endl;
    
  //   const auto _A=getRawDirect(so,si2);
  //   for(size_t t=0;t<T;t++)
  //     cout<<_A(0,t,iConf)/(L*L*L)<<endl;
  // }
  
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
  
  std::vector<djvec_t> res(idOut.max());
  
  for(size_t iBSo=0;iBSo<interpDef.size();iBSo++)
    for(size_t iBSi=0;iBSi<interpDef.size();iBSi++)
      {
	const InterpDef& bSo{interpDef[iBSo]};
	const InterpDef& bSi{interpDef[iBSi]};
	
	const std::string& so{bSo.id};
	const std::string& si1{bSi.rap[0]};
	const std::string& si2{bSi.rap[1]};
	
	const auto _A=getDirect(so,si1);
	const auto _B=getDirect(so,si2);
	
	for(size_t iCombo=0;iCombo<3;iCombo++)
	  res[idOut({iBSo,iBSi,iCombo})]=_A[iCombo]-_B[iCombo];
	
	const djvec_t A=_A[DIR];
	const djvec_t B=_B[DIR];
	const djvec_t D=_A[PAR];
	const djvec_t E=_B[PAR];
	const djvec_t G=_A[SIN];
	const djvec_t H=_B[SIN];
	
	const std::string& repSo{bSo.rep};
	const std::string& repSi{bSi.rep};
	(A-B).ave_err().write("plots/dC_"+repSi+"_"+repSo+".xmg");
	
	(effective_mass(_B[DIR])-2*effective_mass(_B[SIN])).ave_err().write("plots/effDinte_"+repSi+"_"+repSo+".xmg");
	(effective_mass(_B[DIR]-_A[DIR])-2*effective_mass(_B[SIN])).ave_err().write("plots/effDinte2_"+repSi+"_"+repSo+".xmg");
	
	effective_mass(A).ave_err().write("plots/effDA_"+repSi+"_"+repSo+".xmg");
	effective_mass(B).ave_err().write("plots/effDB_"+repSi+"_"+repSo+".xmg");
	effective_mass(D).ave_err().write("plots/effDD_"+repSi+"_"+repSo+".xmg");
	effective_mass(E).ave_err().write("plots/effDE_"+repSi+"_"+repSo+".xmg");
	(2*effective_mass(G)).ave_err().write("plots/effDG_"+repSi+"_"+repSo+".xmg");
	(2*effective_mass(H)).ave_err().write("plots/effDH_"+repSi+"_"+repSo+".xmg");
	effective_mass(A-B).ave_err().write("plots/effDC_"+repSi+"_"+repSo+".xmg");
	effective_mass(H-G).ave_err().write("plots/effDI_"+repSi+"_"+repSo+".xmg");
    }
  
  return res;
}

auto direct()
{
  index_t idOut{{{"bSo",interpDef.size()},{"bSi",interpDef.size()},{"combo",3}}};
  
  return [data=computeOrLoad(idOut,"direct.dat",computeDirect),
	  idOut](const size_t& iBso,
		 const size_t& iBsi,
		 const size_t& iCombo) -> const djvec_t&
  {
    return data[idOut({iBso,iBsi,iCombo})];
  };
}

int main()
{
  njacks=50;
  
  const auto b=box();
  
  const auto d=direct();
  
  for(size_t ibSo=0;ibSo<interpDef.size();ibSo++)
    for(size_t ibSi=0;ibSi<interpDef.size();ibSi++)
      {
	b(ibSo,ibSi).ave_err().write(combine("plots/box_%zu_%zu.xmg",ibSo,ibSi));
	d(ibSo,ibSi,0).ave_err().write(combine("plots/dir_%zu_%zu.xmg",ibSo,ibSi));
      }
  
  return 0;
}
