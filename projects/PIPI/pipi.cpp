#include <tranalisi.hpp>

const string dataPath="data";
vector<string> confs;
size_t nConfs;

size_t T;
size_t nHits=20;
  

inline void setConfs()
{
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
  
  nConfs=confs.size();
}

inline map<string,vector<vector<vector<double>>>> getRawData()
{
  map<string,vector<vector<vector<double>>>> rawData;
  
  for(const filesystem::path conf : confs)
    {
      raw_file_t file(dataPath/conf/"mes_contr_2Pi","r");
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
		  
		  vector<double>& data=confData[tag].emplace_back(2*T);
		  for(size_t t=0;t<T;t++)
		    {
		      if(not readLine())
			CRASH("Unable to read time %zu for contr %s %s",t,a,b);
		      
		      double r,i;
		      if(sscanf(line,"%lg %lg",&r,&i)!=2)
			CRASH("Unable to convert %s to two doubles",line);
		      
		      data[t]=r;
		      data[t+T]=i;
		    }
		}
	    }
	}
      
      for(const auto& [tag,in] : confData)
	{
	  auto& outs=
	    rawData[tag];
	  
	  if(outs.empty() or in.size()==outs.front().size())
	    outs.push_back(in);
	  // else
	  //   cout<<"not considering conf "<<conf<<" "<<in.size()<<"!="<<outs.front().size()<<endl;
	}
    }
  
  if(size_t tmpNConfs=rawData.begin()->second.size();nConfs!=tmpNConfs)
    {
      nConfs=tmpNConfs;
      cout<<"Adjusted nConfs to: "<<nConfs<<endl;
    }

  if(const size_t tmpNHits=rawData.begin()->second.front().size();tmpNHits!=nHits)
    {
      nHits=tmpNHits;
      cout<<"Adjusted nHits to: "<<nHits<<endl;
    }
  
  return rawData;
}

int main()
{
  T=128;
  
  setConfs();
  njacks=50;
  cout<<"NConfs: "<<nConfs<<endl;
  
  const index_t idx({{"hit",nHits},{"T",T},{"conf",nConfs}});
  
  // auto rawData=
  //   getRawData();
  // for(const auto& i : rawData)
  //   cout<<i.first<<endl;
  // return 0;
  
  auto getCorr=
    [rawData=
     getRawData(),
     &idx](const std::string& a,
	   const std::string& b)
    {
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("%s__%s,__P5P5",a.c_str(),b.c_str());
      // cout<<"Searching for "<<what<<endl;
      
      const auto& v=
	rawData.at(what);
      
      for(size_t iHit=0;iHit<nHits;iHit++)
	{
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<T;t++)
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
  
#define DEFINE_FUNCTOR(NAME,					\
		       CODE...)					\
  auto NAME=							\
    [&](const std::string& a,					\
	const std::string& b)					\
    {								\
      return							\
	[&,							\
	 d=getCorr(a,b)](const size_t& copy,			\
			 const size_t& t,			\
			 const size_t conf)			\
	{							\
	  return CODE;						\
	};							\
    }
  
  DEFINE_FUNCTOR(getter,d[idx({copy,t,conf})]);
  
  auto copyCombine=
    [](auto c,
       const size_t& iConf,
       const size_t& t,
       const auto&...obs)
    {
      decltype(c(obs(0,0,0)...)) res{};
      
      for(size_t copy=0;copy<nHits;copy++)
	res+=c(obs(copy,t,iConf)...);
      
      res/=nHits;
      
      return res;
    };
  
  auto jackCombine=
    [&](auto c,
	const auto&...obs)
    {
      djvec_t res(T);
      
      jackknivesFill(nConfs,
		     [&](const size_t& iConf,
			 const size_t& iClust,
			 const double& weight)
		     {
		       for(size_t t=0;t<T;t++)
			 res[t][iClust]+=weight*c(iConf,t,obs...);
		     });
		     
      res.clusterize((double)nConfs/njacks);
      
      return res;
    };
  
  const auto Pi00=
    getter("ph0_sm0_Sr0_sm0_ph0","ph0_sm0_Sr0_sm0_ph0");

  const auto PiP1P1R0=getter("phM_smM_Sr0_smP_phP","phP_smP_Sr0_smM_phM");
  const auto PiP1M1R0=getter("phM_smM_Sr0_smM_phM","phP_smP_Sr0_smP_phP");
  const auto PiM1M1R0=getter("phP_smP_Sr0_smM_phM","phM_smM_Sr0_smP_phP");
  
  const auto PiM1M1R1=getter("phP_smP_Sr1_smM_phM","phM_smM_Sr1_smP_phP");
  const auto PiM1P1R1=getter("phP_smP_Sr1_smP_phP","phM_smM_Sr1_smM_phM");
  const auto PiP1P1R1=getter("phM_smM_Sr1_smP_phP","phP_smP_Sr1_smM_phM");
  
  auto real=
    [&copyCombine](const size_t& iConf,
		   const size_t& t,
		   const auto& f)
    {
      return
	copyCombine([](const complex<double>& x)
	{
	  return x.real();
	},iConf,t,f);
    };
  
  [[maybe_unused]]
  auto imag=
    [&copyCombine](const size_t& iConf,
		   const size_t& t,
		   const auto& f)
    {
      return
	copyCombine([](const complex<double>& x)
	{
	  return x.imag();
	},iConf,t,f);
    };
  
  auto sub=
    [&copyCombine](const size_t& iConf,
		   const size_t& t,
		   const auto& f1,
		   const auto& f2)
    {
      const complex<double> ab=
	copyCombine([](const auto& x,
		       const auto& y)
	{
	  return x*conj(y);
	},iConf,t,f1,f2);
      
      const complex<double> a=
	copyCombine([](const auto& x)
	{
	  return x;
	},iConf,t,f1);
      
      const complex<double> b=
	copyCombine([](const auto& x)
	{
	  return x;
	},iConf,t,f2);
      
      return ((a*conj(b)).real()*nHits-ab.real())/(double)(nHits-1);
    };
  
  const djvec_t P00=jackCombine(real,Pi00);
  const djvec_t P11=jackCombine(real,PiP1P1R0);
  const djvec_t P11b=jackCombine(real,PiM1P1R1);
  const djvec_t P11c=jackCombine(real,PiM1M1R0);
  const djvec_t P11d=jackCombine(real,PiP1P1R1);
  // const complex<djvec_t> P00ri={jackCombine(real,Pi00),jackCombine(imag,Pi00)};
  const djvec_t P0000=jackCombine(sub,Pi00,Pi00);
  const djvec_t PI_UD_DU_mm=jackCombine(sub,PiP1P1R0,PiM1M1R1);
  const djvec_t PI_UD_DU_pm=jackCombine(sub,PiP1M1R0,PiM1P1R1);
  const djvec_t PI_UD_DU_pp=jackCombine(sub,PiM1M1R0,PiP1P1R1);
  
  // const djvec_t P11=
  //   copyCombine(,Pi11);
  
  effective_mass(P00.symmetrized()).ave_err().write("plots/Pi00.xmg");
  effective_mass(P11.symmetrized()).ave_err().write("plots/Pi11.xmg");
  effective_mass(P11b.symmetrized()).ave_err().write("plots/Pi11b.xmg");
  effective_mass(P11c.symmetrized()).ave_err().write("plots/Pi11c.xmg");
  effective_mass(P11d.symmetrized()).ave_err().write("plots/Pi11d.xmg");
  effective_mass(P0000.symmetrized()).ave_err().write("plots/Pi0000.xmg");
  effective_mass(PI_UD_DU_mm.symmetrized()).ave_err().write("plots/Pi1111a.xmg");
  effective_mass(PI_UD_DU_pm.symmetrized()).ave_err().write("plots/Pi1111b.xmg");
  effective_mass(PI_UD_DU_pp.symmetrized()).ave_err().write("plots/Pi1111c.xmg");
  
  PI_UD_DU_mm.symmetrized().ave_err().write("plots/cPi1111a.xmg");
  PI_UD_DU_pm.symmetrized().ave_err().write("plots/cPi1111b.xmg");
  PI_UD_DU_pp.symmetrized().ave_err().write("plots/cPi1111c.xmg");
  
  const djvec_t P=effective_mass(P00.symmetrized());
  const size_t L=64;
  const djvec_t onePm=sqrt(P*P+sqr(4*M_PI/L));
  onePm.ave_err().write("plots/1Pm.xmg");
  const djvec_t twoPm=2*effective_mass(P11.symmetrized());
  twoPm.ave_err().write("plots/2Pm.xmg");
  const djvec_t twoPmp=effective_mass(sqr(P11).symmetrized());
  twoPmp.ave_err().write("plots/2Pmp.xmg");
  const djvec_t twoP=2*P;
  twoP.ave_err().write("plots/2P.xmg");
  
  return 0;
}
