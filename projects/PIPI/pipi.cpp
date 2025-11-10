#include <tranalisi.hpp>

const string dataPath="data";
vector<string> confs;
size_t nConfs;

size_t T;

inline void setConfs()
{
  map<pair<int,int>,string> tmpConfs;
  for(const auto& entry : filesystem::directory_iterator(dataPath))
    if(not entry.is_directory())
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

inline map<string,vector<vector<double>>> getRawData()
{
  map<string,vector<vector<double>>> rawData;
  
  for(const filesystem::path conf : confs)
    {
      raw_file_t file(dataPath/conf,"r");
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
      
      char a[100]{},b[100]{},c[100];
      string baseTag;
      string tag;
      while(readLine())
	{
	  if(sscanf(line," # Contraction of %s ^ \\dag and %s",a,b)==2)
	    baseTag=tag=(string)a+"__"+b;
	  else
	    {
	      if(sscanf(line," # %s",c)!=1)
		CRASH("Unable to get the corr name on line: %s",line);
	      else
		{
		  tag=baseTag+"__"+c;
		  
		  vector<double>& data=rawData[tag].emplace_back(2*T);
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
    }
  
  return rawData;
}

int main()
{
  T=128;
  
  setConfs();
  njacks=70;
  cout<<"NConfs: "<<nConfs<<endl;
  
  const size_t nCopies=20;
  
  const index_t idx({{"copy",nCopies},{"T",T},{"conf",nConfs}});
  
  auto getCorr=
    [rawData=
     getRawData(),
     &idx](const std::string& a,
	   const std::string& b)
    {
      vector<complex<double>> res(idx.max());
      
      for(size_t copy=0;copy<nCopies;copy++)
	{
	  const auto& v=
	    rawData.at(combine("%s_copy%zu__%s_copy%zu,__P5P5",a.c_str(),copy,b.c_str(),copy));
	  
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<T;t++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][t+T*ri];
		
		res[idx({copy,t,iConf})]+=c;
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
      
      for(size_t copy=0;copy<nCopies;copy++)
	res+=c(obs(copy,t,iConf)...);
      
      res/=nCopies;
      
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
    getter("lr0_0_sm","SM_lr0_0_sm");
  const auto Pi11=
    getter("ur0_sm","SM_lr0_sm");
  
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
      
      return ((a*conj(b)).real()*nCopies-ab.real())/(double)(nCopies-1);
    };
  
  const djvec_t P00=jackCombine(real,Pi00);
  const complex<djvec_t> P00ri={jackCombine(real,Pi00),jackCombine(imag,Pi00)};
  const complex<djvec_t> P11ri={jackCombine(real,Pi11),jackCombine(imag,Pi11)};
  const djvec_t P1111=jackCombine(sub,Pi11,Pi11)-(P11ri*conj(P11ri)).real();
  
  const auto getSimple=
    [&getCorr,
     &idx](const std::string& a,
	   const std::string& b)
    {
      const vector<complex<double>> g=
	getCorr(a,b);
      
      djvec_t jpi1111(T);
      
      complex<djvec_t> jpi1;
      djvec_t& jpi1r=*(djvec_t*)&jpi1;
      djvec_t& jpi1i=*((djvec_t*)&jpi1+1);
      jpi1r.resize(T);
      jpi1i.resize(T);
      for(size_t t=0;t<T;t++)
	{
	  std::vector<double> pi11(nConfs,0.0);
	  std::vector<double> pi1r(nConfs,0.0);
	  std::vector<double> pi1i(nConfs,0.0);
	  
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    {
	      double s11=0;
	      complex<double> s1=0;
	      
	      for(size_t iCopy=0;iCopy<nCopies;iCopy++)
		{
		  const complex<double>& x=
		    g[idx({iCopy,t,iConf})];
		  
		  s11+=norm(x);
		  s1+=x;
		}
	      
	      s11-=norm(s1);
	      
	      s11/=nCopies*(nCopies-1);
	      s1/=nCopies;
	      
	      pi11[iConf]=s11;
	      pi1r[iConf]=s1.real();
	      pi1i[iConf]=s1.imag();
	    }
	  
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
			 {
			   jpi1111[t][iClust]+=weight*pi11[iConf];
			   jpi1.real()[t][iClust]+=weight*pi1r[iConf];
			   jpi1.imag()[t][iClust]+=weight*pi1i[iConf];
			 });
	}
      
      jpi1111.clusterize((double)nConfs/njacks);
      jpi1r.clusterize((double)nConfs/njacks);
      jpi1i.clusterize((double)nConfs/njacks);
      
      jpi1111-=norm(jpi1);
      
      return jpi1111;
    };
  
  const djvec_t jpi0000=getSimple("lr0_0_sm","SM_lr0_0_sm");
  const djvec_t jpi1111=getSimple("ur0_sm","SM_lr0_sm");
  
  // const djvec_t P11=
  //   copyCombine(,Pi11);
  
  effective_mass(P00.symmetrized()).ave_err().write("plots/Pi00.xmg");
  effective_mass(P1111.symmetrized()).ave_err().write("plots/Pi1111.xmg");
  //effective_mass(P11.symmetrized()).ave_err().write("plots/Pi11.xmg");
  
  effective_mass(jpi1111.symmetrized()).ave_err().write("plots/Pi1111check.xmg");
  effective_mass(jpi0000.symmetrized()).ave_err().write("plots/Pi0000check.xmg");
  const djvec_t P=effective_mass(P00.symmetrized());
  const size_t L=64;
  const djvec_t twoPm=2*sqrt(P*P+sqr(2*M_PI/L));
  twoPm.ave_err().write("plots/2Pm.xmg");
  
  return 0;
}
