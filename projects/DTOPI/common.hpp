#include "effective.hpp"
#include "fit.hpp"
#include "functions.hpp"
#include "grace.hpp"
#include "index.hpp"
#include "jack.hpp"
#include "math.hpp"
#include "meas_vec.hpp"
#include "oper.hpp"
#include "random.hpp"
#include "raw_file.hpp"
#include <filesystem>
#include <string>
#include <vector>

const string dataPath="data";
inline vector<string> confs;
inline size_t nConfs;

inline size_t T;

inline void setConfs()
{
  map<pair<int,int>,string> tmpConfs;
  for(const auto& entry : filesystem::directory_iterator(dataPath))
    if(entry.is_directory())
      //if(filesystem::exists(entry.path()/"mes_contr_3pUnsme"))
      {
	const string& conf=filesystem::path(entry.path()).filename();
	const size_t iStream=atoi(conf.substr(6,1).c_str());
	const size_t iConf=atoi(conf.substr(0,4).c_str());
	const int par=(iStream%2)*2-1;
	const int offs=iStream/2*2;
	tmpConfs[{offs,iConf*par}]=conf;
      }
  
  for(const auto& [dum,conf] : tmpConfs)
    // {
    confs.emplace_back(conf);
    //   cout<<conf<<endl;
    // }
  
  nConfs=confs.size();
}

inline map<string,vector<vector<double>>> getRawData(const std::string& entry)
{
  map<string,vector<vector<double>>> rawData;
  
  for(const filesystem::path conf : confs)
    {
      raw_file_t file(dataPath/conf/("mes_contr_"+(string)entry),"r");
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
		CRASH("Unable to get the corr name");
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
