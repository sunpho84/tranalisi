#include <tranalisi.hpp>

const size_t T=128;
const size_t L=64;

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
							 const std::vector<std::string> confs,const std::vector<std::string>& corrsName={"P5P5"})
{
  const size_t nConfs=confs.size();
  cout<<"nConfs: "<<confs.size()<<endl;
  
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
			
			if(std::find(corrsName.begin(),corrsName.end(),c)!=corrsName.end())
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

inline std::vector<djvec_t> computeOrLoad(const index_t& idOut,
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

