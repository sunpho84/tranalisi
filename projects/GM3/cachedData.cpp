#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <filesystem>

#include <cachedData.hpp>

bool loadCachedAveCorr()
{
  if(not file_exists(cachedAveCorrPath))
    {
      console<<cachedAveCorrPath<<" does not exist, skipping loading cached average correlators"<<endl;
      
      return false;
    }
  
  if(file_exists(rawDataPackedPath) and std::filesystem::last_write_time(cachedAveCorrPath)<std::filesystem::last_write_time(rawDataPackedPath))
    {
      console<<rawDataPackedPath<<" newer than "<<cachedAveCorrPath<<", skipping loading averaged"<<endl;
      
      return false;
    }
  
  raw_file_t file(cachedAveCorrPath,"r");
  console<<"Loading cached average correlators from file "<<cachedAveCorrPath<<endl;
  const size_t cachedNJacks=file.bin_read<size_t>();
  
  //load only if njacks agree
  if(cachedNJacks!=njacks)
    {
      console<<"NJacks in the file is "<<cachedNJacks<<" expecting "<<njacks<<", not loading"<<endl;
      
      return false;
    }
  
  file.bin_read(nConfs);
  file.bin_read(nSources);
  if(nSourcesMax and nSources>nSourcesMax)
    {
      console<<"NSources in the file is "<<nSources<<" is exceeding nSourcesMax, "<<nSourcesMax<<", not loading"<<endl;
      
      return false;
    }
  
  const size_t nCached=file.bin_read<size_t>();
  console<<"Reading  "<<nCached<<" cached average correlators"<<endl;
  
  for(size_t i=0;i<nCached;i++)
    {
      const AveId id=file.bin_read<AveId>();
      
      const size_t n=file.bin_read<size_t>();
      djvec_t data(n);
      file.bin_read(data);
      aveCorrCache[id]=data;
    }
  
  return true;
}

void storeCachedAveCorr()
{
  const size_t& nCached=aveCorrCache.size();
  
  console<<"Storing "<<nCached<<" cached average correlators to file "<<cachedAveCorrPath<<endl;
  
  raw_file_t file(cachedAveCorrPath,"w");
  file.bin_write(njacks);
  file.bin_write(nConfs);
  file.bin_write(nSources);
  file.bin_write(nCached);
  
  for(const auto& c : aveCorrCache)
    {
      file.bin_write(c.first);
      file.bin_write(c.second.size());
      file.bin_write(c.second);
    }
}

