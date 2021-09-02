#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <filesystem>

#include <cachedData.hpp>
#include <data.hpp>
#include <params.hpp>
#include "rawData.hpp"

using namespace std;

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes)
{
  djvec_t res;
  
   const AveId id{iSourceMin,iSourceMax,iGammaComb,iMes};
  
   const auto ref=aveCorrCache.find(id);
   if(ref!=aveCorrCache.end())
     res=ref->second;
   else
     {
       res=getRawAve(iSourceMin,iSourceMax,iGammaComb,iMes);
       
       if(iSourceMin==0 and iSourceMax==nSources)
	 {
	   aveCorrCache[id]=res;
	   hasTostoreCachedAveCorr=true;
	 }
     }
  
  return res;
}

djvec_t getTMAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb)
{
  djvec_t ave(THp1);
  for(size_t iMes=0;iMes<nMes;iMes+=2)
    ave+=getAve(iSourceMin,iSourceMax,iGammaComb,iMes);
  
  return ave/2;
}

djvec_t getAveForRego(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const RegoType& rego)
{
  djvec_t res;
  
  switch(rego)
    {
    case REGO_TM:
      res=getTMAve(iSourceMin,iSourceMax,iGammaComb);
      break;
    case REGO_OS:
      res=getAve(iSourceMin,iSourceMax,iGammaComb,1);
    }
  
  return res;
}
