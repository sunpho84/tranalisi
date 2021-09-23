#define EXTERN
#define INIT_EXTERN_TO(ARGS...) ARGS

#include <aLaLuscherRepresentation.hpp>
#include <amu.hpp>
#include <cachedData.hpp>
#include <data.hpp>
#include <fit.hpp>
#include <kernels.hpp>
#include <mPCAC.hpp>
#include <params.hpp>
#include <rawData.hpp>
#include <renoConstants.hpp>
#include <VKVKRepresentation.hpp>

using namespace std;

int main(int narg,char **arg)
{
  init(narg,arg);
  
  loadData(narg,arg);
  
  if(MPIrank==0)
    {
      Z=determineRenoConst();
      
      if(canUseRawData)
	{
	  //analyzeRawData();
	  convertForSilvano();
	}
      
      const int nLevels=6;
      const int tMin=5;
      const jack_t<VKVKRep> rep=
	fitVKVK(nLevels,tMin);
      
      computeAmu(REGO_TM,rep);
      
      computeMPCAC();
      
    }
  
  close();
  
  return 0;
}
