#define EXTERN
#define INIT_EXTERN_TO(ARGS...) ARGS

#include <amu.hpp>
#include <cachedData.hpp>
#include <data.hpp>
#include <kernels.hpp>
#include <params.hpp>
#include <rawData.hpp>
#include <renoConstants.hpp>

using namespace std;

int main(int narg,char **arg)
{
  init(narg,arg);
  
  readInput();
  
  set_njacks(30);
  
  setPars();
  
  if(not loadCachedAveCorr())
    loadRawData(narg,arg);
  
  readConfMap();
  
  clustSize=(double)nConfs/njacks;
  
  const djvec_t Z=determineRenoConst();
  
  if(canUseRawData)
    {
      analyzeRawData();
      convertForSilvano();
    }
  
  /////////////////////////////////////////////////////////////////
  
  computeAmu(Z);
  
  close();
  
  return 0;
}
