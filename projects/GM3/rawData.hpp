#ifndef _RAWDATA_HPP
#define _RAWDATA_HPP

#include <externer.hpp>

#include "index.hpp"
#include "meas_vec.hpp"

#include <params.hpp>

EXTERN bool canUseRawData INIT_EXTERN_TO(=false);

EXTERN index_t idData;
EXTERN vector<double> _rawData;
EXTERN vector<int> confMap;
EXTERN vector<string> possibleConfsList;
EXTERN vector<string> sourcesList;

inline double& rawData(const size_t& _iConf,const size_t& iSource,const size_t& igamma_out,const size_t& iMes,const size_t& tOut)
{
  if(not canUseRawData)
    CRASH("Cannot use raw data");
  
  return _rawData[idData({_iConf,iSource,igamma_out,iMes,tOut})];
}

inline void setRawData(const size_t& nConfsToRes)
{
  idData.set_ranges({{"Confs",nConfsToRes},{"Source",nSources},{"GammComb",nGammaComb},{"Mes",nMes},{"T",THp1}});
  _rawData.resize(idData.max(),0.0);
}

djvec_t getRawAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes);
void loadAndPackRawData(int narg,char** arg);
void analyzeRawData();
void loadRawData(int narg,char **arg);
void convertForSilvano();

#endif
