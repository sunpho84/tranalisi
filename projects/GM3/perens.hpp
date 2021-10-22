#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <memory>

#include <index.hpp>
#include <jack.hpp>
#include <meas_vec.hpp>

#include "params.hpp"
#include "VKVKRepresentation.hpp"

struct perens_t
{
  double amq;
  size_t T,TH,THp1;
  double L;
  size_t tMinVKVK,tMaxVKVK;
  size_t tMinP5P5[2],tMaxP5P5[2];
  size_t tMinLusch[2];
  double aAve,aErr,ZaPetros;
  djack_t a;
  djack_t aFromPhysLine;
  string confsPattern;
  string refConfPattern;
  string rawDataPackedPath;
  size_t nConfs,nConfsMax,nSources,nSourcesMax;
  
  index_t idData_loader;
  index_t idOpenData_loader;
  int nLevels;
  
  index_t idData;
  vector<double> _rawData;
  vector<int> confMap;
  vector<string> possibleConfsList;
  vector<string> sourcesList;
  vector<string> confsList;
  
  djvec_t Z;
  
  bool canUseRawData{false};
  
  mutable map<AveId,djvec_t> aveCorrCache;
  
  mutable bool hasTostoreCachedAveCorr{false};
  
  const string basePath;
  
  void readInput()
  {
    raw_file_t input(basePath+"/input.txt","r");
    
    T=input.read<size_t>("T");
    L=input.read<double>("L");
    amq=input.read<double>("amq");
    aAve=input.read<double>("a");
    aErr=input.read<double>();
    ZaPetros=input.read<double>("Za");
    confsPattern=input.read<string>("ConfsPattern");
    refConfPattern=input.read<string>("RefConfPattern");
    rawDataPackedPath=input.read<string>("Output");
    input.expect("TFitP5P5");
    for(size_t iMes=0;iMes<2;iMes++)
      {
	tMinP5P5[iMes]=input.read<size_t>();
	tMaxP5P5[iMes]=input.read<size_t>();
      }
    tMinVKVK=input.read<size_t>("TFitVKVK");
    tMaxVKVK=input.read<size_t>();
    tMinLusch[0]=input.read<size_t>("TMinLusch");
    tMinLusch[1]=input.read<size_t>();
    nSourcesMax=input.read<size_t>("NSourcesMax");
    nConfsMax=input.read<size_t>("NConfsMax");
    nLevels=input.read<size_t>("NLevels");
    
    a.fill_gauss({aAve,aErr,2354235});
    TH=T/2;
    
    THp1=TH+1;
    
    idData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Gamma",16}});
    idOpenData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Id1",4},{"Id2",4},{"Id3",4},{"Id4",4},{"Ri",2}});
  }
  
  void readConfMap()
  {
    if(file_exists(basePath+"/map.txt"))
      {
	raw_file_t confMapFile("map.txt","r");
	for(size_t iConf=0;iConf<nConfs;iConf++)
	  confMap.push_back(confMapFile.read<int>());
      }
    else
      confMap=vector_up_to<int>(nConfs);
  }
  
  perens_t(const string& basePath) :
    basePath(basePath)
  {
    readInput();
  }
  
  ~perens_t()
  {
    if(hasTostoreCachedAveCorr)
    storeCachedAveCorr();
  }
  
  const double& rawData(const size_t& _iConf,const size_t& iSource,const size_t& igamma_out,const size_t& iMes,const size_t& tOut) const;
  
  double& rawData(const size_t& _iConf,const size_t& iSource,const size_t& igamma_out,const size_t& iMes,const size_t& tOut);
  
  void setRawData(const size_t& nConfsToRes);
  
  string sourceName(const size_t& iConf,const size_t& iSource);
  
  void loadAndPackRawData(int narg,char** arg);
  
  void rawDataAn(const size_t& iGammaComb);
  
  void analyzeRawData();
  
  void loadRawData(int narg,char **arg);
  
  djvec_t getRawAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes) const;
  
  void convertForSilvano() const;
  
  void computeAmu(const RegoType& rego,
		  const jack_t<VKVKRepFiniteVol>& rep,
		  const jack_t<VKVKRepInfiniteVol>& infVolRep);
  
  void Gilberto() const;
  
  bool loadCachedAveCorr();
  
  void storeCachedAveCorr() const;
  
  djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes) const;
  
  djvec_t getTMAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb) const;
  
  djvec_t getAveForRego(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const RegoType& rego) const;
  
  void loadData(int narg,char **arg);
  
  void loadDisco();
  
  void computeMPCAC();
  
  djvec_t determineRenoConst() const;
  
  double VKVKInfVol(const VKVKRepInfiniteVol& rep,
		    const ALaLuscherRepresentationInfVol& aLaLusch,
		    const double& t) const;
  
  pair<jack_t<VKVKRepFiniteVol>,
     jack_t<VKVKRepInfiniteVol>> fitVKVK(const RegoType& rego,
					 const int& nLevels) const;
};

#endif
