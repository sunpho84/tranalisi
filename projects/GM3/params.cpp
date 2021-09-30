#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>

#include <params.hpp>
#include <cachedData.hpp>
#include <rawData.hpp>

void readInput()
{
  raw_file_t input("input.txt","r");
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
  nSourcesMax=input.read<size_t>("NSourcesMax");
  nConfsMax=input.read<size_t>("NConfsMax");
  nLevels=input.read<size_t>("NLevels");
}

void readConfMap()
{
  if(file_exists("map.txt"))
    {
      raw_file_t confMapFile("map.txt","r");
      for(size_t iConf=0;iConf<nConfs;iConf++)
	confMap.push_back(confMapFile.read<int>());
    }
  else
    confMap=vector_up_to<int>(nConfs);
}

void init(int narg,char **arg)
{
  MPI_Init(&narg,&arg);
  
  int temp;
  MPI_Comm_size(MPI_COMM_WORLD,&temp);
  nMPIranks=temp;
  MPI_Comm_rank(MPI_COMM_WORLD,&temp);
  MPIrank=temp;
  
  console.open((MPIrank==0)?"/dev/stdout":"/dev/null");
  
  readInput();
  
  set_njacks(30);
  
  a=make_unique<djack_t>();
  aFromPhysLine=make_unique<djack_t>();
  a->fill_gauss({aAve,aErr,2354235});
  TH=T/2;
  
  THp1=TH+1;
  
  idData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Gamma",16}});
  idOpenData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Id1",4},{"Id2",4},{"Id3",4},{"Id4",4},{"Ri",2}});
}

void close()
{
  if(hasTostoreCachedAveCorr)
    storeCachedAveCorr();
  
  MPI_Finalize();
}
