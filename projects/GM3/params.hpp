#ifndef _PARAMS_HPP
#define _PARAMS_HPP

#include <index.hpp>

#include <fstream>
#include <string>

#include <externer.hpp>

using namespace std;

EXTERN size_t nMPIranks,MPIrank INIT_EXTERN_TO(=0);
EXTERN double amq;
EXTERN size_t T,TH,THp1;
EXTERN double L;
EXTERN size_t tMinVKVK,tMaxVKVK;
EXTERN size_t tMinP5P5[2],tMaxP5P5[2];
EXTERN double a,ZaPetros;
EXTERN string confsPattern;
EXTERN string refConfPattern;
EXTERN string rawDataPackedPath;
EXTERN size_t nConfs,nConfsMax,nSources,nSourcesMax;
EXTERN ofstream console;
EXTERN index_t idData_loader;
EXTERN index_t idOpenData_loader;

constexpr size_t nGammaComb=9;
constexpr bool isVK[]=                     {0     ,1     ,1     ,1     ,1     ,0     ,0     ,0     ,0    };
constexpr int parity[]=                    {+1    ,+1    ,+1    ,-1    ,-1    ,-1    ,-1    ,-1    ,-1   };
constexpr char gammaCombTag[nGammaComb][5]={"P5P5","VKVK","TKTK","VKTK","TKVK","A0P5","P5A0","V0P5","P5V0"};
enum CORR_ID                               {idP5P5,idVKVK,idTKTK,idVKTK,idTKVK,idA0P5,idP5A0,idV0P5,idP5V0};

enum RegoType{REGO_TM,REGO_OS};
static const string regoTag[2]={"TM","OS"};

using AveId=array<size_t,4>;

constexpr double corrNorm(const CORR_ID& id)
{
  switch(id)
  {
  case idP5P5:
    return 1.0;
    break;
  case idVKVK:
  case idTKTK:
  case idVKTK:
  case idTKVK:
    return -3.0;
    break;
  case idA0P5:
    return -1.0;
    break;
  case idP5A0:
    return 1.0;
    break;
  case idV0P5:
    return -1.0;
    break;
  case idP5V0:
    return 1.0;
    break;
  }
  
  return 1.0;
}

constexpr size_t nMes=3;
constexpr char mesTag[nMes][3]={"uu","ud","dd"};

void init(int narg,char **arg);
void close();
void readConfMap();

#endif
