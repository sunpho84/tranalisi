#ifndef _PARAMS_HPP
#define _PARAMS_HPP

#include <index.hpp>

#include <fstream>
#include <memory>
#include <string>

#include "externer.hpp"

using namespace std;

EXTERN size_t nMPIranks,MPIrank INIT_EXTERN_TO(=0);
EXTERN ofstream console;

constexpr size_t nGammaComb=9;
constexpr bool isVK[]=                     {0     ,1     ,1     ,1     ,1     ,0     ,0     ,0     ,0    };
constexpr int parity[]=                    {+1    ,+1    ,+1    ,-1    ,-1    ,-1    ,-1    ,-1    ,-1   };
constexpr char gammaCombTag[nGammaComb][5]={"P5P5","VKVK","TKTK","VKTK","TKVK","A0P5","P5A0","V0P5","P5V0"};
enum CORR_ID                               {idP5P5,idVKVK,idTKTK,idVKTK,idTKVK,idA0P5,idP5A0,idV0P5,idP5V0};

enum RegoType{REGO_TM,REGO_OS};
static const string regoTag[2]={"TM","OS"};

constexpr size_t regoZId[2]={1,0};

enum class LoadMethod{Compact,Extended};
EXTERN LoadMethod loadMethod;

inline LoadMethod convertTagToLoadMethod(const string method)
{
  if(method=="compact")
    return LoadMethod::Compact;
  
  if(method=="extended")
    return LoadMethod::Extended;
  
  CRASH("Unknown method %s, use compact or extended",method.c_str());
  
  return {};
}

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

#endif
