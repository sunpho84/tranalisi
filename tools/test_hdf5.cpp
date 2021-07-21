#include <iostream>

#include <filesystem>
#include <glob.h>
#include <H5Cpp.h>
#include <vector>

#include <tranalisi.hpp>

using namespace H5;
using namespace std;

vector<string> getConfsList(const string& confsPattern)
{
  vector<string> confsList;
  glob_t globbuf;
  
  if(glob(confsPattern.c_str(),0,nullptr,&globbuf))
    CRASH("Unable to find pattern %s for conf",confsPattern.c_str());
  else
    for(int j=0;j<(int)globbuf.gl_pathc;j++)
      confsList.push_back(globbuf.gl_pathv[j]);
  globfree(&globbuf);
  
  return confsList;
}

vector<string> getSourcesList(const string& firstConf)
{
  vector<string> sourcesList;
  glob_t globbuf;
  
  const string sourcesPattern=(firstConf+"/twop_id*_st*.h5");
  if(glob(sourcesPattern.c_str(),0,nullptr,&globbuf))
    CRASH("Unable to find pattern %s for source",sourcesPattern.c_str());
  else
    for(int j=0;j<(int)globbuf.gl_pathc;j++)
      {
	const std::filesystem::path confFull=globbuf.gl_pathv[j];
	sourcesList.push_back(confFull.filename());
      }
  globfree(&globbuf);
  
  return sourcesList;
}

struct dataLoader
{
  //const int rank=dataspace.getSimpleExtentNdims();
  static constexpr int rank=4;
  static constexpr int nGamma=16;
  
  hsize_t dimsm[2];
  DataSpace memspace;
  
  hsize_t      offset[rank];   // hyperslab offset in the file
  hsize_t      count[rank];    // size of the hyperslab in the file
  
  H5File file;
  string groupName;
  
  std::vector<float> data_out;
  
  dataLoader(const int T)
  {
    dimsm[0]=T;
    dimsm[1]=nGamma;
    memspace.setExtentSimple(2,dimsm);
    
    for(int i=0;i<rank;i++)
      offset[i]=0;
    
    count[0]=T;
    count[1]=1;
    count[2]=nGamma;
    count[3]=1;
    
    data_out.resize(T*nGamma);
  }
  
  void open(const string& path)
  {
    file.openFile(path,H5F_ACC_RDONLY);
    groupName=file.getObjnameByIdx(0);
  }
  
  void close()
  {
    file.close();
  }
  
  void load(const string& tag)
  {
    const DataSet dataset=file.openDataSet(groupName+"/mesons/"+tag);
    const DataSpace dataspace=dataset.getSpace();
    
    dataspace.selectHyperslab(H5S_SELECT_SET,count,offset);
    
    dataset.read(&data_out[0],PredType::NATIVE_FLOAT,memspace,dataspace);
  }
};

int main()
{
  raw_file_t input("input.txt","r");
  const int T=input.read<size_t>("T");
  const int TH=T/2;
  const string confsPattern=input.read<string>("ConfsPattern");
  const vector<string> confsList=getConfsList(confsPattern);
  cout<<"NConfs: "<<confsList.size()<<endl;
  const vector<string> sourcesList=getSourcesList(confsList.front());
  cout<<"NSources: "<<sourcesList.size()<<endl;
  
  constexpr int nPV=2;
  index_t data({{"T",TH+1},{"PV",nPV}});
  
  dataLoader loader(T);
  
  for(size_t iConf=0;iConf<confsList.size();iConf++)
    for(size_t iSource=0;iSource<sourcesList.size();iSource++)
      {
	const string file=confsList[iConf]+"/"+sourcesList[iSource];
	cout<<file<<endl;
	loader.open(file);
	loader.load("uu");
    }
  
  // cout<<"rank: "<<rank<<endl;
  // hsize_t dims_out[rank];
  // dataspace.getSimpleExtentDims(dims_out,nullptr);
  // for(int i=0;i<rank;i++)
  //   cout<<dims_out[i]<<endl;
  
  // for(int t=0;t<T;t++)
  //   cout<<loader.data_out[0+nGamma*t]<<endl;
  
  // const DataType dataType=dataset.getDataType();
  
  // dataset.q
  // cout<<"NSubobj:"<<dataset.getNumObjs()<<endl;
  // for(ssize_t i=0;i<dataset.getNumObjs();i++)
  //   {
  //     const string subGroupName=dataset.getObjnameByIdx(i);
  //     cout<<"Subobj:"<<subGroupName<<endl;
  //   }
  
  return 0;
}
