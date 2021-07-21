#include <iostream>

#include <glob.h>
#include <H5Cpp.h>
#include <vector>

#include <tranalisi.hpp>

using namespace H5;
using namespace std;

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
  const string confPattern=input.read<string>("ConfsPattern");
  
  dataLoader loader(T);

  vector<string> conf_list;
  glob_t globbuf;
  if(glob(confPattern.c_str(),0,nullptr,&globbuf))
    CRASH("Unable to find pattern %s for conf",confPattern.c_str());
  else
    for(int j=0;j<(int)globbuf.gl_pathc;j++)
      conf_list.push_back(globbuf.gl_pathv[j]);
  globfree(&globbuf);
  
  constexpr int nPV=2;
  index_t data({{"T",TH+1},{"PV",nPV}});
  
  for(size_t iConf=0;iConf<conf_list.size();iConf++)
    {
      loader.open(conf_list[iConf]+"twop_id178_st050.h5");
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
