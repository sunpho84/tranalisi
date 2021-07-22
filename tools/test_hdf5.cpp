#include <iostream>

#include <filesystem>
#include <glob.h>
#include <H5Cpp.h>

#include <mpi.h>

#include <vector>

#include <tranalisi.hpp>

using namespace H5;
using namespace std;

int nMPIranks,MPIrank;
size_t T,L,TH;
string confsPattern;
string output;
size_t nConfs,nSources;
constexpr size_t nPV=3;

index_t idData;
index_t idData_loader;
vector<float> rawData;

void setPars()
{
  TH=T/2;
  cout<<"NConfs: "<<nConfs<<endl;
  cout<<"NSources: "<<nSources<<endl;
  idData.set_ranges({{"T",TH+1},{"PV",nPV},{"Confs",nConfs},{"Source",nSources}});
  idData_loader.set_ranges({{"T",T},{"Gamma",16}});
  rawData.resize(idData.max(),0.0);
  cout<<"Data size: "<<rawData.size()<<endl;
}

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

struct DataLoader
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
  
  std::vector<float> dataIn;
  
  DataLoader(const int T)
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
    
    dataIn.resize(T*nGamma);
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
    
    dataset.read(&dataIn[0],PredType::NATIVE_FLOAT,memspace,dataspace);
  }
};

void loadRawData()
{
  const vector<string> confsList=getConfsList(confsPattern);
  nConfs=confsList.size();
  const vector<string> sourcesList=getSourcesList(confsList.front());
  nSources=sourcesList.size();
  
  setPars();
  
  DataLoader loader(T);
  const array<pair<int,int>,7> map{std::pair<int,int>{0,0},{1,6},{1,7},{1,8},{2,10},{2,11},{2,12}};
  
  const size_t confChunk=(nConfs+nMPIranks-1)/nMPIranks;
  const size_t firstConf=std::min((size_t)MPIrank*confChunk,nConfs);
  const size_t lastConf=std::min(firstConf+confChunk,nConfs);
  for(size_t iConf=firstConf;iConf<lastConf;iConf++)
    {
      cout<<MPIrank<<" "<<iConf<<" ["<<firstConf<<":"<<lastConf<<"] "<<confsList[iConf]<<endl;
      
      for(size_t iSource=0;iSource<nSources;iSource++)
	{
	  const string file=confsList[iConf]+"/"+sourcesList[iSource];
	  loader.open(file);
	  loader.load("dd");
	  
	  for(size_t tIn=0;tIn<T;tIn++)
	    for(const auto& m : map)
	      {
		const size_t tOut=(tIn>=TH)?(T-tIn):tIn;
		const size_t& igamma_out=m.first;
		const size_t& igamma_in=m.second;
		const float& in=loader.dataIn[idData_loader({tIn,igamma_in})];
		float& out=rawData[idData({tOut,igamma_out,iConf,iSource})];
		
		out+=in;
	      }
	}
    }
  
  MPI_Allreduce(MPI_IN_PLACE,&rawData[0],idData.max(),MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  
  for(size_t i=0;i<idData.max();i++)
    {
      const std::vector<size_t> coords=idData(i);
      
      const size_t &t=coords[0];
      const size_t &ig=coords[1];
      
      const double norm=((t==0 or t==TH)?1:2)*((ig==0)?1:-3);
      rawData[i]/=norm;
    }
  
  if(MPIrank==0)
    {
      raw_file_t out(output,"w");
      out.bin_write(nConfs);
      out.bin_write(nSources);
      out.bin_write(rawData);
    }
  
  MPI_Finalize();
}

void loadData()
{
  raw_file_t out(output,"r");
  out.bin_read(nConfs);
  out.bin_read(nSources);
  
  setPars();
  
  out.bin_read(rawData);
}

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t icorr)
{
  const size_t clust_size=nConfs/njacks;
  djvec_t ave(TH+1);
  ave=0.0;
  for(size_t iConf=0;iConf<nConfs;iConf++)
    {
      const size_t iClust=iConf/clust_size;
      for(size_t iSource=iSourceMin;iSource<iSourceMax;iSource++)
	for(size_t t=0;t<TH+1;t++)
	  ave[t][iClust]+=rawData[idData({t,icorr,iConf,iSource})];
    }
  
  ave.clusterize(clust_size);
  ave/=(iSourceMax-iSourceMin)*L*L*L;
  
  return ave;
}

int main(int narg,char **arg)
{
  MPI_Init(&narg,&arg);
  MPI_Comm_size(MPI_COMM_WORLD,&nMPIranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);
  
  raw_file_t input("input.txt","r");
  T=input.read<size_t>("T");
  L=input.read<size_t>("L");
  confsPattern=input.read<string>("ConfsPattern");
  output=input.read<string>("Output");
  
  if(not file_exists(output))
    loadRawData();
  else
    loadData();
  
  set_njacks(15);
  
  const djvec_t aveVK=getAve(0,nSources,1);
  const djvec_t aveTK=getAve(0,nSources,2);
  effective_mass(aveVK).ave_err().write("plots/VKVK_new.xmg");
  effective_mass(aveTK).ave_err().write("plots/TKTK_new.xmg");
  
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
