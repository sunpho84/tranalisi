#include <iostream>

#include <filesystem>
#include <glob.h>
#include <H5Cpp.h>

#include <mpi.h>

#include <vector>

#include <tranalisi.hpp>

using namespace H5;
using namespace std;

int nMPIranks,MPIrank=0;
size_t T,L,TH,THp1;
string confsPattern;
string output;
size_t nConfs,nSources;
constexpr size_t nPV=3;

bool tAve;
index_t idData;
index_t idData_loader;
vector<float> rawData;
vector<int> confMap;

void setPars()
{
  TH=T/2;
  
  if(tAve)
    THp1=TH+1;
  else
    THp1=T;
  
  idData.set_ranges({{"T",THp1},{"PV",nPV},{"Confs",nConfs},{"Source",nSources}});
  idData_loader.set_ranges({{"T",T},{"Gamma",16}});
  rawData.resize(idData.max(),0.0);
  
  if(MPIrank==0)
    {
      cout<<"NConfs: "<<nConfs<<endl;
      cout<<"NSources: "<<nSources<<endl;
      cout<<"Data size: "<<rawData.size()<<endl;
    }
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
  //Here for future memory
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

void loadRawData(int narg,char** arg)
{
  MPI_Init(&narg,&arg);
  MPI_Comm_size(MPI_COMM_WORLD,&nMPIranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);
  
  const vector<string> possibleConfsList=getConfsList(confsPattern);
  const vector<string> sourcesList=getSourcesList(possibleConfsList.front());
  vector<string> confsList;
  for(const string& conf : possibleConfsList)
    {
      bool exists=true;
      string confPath;
      size_t iSource=0;
      while(iSource<nSources and exists)
	{
	  confPath=conf+"/"+sourcesList[iSource];
	  if(not file_exists(confPath))
	    {
	      exists=false;
	      if(MPIrank==0)
		cout<<" "<<confPath<<" not found"<<endl;
	    }
	  iSource++;
	}
      
      if(exists)
	  confsList.push_back(conf);
      if(MPIrank==0)
	cout<<"Conf "<<conf<<(exists?" accepted":(" discarded, "+confPath+" not found"))<<endl;
    }
  
  if(MPIrank==0 and confsList.size()!=possibleConfsList.size())
    cout<<"Confs list resized from "<<possibleConfsList.size()<<"to: "<<confsList.size();
  
  nConfs=confsList.size();
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
		const size_t tOut=
		  tAve?
		  ((tIn>=TH)?(T-tIn):tIn):
		  tIn;
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
      
      const double norm=((t==0 or t==TH or tAve==false)?1:2)*
	((ig==0)?1:-3);
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
  djvec_t ave(THp1);
  ave=0.0;
  for(size_t _iConf=0;_iConf<nConfs;_iConf++)
    {
      const size_t iConf=confMap[_iConf];
      const size_t iClust=iConf/clust_size;
      for(size_t iSource=iSourceMin;iSource<iSourceMax;iSource++)
	for(size_t t=0;t<THp1;t++)
	  ave[t][iClust]+=rawData[idData({t,icorr,_iConf,iSource})];
    }
  
  ave.clusterize(clust_size);
  ave/=(iSourceMax-iSourceMin)*L*L*L;
  
  return ave;
}

int main(int narg,char **arg)
{
  raw_file_t input("input.txt","r");
  T=input.read<size_t>("T");
  L=input.read<size_t>("L");
  confsPattern=input.read<string>("ConfsPattern");
  output=input.read<string>("Output");
  tAve=input.read<bool>("tAve");
  
  if(not file_exists(output))
    loadRawData(narg,arg);
  else
    loadData();
  
  raw_file_t confMapFile("map.txt","r");
  for(size_t iConf=0;iConf<nConfs;iConf++)
    confMap.push_back(confMapFile.read<int>());
  
  set_njacks(195);
  
  const djvec_t aveP5=getAve(0,nSources,0);
  const djvec_t aveVK=getAve(0,nSources,1);
  aveP5.ave_err().write("plots/corr_P5P5.xmg");
  aveVK.ave_err().write("plots/corr_VKVK.xmg");
  if(not tAve)
    aveVK.symmetrized().ave_err().write("plots/corr_VKVK_symm.xmg");
  const djack_t mP5=constant_fit(effective_mass(tAve?aveP5:aveP5.symmetrized()),25,32,"plots/eff_mass_P5P5.xmg");
  const djack_t mVK=constant_fit(effective_mass(tAve?aveVK:aveVK.symmetrized()),25,32,"plots/eff_mass_VKVK.xmg");
  const djack_t E2p=2*sqrt(sqr(mP5)+sqr(2*M_PI/L));
  cout<<"mP5: "<<mP5.ave_err()<<endl;
  cout<<"EVK: "<<E2p.ave_err()<<endl;
  
  grace_file_t errVK_plot("plots/VKVK_err.xmg");
  grace_file_t aveVK_plot("plots/VKVK_ave.xmg");
  vec_ave_err_t y1(THp1);
  for(size_t n=1;n<=nSources;n*=2)
    {
      const size_t nCopies=nSources/n;
      
      vec_ave_err_t copyAveVK(THp1);
      vector<double> s(THp1,0.0),s2(THp1,0.0);
      for(size_t iCopy=0;iCopy<nCopies;iCopy++)
	{
	  const djvec_t aveVK=getAve(iCopy*n,(iCopy+1)*n,1);
	  
	  for(size_t t=0;t<THp1;t++)
	    {
 	      copyAveVK[t].ave()+=aveVK[t].ave();
	      copyAveVK[t].err()+=aveVK[t].err();
	      
	      s[t]+=aveVK[t].err();
	      s2[t]+=sqr(aveVK[t].err());
	    }
	}
      
      for(size_t t=0;t<THp1;t++)
	{
	  copyAveVK[t].ave()/=nCopies;
	  copyAveVK[t].err()/=nCopies;
	}
      
      aveVK_plot.write_vec_ave_err(copyAveVK);
      
      if(nCopies>1)
	{
	  vec_ave_err_t y(THp1);
	  for(size_t t=0;t<THp1;t++)
	    {
	      s[t]/=nCopies;
	      s2[t]/=nCopies;
	      s2[t]-=sqr(s[t]);
	      
	      y[t].ave()=s[t];
	      y[t].err()=sqrt(s2[t]/(nCopies-1));
	    }
	  // aveVK[t]*=t*t*t;
	  
	  // normalize
	  // const djvec_t aveTK=getAve(0,nSources,2);
	  // if(n==1)
	  //   y1=y;
	  // else
	  //   {
	  //     for(size_t t=0;t<THp1;t++)
	  // 	{
	  // 	  y[t].err()=sqrt(sqr(y[t].err()/y[t].ave())+sqr(y1[t].err()/y1[t].ave()))*y[t].ave()/y1[t].ave()*sqrt(n);
	  // 	  y[t].ave()*=sqrt(n)/y1[t].ave();
	  // 	}
	    errVK_plot.write_vec_ave_err(y);
	    // }
	}
      // effective_mass(aveTK).ave_err().write("plots/TKTK_new.xmg");
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
