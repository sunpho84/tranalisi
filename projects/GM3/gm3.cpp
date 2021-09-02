#include <iostream>

#include <filesystem>
#include <glob.h>
#include <H5Cpp.h>

#include <mpi.h>

#include <vector>

#include <tranalisi.hpp>

using namespace H5;
using namespace std;

size_t nMPIranks,MPIrank=0;
double amq;
size_t T,L,TH,THp1;
size_t tMinVKVK,tMaxVKVK;
size_t tMinP5P5[2],tMaxP5P5[2];
double a,ZaPetros;
string confsPattern;
string refConfPattern;
string rawDataPackedPath;
size_t nConfs,nSources,nSourcesMax;
double clustSize;

constexpr size_t nGammaComb=7;
const bool isVK[]=                         {0     ,1     ,1     ,1     ,1     ,0     ,0};
const int parity[]=                         {+1    ,+1    ,+1    ,-1    ,-1    ,-1    ,-1};
constexpr char gammaCombTag[nGammaComb][5]={"P5P5","VKVK","TKTK","VKTK","TKVK","A0P5","P5A0"};
enum CORR_ID{idP5P5,idVKVK,idTKTK,idVKTK,idTKVK,idA0P5,idP5A0};

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
  }
  
  return 1.0;
}

constexpr size_t nMes=3;
constexpr char mesTag[nMes][3]={"uu","ud","dd"};

index_t idData;
index_t idData_loader;
index_t idOpenData_loader;
vector<double> _rawData;
bool canUseRawData=false;

double& rawData(const size_t& _iConf,const size_t& iSource,const size_t& igamma_out,const size_t& iMes,const size_t& tOut)
{
  if(not canUseRawData)
    CRASH("Cannot use raw data");
  
  return _rawData[idData({_iConf,iSource,igamma_out,iMes,tOut})];
}

vector<int> confMap;

vector<string> possibleConfsList;
vector<string> sourcesList;

ofstream console;

//! parameters to solve
struct params_t
{
  size_t t;
  double a;
  params_t(size_t t,double a) : t(t),a(a) {}
};

//! compute the kernel f(Q)
double f_Q(double Q,double a)
{
  const double mass_muon=0.1056583745;
  double am=mass_muon*a;
  double w=Q/am;
  double s=sqr(w);
  double A=sqrt(4+s);
  
  return 1/(sqr(am)*w*A)*sqr((A-w)/(A+w));
  //double Z=(sqrt(1+4/s)-1)/2;
  //return 1/sqr(am)*s*Z*Z*Z*(1-s*Z)/(1+s*Z*Z);
}

//! kernel of eq.10
double kern_Q(double Q,void *params)
{
  double t=((params_t*)params)->t;
  double a=((params_t*)params)->a;
  return 4*Q*f_Q(Q,a)*((cos(Q*t)-1)/(Q*Q)+t*t/2);
}

//! parameters of LO corr
struct params_LO_t
{
  double Z2,M,a;
  params_LO_t(double Z2,double M,double a) : Z2(Z2),M(M),a(a) {}
  template <class T> params_LO_t(const vector<T> &p,size_t i) : Z2(p[0][i]),M(p[1][i]),a(p[2][i]) {}
};

//! compute tilde for double
map<pair<size_t,double>,double> looktab;
inline double ftilde_t(size_t t,double a)
{
  pair<size_t,double> key(t,a);
  auto it=looktab.find(key);
  if(it!=looktab.end()) return it->second;
  
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  params_t param(t,a);
  
  //! function structure
  gsl_function f;
  f.function=kern_Q;
  f.params=&param;
  
  //integrate
  double result;
  double abserr;
  double start=0,epsabs=0,epsrel=1e-6;
  gsl_integration_qagiu(&f,start,epsabs,epsrel,workspace_size,workspace,&result,&abserr);
  
  gsl_integration_workspace_free(workspace);
  
  looktab[key]=result;
  
  //cout<<"t: "<<t<<" "<<"ftilde: "<<result<<endl;
  
  return result;
}

//! return the LO
double kern_LO_reco_gsl(double t,void *_params)
{
  params_LO_t *params=(params_LO_t*)_params;
  double &M=params->M;
  double &Z2=params->Z2;
  double &a=params->a;
  
  return ftilde_t(t,a)*nonperiodic_two_pts_corr_fun(Z2,M,t);
}

template <class T> T kern_num(const T &corr_t,double t,const double &a)
{
  return corr_t*ftilde_t(t,a);
}

//! integrate the kernel
template <class TV,class TS=typename TV::base_type> TS integrate_corr_times_kern_up_to(const TV &corr,size_t T,const double &a,size_t &upto,size_t ord=1)
{
  //store the kernel
  TV kern(T/2);
  for(size_t t=1;t<T/2;t++) kern[t]=kern_num(corr[t],t,a);
  kern[0]=0.0;
  
  const double eu=2.0/3,ed=-1.0/3;
  return 4*sqr(alpha_em)*(sqr(eu)+sqr(ed))*integrate_corr_up_to(kern,upto);
}

void setPars()
{
  TH=T/2;
  
  THp1=TH+1;
  
  idData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Gamma",16}});
  idOpenData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Id1",4},{"Id2",4},{"Id3",4},{"Id4",4},{"Ri",2}});
}

void setRawData(const size_t& nConfsToRes)
{
  idData.set_ranges({{"Confs",nConfsToRes},{"Source",nSources},{"GammComb",nGammaComb},{"Mes",nMes},{"T",THp1}});
  _rawData.resize(idData.max(),0.0);
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
  
  console<<"Found: "<<sourcesList.size()<<" sources"<<endl;
  
  return sourcesList;
}

// const Gamma_data_t Gamma_data[nGamma]=
//   {{{0,1,2,3},{1,01,1,1},{0,0,0,0}},
//    {{3,2,1,0},{0,0,0,0},{1,1,-1,-1}},
//    {{3,2,1,0},{1,-1,-1,1},{0,0,0,0}},
//    {{2,3,0,1},{0,0,0,0},{1,-1,-1,1}},
//    {{0,1,2,3},{1,1,-1,-1},{0,0,0,0}},
//    {{2,3,0,1},{1,1,1,1},{0,0,0,0}}};

//The structure for gamma matrix
struct dirac_matr
{
  int pos[4];
  std::complex<double> entr[4];
  
  dirac_matr operator*(const dirac_matr& oth) const
  {
    dirac_matr out;
    
    //This is the line on the first matrix
    for(int ig1=0;ig1<4;ig1++)
      {
	//This is the line to be taken on the second matrix
	const int ig2=this->pos[ig1];
	
	//For each line, the column of the output matrix which is
	//different from 0 is the column of the second matrix different
	//from 0 on the line with index equal to the column of the first
	//matrix which is different from 0 (that is, ig2)
	out.pos[ig1]=oth.pos[ig2];
	
	//The entries of the output is, on each line, the complex
	//product of the entries of the first matrix on that line, for
	//the entries of the second matrix on the line with the index
	//equal to the column of the first matrix which is different
	//from 0 (which again is ig2)
	out.entr[ig1]=this->entr[ig1]*oth.entr[ig2];
	
      }
    return out;
  }
};

dirac_matr gam[16]={dirac_matr{{0,1,2,3},std::complex<double>{1,0},{1,0},{1,0},{1,0}},
		    {{3,2,1,0},std::complex<double>{0,1},{0,1},{0,-1},{0,-1}},
		    {{3,2,1,0},std::complex<double>{1,0},{-1,0},{-1,0},{1,0}},
		    {{2,3,0,1},std::complex<double>{0,1},{0,-1},{0,-1},{0,1}},
		    {{0,1,2,3},std::complex<double>{1,0},{1,0},{-1,0},{-1,0}},
		    {{2,3,0,1},std::complex<double>{1,0},{1,0},{1,0},{1,0}},
		    gam[1]*gam[5],
		    gam[2]*gam[5],
		    gam[3]*gam[5],
		    gam[4]*gam[5],
		    gam[4]*gam[1],
		    gam[4]*gam[2],
		    gam[4]*gam[3],
		    gam[2]*gam[3],
		    gam[3]*gam[1],
		    gam[1]*gam[2]};

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
  
  static constexpr int openRank=7;
  
  hsize_t openDimsm[6];
  DataSpace openMemspace;
  
  hsize_t      openOffset[openRank];   // hyperslab offset in the file
  hsize_t      openCount[openRank];    // size of the hyperslab in the file
  
  std::vector<double> dataIn;
  
  std::vector<double> openDataIn;
  
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
    
    dataIn.resize(idData_loader.max(),0.0);
    
    openDimsm[0]=T;
    openDimsm[1]=4;
    openDimsm[2]=4;
    openDimsm[3]=4;
    openDimsm[4]=4;
    openDimsm[5]=2;
    openMemspace.setExtentSimple(6,openDimsm);
    
    for(int i=0;i<openRank;i++)
      openOffset[i]=0;
    
    openCount[0]=T;
    openCount[1]=1;
    openCount[2]=4;
    openCount[3]=4;
    openCount[4]=4;
    openCount[5]=4;
    openCount[6]=2;
    
    openDataIn.resize(idOpenData_loader.max(),0.0);
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
  
  void load()
  {
    for(size_t iMes=0;iMes<nMes;iMes++)
      {
	DataSet dataset;
	
	const string fullGroupName=groupName+"/mesons/"+mesTag[iMes];
	try
	  {
	    dataset=file.openDataSet(fullGroupName);
	  }
	catch(const H5::FileIException& exc)
	  {
	    CRASH("Unable to open group %s",fullGroupName.c_str());
	  }
	
	const DataSpace dataspace=dataset.getSpace();
	
	dataspace.selectHyperslab(H5S_SELECT_SET,count,offset);
	
	dataset.read(&dataIn[idData_loader({iMes,0,0})],PredType::NATIVE_DOUBLE,memspace,dataspace);
      }
  }
  
  void openLoad()
  {
    const string tag[3]={"uu_open","ud_open","dd_open"};
    
    for(size_t iMes=0;iMes<nMes;iMes++)
      {
	const DataSet dataset=file.openDataSet(groupName+"/mesons/"+tag[iMes]);
	const DataSpace dataspace=dataset.getSpace();
	
	dataspace.selectHyperslab(H5S_SELECT_SET,openCount,openOffset);
	
	dataset.read(&openDataIn[idOpenData_loader({iMes,0,0,0,0,0,0})],PredType::NATIVE_DOUBLE,openMemspace,dataspace);
      }
  }
};

string sourceName(const size_t& iConf,const size_t& iSource)
{
  return possibleConfsList[iConf]+"/"+sourcesList[iSource];
}

void loadAndPackRawData(int narg,char** arg)
{
  console<<"Loading raw data from scratch"<<endl;
  
  possibleConfsList=getConfsList(confsPattern);
  sourcesList=getSourcesList(refConfPattern);
  nSources=sourcesList.size();
  console<<"NSources: "<<nSources<<endl;
  if(nSourcesMax and nSourcesMax<nSources)
    {
      console<<"Reducing nSources from "<<nSources<<" to "<<nSourcesMax<<endl;
      sourcesList.resize(nSourcesMax);
      nSources=nSourcesMax;
    }
  
  const size_t refSize=
    std::filesystem::file_size(sourceName(0,0));
  console<<"Reference size: "<<refSize<<endl;
  
  const size_t nPossibleConfs=possibleConfsList.size();
  vector<int> accepted(nPossibleConfs,true);
  {
    const size_t confChunk=(nPossibleConfs+nMPIranks-1)/nMPIranks;
    const size_t firstConf=std::min((size_t)MPIrank*confChunk,nPossibleConfs);
    const size_t lastConf=std::min(firstConf+confChunk,nPossibleConfs);
    for(size_t iConf=firstConf;iConf<lastConf;iConf++)
      {
	const string& conf=possibleConfsList[iConf];
	
	string sourcePath;
	size_t iSource=0;
	while(iSource<nSources and accepted[iConf])
	  {
	    sourcePath=sourceName(iConf,iSource);
	    
	    const bool exists=
	      std::filesystem::exists(sourcePath);
	    
	    bool correctSize=false;
	    size_t size=refSize;
	    if(exists)
	      {
		size=std::filesystem::file_size(sourcePath);
		correctSize=(size==refSize);
	      }
	    
	    accepted[iConf]=exists and correctSize;
	    
	    if(not accepted[iConf])
	      {
		ostringstream os;
		os<<"Conf "<<conf<<" discarded since "<<sourcePath;
		if(not exists)
		  os<<" does not exists"<<endl;
		else
		  if(not correctSize)
		    os<<" has wrong size "<<size<<" against reference size "<<refSize;
		cout<<os.str()<<endl;
	      }
	    iSource++;
	  }
      }
  }
  
  for(size_t loopRank=0;loopRank<nMPIranks;loopRank++)
    {
      const size_t confChunk=(nPossibleConfs+nMPIranks-1)/nMPIranks;
      const size_t firstConf=std::min((size_t)loopRank*confChunk,nPossibleConfs);
      const size_t lastConf=std::min(firstConf+confChunk,nPossibleConfs);
      MPI_Bcast(&accepted[firstConf],lastConf-firstConf,MPI_INT,loopRank,MPI_COMM_WORLD);
    }
  
  vector<string> confsList;
  for(size_t iConf=0;iConf<nPossibleConfs;iConf++)
    if(accepted[iConf])
      confsList.push_back(possibleConfsList[iConf]);
  
  if(confsList.size()!=possibleConfsList.size())
    console<<"Confs list resized from "<<possibleConfsList.size()<<"to: "<<confsList.size()<<endl;
  
  nConfs=confsList.size();
  console<<"NConfs: "<<nConfs<<endl;
  
  setPars();
  
  DataLoader loader(T);
  const array<pair<int,int>,7> map{std::pair<int,int>{0,0},{1,6},{1,7},{1,8},{2,10},{2,11},{2,12}};
  
  const size_t confChunkSize=(nConfs+nMPIranks-1)/nMPIranks;
  const size_t firstConf=std::min((size_t)MPIrank*confChunkSize,nConfs);
  const size_t lastConf=std::min(firstConf+confChunkSize,nConfs);
  const size_t nConfsPerRank=lastConf-firstConf;
  
  const size_t nConfsToStore=
    (MPIrank!=0)?nConfsPerRank:nConfs;
  
  setRawData(nConfsToStore);
  
  for(size_t _iConf=0;_iConf<nConfsPerRank;_iConf++)
    {
      const size_t iConf=_iConf+firstConf;
      
      cout<<MPIrank<<" "<<iConf<<" ["<<firstConf<<":"<<lastConf<<"] "<<confsList[iConf]<<endl;
      
      for(size_t iSource=0;iSource<nSources;iSource++)
	{
	  const string file=confsList[iConf]+"/"+sourcesList[iSource];
	  loader.open(file);
	  loader.load();
	  loader.openLoad();
	  
	  for(size_t iMes=0;iMes<nMes;iMes++)
	    for(size_t tIn=0;tIn<T;tIn++)
	      for(const auto& m : map)
		{
		  const size_t tOut=
		    ((tIn>=TH)?(T-tIn):tIn);
		  const size_t& igamma_out=m.first;
		  const size_t& igamma_in=m.second;
		  const double& in=loader.dataIn[idData_loader({iMes,tIn,igamma_in})];
		  double& out=rawData(_iConf,iSource,igamma_out,iMes,tOut);
		  
		  out+=in;
		}
	  
	  // A(i)=(GSO)_{ij(i)} (G5)_{j(i)}
	  // B(k)=(G5)_k (GSI)_{kl(k)}
	  const array<array<int,4>,8> map{array<int,4>
					  {idVKTK,10,1,-1},{idVKTK,11,2,-1},{idVKTK,12,3,-1},
					  {idTKVK,1,10,-1},{idTKVK,2,11,-1},{idTKVK,3,12,-1},
					  {idA0P5,9,5,-1},
					  {idP5A0,5,9,-1}};
	  for(size_t iMes=0;iMes<nMes;iMes++)
	    for(size_t tIn=0;tIn<T;tIn++)
	      {
		const size_t tOut=
		  ((tIn>=TH)?(T-tIn):tIn);
		
		for(const auto& m : map)
		  {
		    const double s=
		      ((tIn>=TH)?m[3]:+1.0);
		    
		    const size_t iGammaOut=m[0];
		    
		    const size_t iGammaIn1=m[1];
		    const size_t iGammaIn2=m[2];
		    
		    const dirac_matr g1=gam[iGammaIn1]*gam[5];
		    const dirac_matr g2=gam[5]*gam[iGammaIn2];
		    
		    double& out=rawData(_iConf,iSource,iGammaOut,iMes,tOut);
		    
		    for(size_t nu=0;nu<4;nu++)
		      for(size_t rh=0;rh<4;rh++)
			{
			  const size_t si=g1.pos[nu];
			  const size_t mu=g2.pos[rh];
			  const auto& c1=g1.entr[nu];
			  const auto& c2=g2.entr[rh];
			  
			  const complex<double> c=c1*c2;
			  
			  for(size_t ri=0;ri<2;ri++)
			    {
			      const double& in=loader.openDataIn[idOpenData_loader({iMes,tIn,mu,nu,si,rh,ri})];
			      
			      out+=in*((double*)&c)[ri]*s;
			    }
			}
		  }
	      }
	}
    }
  
  if(MPIrank!=0)
    {
      if(nConfsPerRank)
	MPI_Send(&_rawData[0],idData.max(),MPI_DOUBLE,0,MPIrank,MPI_COMM_WORLD);
    }
  else
    {
      for(size_t iRank=1;iRank<nMPIranks;iRank++)
	{
	  const size_t firstConf=std::min(iRank*confChunkSize,nConfs);
	  const size_t lastConf=std::min(firstConf+confChunkSize,nConfs);
	  const size_t nConfsToRecv=lastConf-firstConf;
	  
	  if(nConfsToRecv>0)
	    {
	      const size_t beg=idData({firstConf,0,0,0,0}),size=idData({nConfsToRecv,0,0,0,0});
	      cout<<"Receiving from rank "<<iRank<<" from conf "<<firstConf<<" for a block of confs "<<nConfsToRecv<<endl;
	      MPI_Recv(&_rawData[beg],size,MPI_DOUBLE,iRank,iRank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }
	}
      
      for(size_t i=0;i<idData.max();i++)
	{
	  const std::vector<size_t> coords=idData(i);
	  
	  const size_t &t=coords.back();
	  const size_t &ig=coords[2];
	  
	  const double norm=((t==0 or t==TH)?1:2)*
	    corrNorm(CORR_ID(ig));
	  
	  _rawData[i]/=norm;
	}
      
      raw_file_t out(rawDataPackedPath,"w");
      out.bin_write(nConfs);
      out.bin_write(nSources);
      out.bin_write(_rawData);
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
}

const string cachedAveCorrPath="cachedAveCorr.dat";
using AveId=std::array<size_t,4>;
std::map<AveId,djvec_t> aveCorrCache;
bool hasTostoreCachedAveCorr=false;;

bool loadCachedAveCorr()
{
  if(not file_exists(cachedAveCorrPath))
    {
      console<<cachedAveCorrPath<<" does not exist, skipping loading cached average correlators"<<endl;
      
      return false;
    }
  
  if(file_exists(rawDataPackedPath) and std::filesystem::last_write_time(cachedAveCorrPath)<std::filesystem::last_write_time(rawDataPackedPath))
    {
      console<<rawDataPackedPath<<" newer than "<<cachedAveCorrPath<<", skipping loading averaged"<<endl;
      
      return false;
    }
  
  raw_file_t file(cachedAveCorrPath,"r");
  console<<"Loading cached average correlators from file "<<cachedAveCorrPath<<endl;
  const size_t cachedNJacks=file.bin_read<size_t>();
  
  //load only if njacks agree
  if(cachedNJacks!=njacks)
    {
      console<<"NJacks in the file is "<<cachedNJacks<<" expecting "<<njacks<<", not loading"<<endl;
      
      return false;
    }
  
  file.bin_read(nConfs);
  file.bin_read(nSources);
  if(nSourcesMax and nSources>nSourcesMax)
    {
      console<<"NSources in the file is "<<nSources<<" is exceeding nSourcesMax, "<<nSourcesMax<<", not loading"<<endl;
      
      return false;
    }
  
  const size_t nCached=file.bin_read<size_t>();
  console<<"Reading  "<<nCached<<" cached average correlators"<<endl;
  
  for(size_t i=0;i<nCached;i++)
    {
      const AveId id=file.bin_read<AveId>();
      
      const size_t n=file.bin_read<size_t>();
      djvec_t data(n);
      file.bin_read(data);
      aveCorrCache[id]=data;
    }
  
  return true;
}

void storeCachedAveCorr()
{
  const size_t& nCached=aveCorrCache.size();
  
  console<<"Storing "<<nCached<<" cached average correlators to file "<<cachedAveCorrPath<<endl;
  
  raw_file_t file(cachedAveCorrPath,"w");
  file.bin_write(njacks);
  file.bin_write(nConfs);
  file.bin_write(nSources);
  file.bin_write(nCached);
  
  for(const auto& c : aveCorrCache)
    {
      file.bin_write(c.first);
      file.bin_write(c.second.size());
      file.bin_write(c.second);
    }
}

void loadRawData(int narg,char **arg)
{
  // Can use raw data
  canUseRawData=true;
  
  if(not file_exists(rawDataPackedPath))
    {
      console<<rawDataPackedPath<<" not found, computing everything from scratch"<<endl;
      loadAndPackRawData(narg,arg);
    }
  console<<"Reading packed raw data"<<endl;
  
  raw_file_t out(rawDataPackedPath,"r");
  out.bin_read(nConfs);
  console<<"NConfs: "<<nConfs<<endl;
  out.bin_read(nSources);
  console<<"NSorces: "<<nSources<<endl;
  
  setRawData(nConfs);
  console<<"Data size: "<<_rawData.size()<<endl;
  
  out.bin_read(_rawData);
}

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes)
{
  const AveId id{iSourceMin,iSourceMax,iGammaComb,iMes};
  
  const auto ref=aveCorrCache.find(id);
  if(ref!=aveCorrCache.end())
    return ref->second;
  else
    if(not canUseRawData)
      CRASH(" cannot reconstruct average as rawdata not present!");
  
  djvec_t ave(THp1);
  ave=0.0;
  jackknivesFill(nConfs,[iSourceMin,iSourceMax,iGammaComb,iMes,&ave](const size_t& iConf,const size_t& iClust,const double& w)
  {
    for(size_t iSource=iSourceMin;iSource<iSourceMax;iSource++)
      for(size_t t=0;t<THp1;t++)
	ave[t][iClust]+=rawData(iConf,iSource,iGammaComb,iMes,t)*w;
  });
  ave.clusterize(clustSize);
  ave/=(iSourceMax-iSourceMin)*L*L*L;
  
  if(iSourceMin==0 and iSourceMax==nSources)
    {
      aveCorrCache[id]=ave;
      hasTostoreCachedAveCorr=true;
    }
  
  return ave;
}

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb)
{
  djvec_t ave(THp1);
  for(size_t iMes=0;iMes<nMes;iMes+=2)
    ave+=getAve(iSourceMin,iSourceMax,iGammaComb,iMes);
  
  return ave/2;
}

enum RegoType{REGO_TM,REGO_OS};

djvec_t getAveForRego(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const RegoType& rego)
{
  djvec_t res;
  
  switch (rego)
    {
    case REGO_TM:
      res=getAve(iSourceMin,iSourceMax,iGammaComb);
      break;
    case REGO_OS:
      res=getAve(iSourceMin,iSourceMax,iGammaComb,1);
    }
  
  return res;
}

void rawDataAn(const size_t& iGammaComb)
{
  const string tag=gammaCombTag[iGammaComb];
  const size_t nCopies=2;
  const size_t nSourcesToBeUsed=nSources;
  const size_t nSourcesPerCopy=nSourcesToBeUsed/nCopies;
  // index_t id({{"copy",ncopies},{"conf",nConfs},{"T",T/2+1}});
  
  vector<double> sourceAveData;
  index_t idDataAve({{"Confs",nConfs},{"Copy",nCopies},{"T",THp1}});
  sourceAveData.resize(idDataAve.max(),0.0);
  for(size_t t=0;t<=T/2;t++)
    for(size_t iConf=0;iConf<nConfs;iConf++)
      for(size_t iCopy=0;iCopy<nCopies;iCopy++)
	{
	  double& s=sourceAveData[idDataAve({iConf,iCopy,t})];
	  s=0.0;
	  for(size_t iSource=nSourcesPerCopy*iCopy;iSource<nSourcesPerCopy*(iCopy+1);iSource++)
	    for(size_t iMes=0;iMes<nMes;iMes+=2)
	      s+=rawData(iConf,iSource,iGammaComb,iMes,t);
	  s/=nSourcesPerCopy*L*L*L*2;
	}
  
  grace_file_t copyAve("plots/copy_ave_"+tag+".xmg");
  vector<djvec_t> copyAveData(2,djvec_t(THp1));
  for(size_t iCopy=0;iCopy<nCopies;iCopy++)
    {
      djvec_t& c=copyAveData[iCopy];
      for(size_t t=0;t<=T/2;t++)
	jackknivesFill(nConfs,[&c,&sourceAveData,&idDataAve,t,iCopy](const size_t& iConf,const size_t& iClust,const double& w)
	{
	    c[t][iClust]+=sourceAveData[idDataAve({iConf,iCopy,t})]*w;
	});
      c.clusterize(clustSize);
      c/=nSourcesPerCopy*L*L*L*2;
      copyAve.write_vec_ave_err(c.ave_err());
    }
  
  // vector<double> covData(THp1);
  // for(size_t t=0;t<=TH;t++)
  //   {
  //     covData[t]=corr(copyAveData[0][t],copyAveData[1][t]);
  //     cout<<"covData["<<t<<"] :"<<covData[t]<<endl;
  //   }
  
  djvec_t c((THp1)*nCopies*nCopies),a((THp1)*nCopies),aa((THp1)),cc(THp1);
  jackknivesFill(nConfs,
		 [&sourceAveData,&idDataAve,&a,&c,&aa,&cc](const size_t& iConf,const size_t& iClust,const double& w)
		 {
		   for(size_t t=0;t<=T/2;t++)
		     {
		       double copy_ave=0;
		       for(size_t iCopy=0;iCopy<nCopies;iCopy++)
			 copy_ave+=sourceAveData[idDataAve({iConf,iCopy,t})];
		       copy_ave/=nCopies;
		       
		       aa[t][iClust]+=copy_ave*w;
		       cc[t][iClust]+=copy_ave*copy_ave*w;
		       
		       for(size_t iCopy=0;iCopy<nCopies;iCopy++)
			 {
			   const double& x=sourceAveData[idDataAve({iConf,iCopy,t})];
			   a[iCopy+nCopies*t][iClust]+=x*w;
			   
			   for(size_t jCopy=0;jCopy<nCopies;jCopy++)
			     {
			       const double& y=sourceAveData[idDataAve({iConf,jCopy,t})];
			       c[jCopy+nCopies*(iCopy+nCopies*t)][iClust]+=x*y*w;
			     }
			 }
		     }
		 });
  
  for(auto& x : {&a,&c,&aa,&cc})
    x->clusterize(clustSize);
  
  djvec_t ave(THp1),err(THp1),corr(THp1);
  for(size_t t=0;t<=T/2;t++)
    {
      cc[t]-=aa[t]*aa[t];
      cc[t]=sqrt(cc[t]/(nConfs-1));
      
      for(size_t icopy=0;icopy<nCopies;icopy++)
	for(size_t jcopy=0;jcopy<nCopies;jcopy++)
	  c[jcopy+nCopies*(icopy+nCopies*t)]-=a[icopy+nCopies*t]*a[jcopy+nCopies*t];
      
      ave[t]+=a[0+nCopies*t];
      err[t]+=c[0+nCopies*(0+nCopies*t)];
      
      err[t]=sqrt(err[t]/(nConfs-1));
      corr[t]=0;
      for(size_t icopy=0;icopy<nCopies;icopy++)
	for(size_t jcopy=icopy+1;jcopy<nCopies;jcopy++)
	  corr[t]+=c[icopy+nCopies*(jcopy+nCopies*t)]/sqrt(c[icopy+nCopies*(icopy+nCopies*t)]*c[jcopy+nCopies*(jcopy+nCopies*t)]);
      corr[t]/=nCopies*(nCopies-1)/2.0;
    }
  
  // ave.ave_err().write(combine("plots/%s_ave.xmg",tag));
  // aa.ave_err().write(combine("plots/%s_combo_ave.xmg",tag));
  // cc.ave_err().write(combine("plots/%s_combo_err.xmg",tag));
  err.ave_err().write("plots/err_with_err_"+tag+".xmg");
  corr.ave_err().write("plots/stat_corr_"+tag+".xmg");
}

void readInput()
{
  raw_file_t input("input.txt","r");
  T=input.read<size_t>("T");
  L=input.read<size_t>("L");
  amq=input.read<double>("amq");
  a=input.read<double>("a");
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

djvec_t determineRenoConst()
{
  djvec_t Z(2),Zsilv(2);
  
  djack_t mP[2],ZA0[2],ZP5[2];
  const djvec_t corrP5P5[2]=
    {(getAve(0,nSources,idP5P5,0)+
      getAve(0,nSources,idP5P5,2))/2.0,
      getAve(0,nSources,idP5P5,1)};
  // const djvec_t corrA0P5[2]=
  //   {(getAve(0,nSources,idA0P5,0)+
  //     getAve(0,nSources,idA0P5,2))/2.0,
  //     getAve(0,nSources,idA0P5,1)};
  const djvec_t corrP5A0[2]=
    {(-getAve(0,nSources,idP5A0,0)+
      -getAve(0,nSources,idP5A0,2))/2.0,
      -getAve(0,nSources,idP5A0,1)};
  
  for(size_t iMes=0;iMes<2;iMes++)
    {
      djack_t Z2P5separated,MP5separated;
      two_pts_fit(Z2P5separated,MP5separated,corrP5P5[iMes],TH,tMinP5P5[iMes],tMaxP5P5[iMes],"plots/fitP5P5forZmes"+to_string(iMes)+".xmg");
      djvec_t testA0=corrP5A0[iMes];
      const djack_t ZP5separated=sqrt(Z2P5separated);
      for(size_t t=0;t<=TH;t++)
	testA0[t]/=two_pts_corr_fun(ZP5separated,MP5separated,TH,t,-1);
      const djack_t ZA0separated=constant_fit(testA0,tMinP5P5[iMes],tMaxP5P5[iMes],"plots/fitP5A0forZmes"+to_string(iMes)+".xmg");
      console<<"from separated fit, mP: "<<MP5separated.ave_err()<<" , ZA0: "<<ZA0separated.ave_err()<<" , ZP5: "<<ZP5separated.ave_err()<<endl;
      
      two_pts_SL_fit(ZP5[iMes],ZA0[iMes],mP[iMes],corrP5A0[iMes],corrP5P5[iMes],TH,tMinP5P5[iMes],tMaxP5P5[iMes],combine("plots/A0P5FitMes%zu.xmg",iMes),-1,+1);
      console<<"from combined fit, mP: "<<mP[iMes].ave_err()<<" , ZA0: "<<ZA0[iMes].ave_err()<<" , ZP5: "<<ZP5[iMes].ave_err()<<endl;
      
      const djack_t fPfromP=2.0*ZP5[0]*amq/sqr(mP[0]);
      const djack_t fPfromA=ZA0[iMes]/mP[iMes];
      Z[iMes]=fPfromP/fPfromA;
      
      console<<"Z"<<((iMes==0)?"V":"A")<<": "<<Z[iMes].ave_err()<<endl;
    }
  
  const djack_t ZV_fr_ZA=Z[0]/Z[1];
  console<<"Zv/Za: "<<ZV_fr_ZA.ave_err()<<endl;

  const djvec_t derP5A0TM=-symmetric_derivative(corrP5A0[0]);
  derP5A0TM.ave_err().write("plots/derP5A0_regoTM.xmg");
  const djvec_t ZvSilvCorr=2*amq*corrP5P5[0]/derP5A0TM;
  const djack_t ZvSilv=constant_fit(ZvSilvCorr,18,54,"plots/ZvSilv.xmg");
  console<<"ZvSilv: "<<ZvSilv.ave_err()<<endl;
  Zsilv[0]=ZvSilv;
  
  const djack_t ZaSilvCorrectingFactor=ZP5[0]/ZP5[1]*mP[1]*sinh(mP[1])/(mP[0]*sinh(mP[0]));
  console<<"ZaSilvCorrectingFactor: "<<ZaSilvCorrectingFactor.ave_err()<<endl;
  
  const djvec_t derP5A0OS=-symmetric_derivative(corrP5A0[1]);
  derP5A0OS.ave_err().write("plots/derP5A0_regoOS.xmg");
  const djvec_t ZaSilvCorr=2*amq*corrP5P5[1]/derP5A0OS*ZaSilvCorrectingFactor;
  const djack_t ZaSilv=constant_fit(ZaSilvCorr,18,54,"plots/ZaSilv.xmg");
  console<<"ZaSilv: "<<ZaSilv.ave_err()<<endl;
  Zsilv[1]=ZaSilv;
  
  return Zsilv;
}

void convertForSilvano()
{
  for(size_t iConf=0;iConf<nConfs;iConf++)
    {
      const string confPath=combine("reordered/%04zu/",iConf+1);
      const string filePath=combine("reordered/%04zu/mes_contr_2pts_ll",iConf+1);
      
      if(not file_exists(filePath))
	{
	  mkdir(confPath);
	  
	  ofstream out(filePath);
	  out.precision(16);
	  const size_t mapMes[4]={0,1,1,2};
	  for(size_t _iMes=0;_iMes<4;_iMes++)
	    {
	      const size_t iMes=mapMes[_iMes];
	      out<<endl;
	      out<<"  # Contraction of S0_th0_m0_r"<<(_iMes&1)<<"_ll ^ \\dag and S0_th0_m0_r"<<((_iMes>>1)&1)<<"_ll:"<<endl;
	      out<<endl;
	      
	      const size_t mapCorr[6]={0,1,1,1,5,6};
	      const char corrName[6][5]={"P5P5","V1V1","V2V2","V3V3","A0P5","P5A0"};
	      for(size_t _iCorr=0;_iCorr<6;_iCorr++)
		{
		  out<<endl<<"  # "<<corrName[_iCorr]<<endl;
		  for(size_t _t=0;_t<T;_t++)
		    {
		      size_t t=_t;
		      if(t>=TH)
			t=T-t;
		      
		      double s=0;
		      for(size_t iSource=0;iSource<nSources;iSource++)
			s+=rawData(iConf,iSource,mapCorr[_iCorr],iMes,t);
		      s/=nSources*L*L*L;
		      out<<s<<" "<<0<<endl;
		    }
		}
	    }
	}
    }
}

void computeAmu(const djvec_t& Z)
{
  const string regoTag[2]={"TM","OS"};
  const size_t regoZId[2]={1,0};
  
  for(const RegoType& rego : {REGO_TM,REGO_OS})
    {
      console<<"Computing amu fore rego "<<regoTag[rego]<<endl;
      
      const djack_t Zrego=Z[regoZId[rego]];
      
      vector<djvec_t> eig;
      vector<djvec_t> recastEigvec;
      vector<djvec_t> origEigvec;
      
      const djvec_t c00=getAveForRego(0,nSources,1,rego);
      const djvec_t c01=-getAveForRego(0,nSources,4,rego);
      const djvec_t c11=-getAveForRego(0,nSources,2,rego);
      const size_t t0=3;
      
      tie(eig,recastEigvec,origEigvec)=gevp({c00,c01,c01,c11},t0);
      
      eig[0].ave_err().write("plots/eig1Rego"+regoTag[rego]+".xmg");
      eig[1].ave_err().write("plots/eig2Reno"+regoTag[rego]+".xmg");
      
      // const djack_t eig0MDiagFit=constant_fit(effective_mass(eig[0]),tMinFit,tMaxFit,"plots/eff_eig1.xmg");
      // console<<"eig0 mass: "<<eig0MDiagFit.ave_err()<<endl;
      
      // const djack_t eig1MDiagFit=constant_fit(effective_mass(eig[1]),13,18,"plots/eff_eig2.xmg");
      // out<<"eig1 mass: "<<eig1MDiagFit.ave_err()<<endl;
      
      djvec_t SL0(THp1),SS0(THp1);
      djvec_t SL1(THp1),SS1(THp1);
      for(size_t t=0;t<=TH;t++)
	{
	  for(size_t ijack=0;ijack<=njacks;ijack++)
	    {
	      typedef Matrix<double,2,2> Matr;
	      
	      Matr e;
	      const auto& ei=origEigvec;
	      e(0,0)=ei[0][t][ijack];
	      e(0,1)=ei[2][t][ijack];
	      e(1,0)=ei[1][t][ijack];
	      e(1,1)=ei[3][t][ijack];
	      
	      Matr c;
	      c(0,0)=c00[t][ijack];
	      c(0,1)=c01[t][ijack];
	      c(1,0)=c01[t][ijack];
	      c(1,1)=c11[t][ijack];
	      
	      const Matr sl=c*e.transpose();
	      const Matr ss=e*c*e.transpose();
	      
	      SL0[t][ijack]=sl(0,0);
	      SS0[t][ijack]=ss(0,0);
	      SL1[t][ijack]=sl(0,1);
	      SS1[t][ijack]=ss(1,1);
	    }
	}
      
      djack_t eig0ZS,eig0ZL,eig0M;
      two_pts_SL_fit(eig0ZS,eig0ZL,eig0M,SL0,SS0,TH,tMinVKVK,tMaxVKVK,"plots/SL0Rego"+regoTag[rego]+".xmg");
      const djack_t eig0Z2L=eig0ZL*eig0ZL;
      console<<"Z20L: "<<eig0Z2L<<endl;
      
      djvec_t subCorr1=c00;
      for(size_t t=0;t<=TH;t++)
	subCorr1[t]-=two_pts_corr_fun(eig0Z2L,eig0M,TH,t,+1);
      
      djack_t eig1Z2L,eig1M;
      two_pts_fit(eig1Z2L,eig1M,subCorr1,TH,8,19,"plots/SL1Rego"+regoTag[rego]+".xmg");
      console<<"Z21L: "<<eig1Z2L<<endl;
      console<<"M1L: "<<eig1M<<endl;
      
      const djvec_t corr=getAveForRego(0,nSources,1,rego);
      djack_t mVK1,Z2VK1;
      two_pts_fit(Z2VK1,mVK1,corr,TH,tMinVKVK,tMaxVKVK,"plots/eff_mass_VKVK_twopts_fit_rego"+regoTag[rego]+".xmg");
      console<<"Z2: "<<Z2VK1<<endl;
      // djack_t mVK2,Z2VK2;
      // two_pts_fit(Z2VK2,mVK2,corr,TH,22,32);
      
      grace_file_t amu("plots/amuRego"+regoTag[rego]+".xmg");
      djvec_t amuInt(TH),amuSubs1(TH),amuSubs2(TH);
      for(size_t upto=0;upto<TH;upto++)
	{
	  djvec_t corrRefatta1=corr;
	  djvec_t corrRefatta2=corr;
	  for(size_t t=upto;t<TH;t++)
	    {
	      corrRefatta1[t]=two_pts_corr_fun(Z2VK1,mVK1,TH,t,0);
	      corrRefatta2[t]=two_pts_corr_fun(eig0Z2L,eig0M,TH,t,0);
	      corrRefatta2[t]+=two_pts_corr_fun(eig1Z2L,eig1M,TH,t,0);
	    }
	  size_t THm1=TH-1;
	  const djack_t cInt=integrate_corr_times_kern_up_to(corr,T,a,upto)*sqr(Zrego)*1e10;
	  const djack_t cSubs1=integrate_corr_times_kern_up_to(corrRefatta1,T,a,THm1)*sqr(Zrego)*1e10;
	  const djack_t cSubs2=integrate_corr_times_kern_up_to(corrRefatta2,T,a,THm1)*sqr(Zrego)*1e10;
	  amuInt[upto]=cInt;
	  amuSubs1[upto]=cSubs1;
	  amuSubs2[upto]=cSubs2;
	}
      amu.set_xaxis_label("t");
      
      djvec_t corrRefatta(THp1);
      for(size_t t=0;t<TH;t++)
	corrRefatta[t]=
	  two_pts_corr_fun(eig0Z2L,eig0M,TH,t,0)+
	  two_pts_corr_fun(eig1Z2L,eig1M,TH,t,0);
      corrRefatta.ave_err().write("plots/corr_VKVK_refatta"+regoTag[rego]+".xmg");
      effective_mass(corrRefatta,TH,0).ave_err().write("plots/eff_mass_VKVK_refatta_Rego"+regoTag[rego]+".xmg");
      
      amu.write_vec_ave_err(amuInt.ave_err());
      amu.set_no_line();
      amu.set_legend("Pure integration");
      amu.set_all_colors(grace::RED);
      
      amu.write_vec_ave_err(amuSubs1.ave_err());
      amu.set_no_line();
      amu.set_all_colors(grace::BLUE);
      amu.set_legend("Analytic continuation of ground state");
      
      amu.write_vec_ave_err(amuSubs2.ave_err());
      amu.set_no_line();
      amu.set_all_colors(grace::ORANGE);
      amu.set_legend("Ground and first exc.state (GEVP)");
      
      amu.new_data_set();
      amu.set_legend("BMW light connected");
      amu.set_all_colors(grace::GREEN4);
      amu.write_constant_band(0,TH,djack_t(gauss_filler_t{633.7,5.0,23423}));
      
      djvec_t corrRefatta1=corr;
      djvec_t corrRefatta2=corr;
      for(size_t t=0;t<=TH;t++)
	{
	  corrRefatta1[t]=two_pts_corr_fun(Z2VK1,mVK1,TH,t,0);
	  corrRefatta2[t]=two_pts_corr_fun(eig0Z2L,eig0M,TH,t,0);
	}
      corrRefatta1.ave_err().write("plots/corr_refatta1_rego"+regoTag[rego]+".xmg");
      corrRefatta2.ave_err().write("plots/corr_refatta2_rego"+regoTag[rego]+".xmg");
    }
}

void analyzeRawData()
{
  if(not canUseRawData)
    CRASH("cannot do the raw data analysis!");
  
  djack_t mP5;
  for(size_t iGammaComb=0;iGammaComb<nGammaComb;iGammaComb++)
    {
      const string cTag=gammaCombTag[iGammaComb];
      
      djvec_t aveCorr(THp1,0.0);
      for(size_t iMes=0;iMes<nMes;iMes++)
	{
	  const djvec_t corr=getAve(0,nSources,iGammaComb,iMes);
	  corr.ave_err().write("plots/corr_"+cTag+"_"+(string)mesTag[iMes]+".xmg");
	  if(iMes!=1)
	    aveCorr+=corr;
	}
      aveCorr/=2.0;
      aveCorr.ave_err().write("plots/corr_"+cTag+".xmg");
      
      const size_t iMes=0;
      const size_t tMinFit[2]={tMinP5P5[iMes],tMinVKVK};
      const size_t tMaxFit[2]={tMaxP5P5[iMes],tMaxVKVK};
      const size_t is=isVK[iGammaComb];
      const djack_t m=constant_fit(effective_mass(aveCorr,TH,parity[iGammaComb]),tMinFit[is],tMaxFit[is],"plots/eff_mass_"+cTag+".xmg");
      
      switch(iGammaComb)
	{
	case 0:
	  mP5=m;
	  console<<"amP5: "<<m.ave_err()<<endl;
	  console<<"mP5: "<<djack_t(m/a).ave_err()<<endl;
	  break;
	case 1:
	  {
	    const djack_t E2p=2*sqrt(sqr(mP5)+sqr(2*M_PI/L));
	    console<<"EVK: "<<E2p.ave_err()<<endl;
	    console<<"mVK: "<<m.ave_err()<<endl;
	  }
	  break;
	case 2:
	  console<<"mTK: "<<m.ave_err()<<endl;
	  break;
	}
      
      rawDataAn(iGammaComb);
      
      grace_file_t err_plot("plots/err_scaling_"+cTag+".xmg");
      grace_file_t ave_plot("plots/ave_"+cTag+".xmg");
      vector<double> y1(THp1);
      for(size_t n=1;n<=nSources;n*=2)
	{
	  const size_t nCopies=nSources/n;
	  
	  vec_ave_err_t copyAve(THp1);
	  vector<double> s(THp1,0.0),s2(THp1,0.0);
	  for(size_t iCopy=0;iCopy<nCopies;iCopy++)
	    {
	      const djvec_t ave=
		(getAve(iCopy*n,(iCopy+1)*n,iGammaComb,0)+
		 getAve(iCopy*n,(iCopy+1)*n,iGammaComb,2))/2.0;
	      
	      for(size_t t=0;t<THp1;t++)
		{
		  copyAve[t].ave()+=ave[t].ave();
		  copyAve[t].err()+=ave[t].err();
		  
		  s[t]+=ave[t].err();
		  s2[t]+=sqr(ave[t].err());
		}
	    }
	  
	  for(size_t t=0;t<THp1;t++)
	    {
	      copyAve[t].ave()/=nCopies;
	      copyAve[t].err()/=nCopies;
	    }
	  
	  ave_plot.write_vec_ave_err(copyAve);
	  
	  vec_ave_err_t y(THp1);
	  for(size_t t=0;t<THp1;t++)
	    {
	      s[t]/=nCopies;
	      s2[t]/=nCopies;
	      s2[t]-=sqr(s[t]);
	      
	      const double norm=sqrt(n)*exp(mP5.ave()*t);
	      y[t].ave()=s[t]*norm;
	      if(n==1)
		y1[t]=y[t].ave();
	      
	      y[t].err()=sqrt(s2[t]*norm/std::max(1lu,nCopies-1));
	      
	      y[t].ave()/=y1[t];
	      y[t].err()/=y1[t];
	    }
	  
	  err_plot.write_vec_ave_err(y);
	  err_plot.set_legend(to_string(n)+"hits");
	}
    }
}

int main(int narg,char **arg)
{
  MPI_Init(&narg,&arg);
  {
    int temp;
    MPI_Comm_size(MPI_COMM_WORLD,&temp);
    nMPIranks=temp;
    MPI_Comm_rank(MPI_COMM_WORLD,&temp);
    MPIrank=temp;
  }
  console.open((MPIrank==0)?"/dev/stdout":"/dev/null");
  
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
  
  if(hasTostoreCachedAveCorr)
    storeCachedAveCorr();
  
  MPI_Finalize();
  
  return 0;
}
