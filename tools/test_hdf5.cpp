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
size_t tMinFit,tMaxFit;
double a,ZaPetros;
djack_t Za,Zv;
string confsPattern;
string refConfPattern;
string output;
size_t nConfs,nSources;
constexpr size_t nGammaComb=7;
constexpr char gammaCombTag[nGammaComb][5]={"P5P5","VKVK","TKTK","VKTK","TKVK","A0P5","P5A0"};
enum{idP5P5,idVKVK,idTKTK,idVKTK,idTKVK,idA0P5,idP5A0};
constexpr size_t nMes=3;
constexpr char mesTag[nMes][3]={"uu","ud","dd"};

index_t idData;
index_t idData_loader;
index_t idOpenData_loader;
vector<double> rawData;
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
  
  console<<"NConfs: "<<nConfs<<endl;
  console<<"NSources: "<<nSources<<endl;
}

void setRawData(const size_t& nConfsToRes)
{
  idData.set_ranges({{"Confs",nConfsToRes},{"Source",nSources},{"GammComb",nGammaComb},{"Mes",nMes},{"T",THp1}});
  rawData.resize(idData.max(),0.0);
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
	const DataSet dataset=file.openDataSet(groupName+"/mesons/"+mesTag[iMes]);
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

void loadRawData(int narg,char** arg)
{
  possibleConfsList=getConfsList(confsPattern);
  sourcesList=getSourcesList(refConfPattern);
  nSources=sourcesList.size();
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
		  double& out=rawData[idData({_iConf,iSource,igamma_out,iMes,tOut})];
		  
		  out+=in;
		}
	  
	  // A(i)=(GSO)_{ij(i)} (G5)_{j(i)}
	  // B(k)=(G5)_k (GSI)_{kl(k)}
	  const array<array<int,4>,8> map{array<int,4>{3,10,1,-1},{3,11,2,-1},{3,12,3,-1},{4,1,10,-1},{4,2,11,-1},{4,3,12,-1},{5,9,5,-1},{6,5,9,-1}};
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
		    
		    double& out=rawData[idData({_iConf,iSource,iGammaOut,iMes,tOut})];
		    
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
	MPI_Send(&rawData[0],idData.max(),MPI_DOUBLE,0,MPIrank,MPI_COMM_WORLD);
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
	      MPI_Recv(&rawData[beg],size,MPI_DOUBLE,iRank,iRank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }
	}
      
      for(size_t i=0;i<idData.max();i++)
	{
	  const std::vector<size_t> coords=idData(i);
	  
	  const size_t &t=coords.back();
	  const size_t &ig=coords[2];
	  
	  const double norm=((t==0 or t==TH)?1:2)*
	    ((ig==0 or ig==6 or ig==7)?1:-3);
	  rawData[i]/=norm;
	}
      
      raw_file_t out(output,"w");
      out.bin_write(nConfs);
      out.bin_write(nSources);
      out.bin_write(rawData);
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
}

void loadData()
{
  raw_file_t out(output,"r");
  out.bin_read(nConfs);
  out.bin_read(nSources);
  
  setPars();
  setRawData(nConfs);
  console<<"Data size: "<<rawData.size()<<endl;
  
  out.bin_read(rawData);
}

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes)
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
	  ave[t][iClust]+=rawData[idData({_iConf,iSource,iGammaComb,iMes,t})];
    }
  
  ave.clusterize(clust_size);
  ave/=(iSourceMax-iSourceMin)*L*L*L;
  
  return ave;
}

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb)
{
  djvec_t ave(THp1);
  for(size_t iMes=0;iMes<nMes;iMes+=2)
    ave+=getAve(iSourceMin,iSourceMax,iGammaComb,iMes);
  
  return ave/2;
}

void an(const size_t& iGammaComb)
{
  const string tag=gammaCombTag[iGammaComb];
  const size_t clust_size=nConfs/njacks;
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
	      s+=rawData[idData({iConf,iSource,iGammaComb,iMes,t})];
	  s/=nSourcesPerCopy*L*L*L*2;
	}
  
  grace_file_t copyAve("plots/copy_ave_"+tag+".xmg");
  vector<djvec_t> copyAveData(2,djvec_t(THp1));
  for(size_t iCopy=0;iCopy<nCopies;iCopy++)
    {
      djvec_t& c=copyAveData[iCopy];
      for(size_t t=0;t<=T/2;t++)
	for(size_t iConf=0;iConf<nConfs;iConf++)
	  {
	    const size_t iClust=iConf/clust_size;
	    c[t][iClust]+=sourceAveData[idDataAve({iConf,iCopy,t})];
	  }
      c.clusterize(clust_size);
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
  for(size_t t=0;t<=T/2;t++)
    for(size_t ijack=0;ijack<njacks+1;ijack++)
      for(size_t iConf=0;iConf<nConfs;iConf++)
	if(iConf<clust_size*ijack or iConf>=clust_size*(ijack+1))
	  {
	    double copy_ave=0;
	    for(size_t iCopy=0;iCopy<nCopies;iCopy++)
	      copy_ave+=sourceAveData[idDataAve({iConf,iCopy,t})];
	    copy_ave/=nCopies;
	    
	    aa[t][ijack]+=copy_ave;
	    cc[t][ijack]+=copy_ave*copy_ave;
	    
	    for(size_t iCopy=0;iCopy<nCopies;iCopy++)
	      {
		const double& x=sourceAveData[idDataAve({iConf,iCopy,t})];
		a[iCopy+nCopies*t][ijack]+=x;
		
		for(size_t jCopy=0;jCopy<nCopies;jCopy++)
		  {
		    const double& y=sourceAveData[idDataAve({iConf,jCopy,t})];
		    c[jCopy+nCopies*(iCopy+nCopies*t)][ijack]+=x*y;
		  }
	      }
	  }
  
  djvec_t ave(THp1),err(THp1),corr(THp1);
  for(size_t t=0;t<=T/2;t++)
    {
      for(size_t ijack=0;ijack<njacks;ijack++)
	{
	  aa[t][ijack]/=nConfs-clust_size;
	  cc[t][ijack]/=nConfs-clust_size;
	}
      aa[t][njacks]/=nConfs;
      cc[t][njacks]/=nConfs;
      cc[t]-=aa[t]*aa[t];
      cc[t]=sqrt(cc[t]/(nConfs-1));
      
      for(size_t icopy=0;icopy<nCopies;icopy++)
	{
	  for(size_t ijack=0;ijack<njacks;ijack++)
	    a[icopy+nCopies*t][ijack]/=nConfs-clust_size;
	  a[icopy+nCopies*t][njacks]/=nConfs;
	}
      
      for(size_t icopy=0;icopy<nCopies;icopy++)
	for(size_t jcopy=0;jcopy<nCopies;jcopy++)
	  {
	    for(size_t ijack=0;ijack<njacks;ijack++)
	      c[jcopy+nCopies*(icopy+nCopies*t)][ijack]/=nConfs-clust_size;
	    c[jcopy+nCopies*(icopy+nCopies*t)][njacks]/=nConfs;
	  }
      
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
  output=input.read<string>("Output");
  tMinFit=input.read<size_t>("TFit");
  tMaxFit=input.read<size_t>();
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

void determineRenoConst()
{
  djack_t mP[2],ZA0[2],ZP5[2];
  const djvec_t corrP5P5[2]=
    {(getAve(0,nSources,idP5P5,0)+
      getAve(0,nSources,idP5P5,2))/2.0,
      getAve(0,nSources,idP5P5,1)};
  const djvec_t corrA0P5[2]=
    {(-getAve(0,nSources,idA0P5,0)+
      -getAve(0,nSources,idA0P5,2))/2.0,
      -getAve(0,nSources,idA0P5,1)};
  const djvec_t corrP5A0[2]=
    {(-getAve(0,nSources,idP5A0,0)+
      -getAve(0,nSources,idP5A0,2))/2.0,
      -getAve(0,nSources,idP5A0,1)};
      
  for(size_t iMes=0;iMes<2;iMes++)
    {
      two_pts_SL_fit(ZP5[iMes],ZA0[iMes],mP[iMes],corrP5A0[iMes],corrP5P5[iMes],TH,28,36,combine("plots/A0P5FitMes%zu.xmg",iMes),-1,+1);
      console<<"mP: "<<mP[iMes].ave_err()<<" , ZA0: "<<ZA0[iMes].ave_err()<<" , ZP5: "<<ZP5[iMes].ave_err()<<endl;
      
      const djack_t fPfromP=2.0*ZP5[0]*amq/sqr(mP[0]);
      const djack_t fPfromA=ZA0[iMes]/mP[iMes];
      const djack_t Z=fPfromP/fPfromA;
      
      console<<"Z"<<((iMes==0)?"V":"A")<<": "<<Z.ave_err()<<endl;
      if(iMes==0)
	Zv=Z;
      else
	Za=Z;
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
  
  if(not file_exists(output))
    loadRawData(narg,arg);
  
  loadData();
  
  readConfMap();
  
  set_njacks(nConfs);
  
  determineRenoConst();
  
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
      
      const djack_t m=constant_fit(effective_mass(aveCorr),(iGammaComb==0)?22:20,32,"plots/eff_mass_"+cTag+".xmg");
      
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
      
      an(iGammaComb);
      
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
  
  // {
  //   grace_file_t sil("/tmp/sil.xmg");
  //   const djvec_t a120=(getAve(0,120,1,0)+
  // 			getAve(0,120,1,2))/2.0;
  //   const djvec_t a512=(getAve(0,nSources,1,0)+
  // 			getAve(0,nSources,1,2))/2.0;
  
  //   sil.write_vec_ave_err(a120.ave_err());
  //   sil.set_legend("120 hits (simone");
  //   sil.write_vec_ave_err(a512.ave_err());
  //   sil.set_legend("512 hits (simone");
  // }
  
  for(size_t iConf=0;iConf<nConfs;iConf++)
    {
      mkdir(combine("reordered/%04zu/",iConf+1));
      ofstream out(combine("reordered/%04zu/mes_contr_2pts_ll",iConf+1));
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
		    s+=rawData[idData({iConf,iSource,mapCorr[_iCorr],iMes,t})];
		  s/=nSources*L*L*L;
		  out<<s<<" "<<0<<endl;
		}
	    }
	}
    }
  
  
  /////////////////////////////////////////////////////////////////
  
  vector<djvec_t> eig;
  vector<djvec_t> recastEigvec;
  vector<djvec_t> origEigvec;
  
  const djvec_t c00=getAve(0,nSources,1);
  const djvec_t c01=-getAve(0,nSources,4);
  const djvec_t c11=-getAve(0,nSources,2);
  const size_t t0=3;
  
  tie(eig,recastEigvec,origEigvec)=gevp({c00,c01,c01,c11},t0);
  
  eig[0].ave_err().write("plots/eig1.xmg");
  eig[1].ave_err().write("plots/eig2.xmg");
  
  const djack_t eig0MDiagFit=constant_fit(effective_mass(eig[0]),tMinFit,tMaxFit,"plots/eff_eig1.xmg");
  console<<"eig0 mass: "<<eig0MDiagFit.ave_err()<<endl;
  
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
  two_pts_SL_fit(eig0ZS,eig0ZL,eig0M,SL0,SS0,TH,tMinFit,tMaxFit,"plots/SL0.xmg");
  const djack_t eig0Z2L=eig0ZL*eig0ZL;
  console<<"Z20L: "<<eig0Z2L<<endl;
  
  djvec_t subCorr1=c00;
  for(size_t t=0;t<=TH;t++)
      subCorr1[t]-=two_pts_corr_fun(eig0Z2L,eig0M,TH,t,+1);
  
  djack_t eig1Z2L,eig1M;
  two_pts_fit(eig1Z2L,eig1M,subCorr1,TH,8,19,"plots/SL1.xmg");
  console<<"Z21L: "<<eig1Z2L<<endl;
  console<<"M1L: "<<eig1M<<endl;
  
  const djvec_t corr=getAve(0,nSources,1);
  djack_t mVK1,Z2VK1;
  two_pts_fit(Z2VK1,mVK1,corr,TH,tMinFit,tMaxFit,"plots/eff_mass_VKVK_twopts_fit.xmg");
  console<<"Z2: "<<Z2VK1<<endl;
  // djack_t mVK2,Z2VK2;
  // two_pts_fit(Z2VK2,mVK2,corr,TH,22,32);
  
  grace_file_t amu("plots/amu.xmg");
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
      const djack_t cInt=integrate_corr_times_kern_up_to(corr,T,a,upto)*sqr(Za)*1e10;
      const djack_t cSubs1=integrate_corr_times_kern_up_to(corrRefatta1,T,a,THm1)*sqr(Za)*1e10;
      const djack_t cSubs2=integrate_corr_times_kern_up_to(corrRefatta2,T,a,THm1)*sqr(Za)*1e10;
      amuInt[upto]=cInt;
      amuSubs1[upto]=cSubs1;
      amuSubs2[upto]=cSubs2;
    }
  amu.set_xaxis_label("t");
  
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
  corrRefatta1.ave_err().write("plots/corr_refatta1.xmg");
  corrRefatta2.ave_err().write("plots/corr_refatta2.xmg");
  
  MPI_Finalize();
  
  return 0;
}
