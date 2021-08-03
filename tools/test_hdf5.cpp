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
constexpr size_t nGammaComb=4;
constexpr char gammaCombTag[nGammaComb][5]={"P5P5","VKVK","TKTK","_KVK"};
constexpr size_t nMes=3;
constexpr char mesTag[nMes][3]={"uu","ud","dd"};

index_t idData;
index_t idData_loader;
index_t idOpenData_loader;
vector<double> rawData;
vector<int> confMap;

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
  
  idData.set_ranges({{"Confs",nConfs},{"Source",nSources},{"GammComb",nGammaComb},{"Mes",nMes},{"T",THp1}});
  idData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Gamma",16}});
  idOpenData_loader.set_ranges({{"Mes",nMes},{"T",T},{"Id1",4},{"Id2",4},{"Id3",4},{"Id4",4},{"Ri",2}});
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
    
    dataIn.resize(idData.max(),0.0);
    
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

void loadRawData(int narg,char** arg)
{
  MPI_Init(&narg,&arg);
  MPI_Comm_size(MPI_COMM_WORLD,&nMPIranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);
  
  const vector<string> possibleConfsList=getConfsList(confsPattern);
  const vector<string> sourcesList=getSourcesList(possibleConfsList.front());
  nSources=sourcesList.size();
  
  vector<string> confsList;
  for(const string& conf : possibleConfsList)
    {
      bool exists=true;
      
      if(MPIrank==0)
	{
	  string confPath;
	  size_t iSource=0;
	  while(iSource<nSources and exists)
	    {
	      confPath=conf+"/"+sourcesList[iSource];
	      if(not std::filesystem::exists(confPath))
		exists=false;
	      iSource++;
	    }
	  
	  cout<<"Conf "<<conf<<(exists?" accepted":(" discarded, "+confPath+" not found"))<<endl;
	}
      MPI_Bcast(&exists,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
      
      if(exists)
	  confsList.push_back(conf);
    }
  
  if(MPIrank==0 and confsList.size()!=possibleConfsList.size())
    cout<<"Confs list resized from "<<possibleConfsList.size()<<"to: "<<confsList.size()<<endl;
  
  nConfs=confsList.size();
  
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
		  double& out=rawData[idData({iConf,iSource,igamma_out,iMes,tOut})];
		  
		  out+=in;
		}

  const Gamma_data_t Gamma_data[nGamma]=
    {{{0,1,2,3},{1,1,1,1},{0,0,0,0}},
     {{3,2,1,0},{0,0,0,0},{-1,-1,1,1}},
     {{3,2,1,0},{-1,1,1,-1},{0,0,0,0}},
     {{2,3,0,1},{0,0,0,0},{-1,1,1,-1}},
     {{2,3,0,1},{-1,-1,-1,-1},{0,0,0,0}},
     {{0,1,2,3},{1,1,-1,-1},{0,0,0,0}},
     {{3,2,1,0},{0,0,0,0},{1,1,1,1}},
     {{3,2,1,0},{1,-1,1,-1},{0,0,0,0}},
     {{2,3,0,1},{0,0,0,0},{1,-1,1,-1}},
     {{2,3,0,1},{1,1,-1,-1},{0,0,0,0}},
     {{1,0,3,2},{0,0,0,0},{-1,-1,1,1}},
     {{1,0,3,2},{-1,1,1,-1},{0,0,0,0}},
     {{0,1,2,3},{0,0,0,0},{-1,1,1,-1}},
     {{1,0,3,2},{0,0,0,0},{1,1,1,1}},
     {{1,0,3,2},{1,-1,1,-1},{0,0,0,0}},
     {{0,1,2,3},{0,0,0,0},{1,-1,1,-1}}};
  
    // A(i)=(GSO)_{ij(i)} (G5)_{j(i)}
    // B(k)=(G5)_k (GSI)_{kl(k)}
   const array<array<size_t,3>,3> map{array<size_t,3>{3,1,1},{3,2,2},{3,3,3}};
  //const array<array<size_t,3>,1> map{array<size_t,3>{3,5,5}};
	  for(size_t iMes=0;iMes<nMes;iMes++)
	    for(size_t tIn=0;tIn<T;tIn++)
	      {
		const size_t tOut=
		  ((tIn>=TH)?(T-tIn):tIn);

		for(const auto& m : map)
		  {
		    const size_t iGammaOut=m[0];
		  
		  const size_t iGammaIn1=m[1];
		  const size_t iGammaIn2=m[2];
		  
		  const auto& g1=Gamma_data[iGammaIn1];
		  const auto& g2=Gamma_data[iGammaIn2];
		  const auto& g5=Gamma_data[5];
		  
		  double& out=rawData[idData({iConf,iSource,iGammaOut,iMes,tOut})];
		  
		  for(size_t nu=0;nu<4;nu++)
		    for(size_t rh=0;rh<4;rh++)
		      {
			const size_t si=g1[0][nu];
			const size_t mu=g2[0][rh];
			complex<double> c1,c2,c15,c25;
			for(size_t ri=0;ri<2;ri++)
			  {
			    ((double*)&c1)[ri]=g1[ri+1][nu];
			    ((double*)&c2)[ri]=g2[ri+1][rh];
			    ((double*)&c15)[ri]=g5[ri+1][nu];
			    ((double*)&c25)[ri]=g5[ri+1][mu];
			  }
			const complex<double> c=c1*c2*c15*c25;
			
			// if(tIn==0 and iSource==0 and iMes==0 and iConf==0)
			//   printf("%zu %zu %zu %zu ",mu,nu,si,rh);
			for(size_t ri=0;ri<2;ri++)
			  {
			    const double& in=loader.openDataIn[idOpenData_loader({iMes,tIn,mu,nu,si,rh,ri})];
			// if(tIn==0 and iSource==0 and iMes==0 and iConf==0)
			//   printf(" %lg",in);
			    
			    out+=in*((double*)&c)[ri];
			  }
			// if(tIn==0 and iSource==0 and iMes==0 and iConf==0)
			//   printf("\n");
		      }
		  }
	      }
      cout<<rawData[idData({iConf,iSource,3,0,0})]<<" "<<rawData[idData({iConf,iSource,0,0,0})]<<endl;
	}
    }
  
  MPI_Allreduce(MPI_IN_PLACE,&rawData[0],idData.max(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  for(size_t i=0;i<idData.max();i++)
    {
      const std::vector<size_t> coords=idData(i);
      
      const size_t &t=coords.back();
      const size_t &ig=coords[2];
      
      const double norm=((t==0 or t==TH)?1:2)*
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

int main(int narg,char **arg)
{
  raw_file_t input("input.txt","r");
  T=input.read<size_t>("T");
  L=input.read<size_t>("L");
  confsPattern=input.read<string>("ConfsPattern");
  output=input.read<string>("Output");
  
  if(not file_exists(output))
    loadRawData(narg,arg);
  else
    loadData();
  
  if(file_exists("map.txt"))
    {
      raw_file_t confMapFile("map.txt","r");
      for(size_t iConf=0;iConf<nConfs;iConf++)
	confMap.push_back(confMapFile.read<int>());
    }
  else
    confMap=vector_up_to<int>(nConfs);
  
  set_njacks(nConfs);
  
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
      
      const djack_t m=constant_fit(effective_mass(aveCorr),(iGammaComb==0)?25:20,32,"plots/eff_mass_"+cTag+".xmg");
      
      switch(iGammaComb)
	{
	case 0:
	  mP5=m;
	  cout<<"mP5: "<<m.ave_err()<<endl;
	  break;
	case 1:
	  {
	    const djack_t E2p=2*sqrt(sqr(mP5)+sqr(2*M_PI/L));
	    cout<<"EVK: "<<E2p.ave_err()<<endl;
	    cout<<"mVK: "<<m.ave_err()<<endl;
	  }
	  break;
	case 2:
	  cout<<"mTK: "<<m.ave_err()<<endl;
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
	  
	  const size_t mapCorr[4]={0,1,1,1};
	  const char corrName[4][5]={"P5P5","V1V1","V2V2","V3V3"};
	  for(size_t _iCorr=0;_iCorr<4;_iCorr++)
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
  
  const djvec_t corr=(getAve(0,nSources,1,0)+getAve(0,nSources,1,2))/2;
  djack_t mVK1,Z2VK1;
  two_pts_fit(Z2VK1,mVK1,corr,TH,18,32);
  djack_t mVK2,Z2VK2;
  two_pts_fit(Z2VK2,mVK2,corr,TH,22,32);
  
  const double a=0.414,Za=0.746;
  
  grace_file_t amu("plots/amu.xmg");
  djvec_t amuInt(TH),amuSubs1(TH),amuSubs2(TH);
  for(size_t upto=0;upto<TH;upto++)
    {
      djvec_t corrRefatta1=corr;
      djvec_t corrRefatta2=corr;
      for(size_t t=upto;t<TH;t++)
	{
	  corrRefatta1[t]=two_pts_corr_fun(Z2VK1,mVK1,TH,t,+1);
	  corrRefatta2[t]=two_pts_corr_fun(Z2VK2,mVK2,TH,t,+1);
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
  amu.set_legend("Substituting 2pts beyond t");
  
  amu.write_vec_ave_err(amuSubs2.ave_err());
  amu.set_no_line();
  amu.set_all_colors(grace::ORANGE);
  amu.set_legend("Substituting 2pts beyond t");
  
  amu.new_data_set();
  amu.set_legend("BMW light connected");
  amu.set_all_colors(grace::GREEN4);
  amu.write_constant_band(0,TH,djack_t(gauss_filler_t{633.7,5.0,23423}));
  
  return 0;
}
