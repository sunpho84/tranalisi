#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <filesystem>
#include <H5Cpp.h>

#include <mpi.h>
#include <glob.h>

#include "effective.hpp"
#include "fit.hpp"
#include "meas_vec.hpp"

#include <data.hpp>
#include <params.hpp>
#include <rawData.hpp>

using namespace H5;

namespace
{
  double clustSize;
}

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
  
  clustSize=(double)nConfs/njacks;
  console<<"Cluster size: "<<endl;
  
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
  
  err.ave_err().write("plots/err_with_err_"+tag+".xmg");
  corr.ave_err().write("plots/stat_corr_"+tag+".xmg");
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

djvec_t getRawAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes)
{
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
  
  return ave;
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