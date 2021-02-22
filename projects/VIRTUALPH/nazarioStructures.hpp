#ifndef _NAZARIO_STRUCTURES_HPP
#define _NAZARIO_STRUCTURES_HPP

#include <index.hpp>
#include <jack.hpp>
#include <meas_vec.hpp>

/// Number of dimension
constexpr int ndim=4;

/// Enumerates the different chirality of the weak current
enum HAV{HA=0,HV=1};

/// Tag of the chirality
constexpr char havTag[2][2]={"A","V"};

/// Flavor of the Nazario header
struct flavour_t
{
  int qhat;
  double kappa;
  double mu;
  double su3csw;
  double u1csw;
  double cF;
  double cF_prime;
  double th1;
  double th2;
  double th3;
};

/// Inverter of the Nazario header
struct inv_t
{
  int isolv;
  double mu,th[3];
};

/// Extended inverstion of the Nazario header
struct einv_t
{
  int isolv;
  double mu,th0[3],tht[3],off;
};

/// Compute combination description
struct combination_t
{
  int i0,it,is;
  double mu1,mu2,off;
  double th0[3],tht[3],ths[3];
};

/// Prints a given combination
ostream& operator<<(ostream& os,const combination_t& c)
{
#define PRINT(A)				\
  os<<" "<<#A<<" "<<c.A<<endl
  
  PRINT(i0);
  PRINT(it);
  PRINT(is);
  PRINT(mu1);
  PRINT(mu2);
  PRINT(off);
  PRINT(th0[0]);
  PRINT(th0[1]);
  PRINT(th0[2]);
  PRINT(tht[0]);
  PRINT(tht[1]);
  PRINT(tht[2]);
  PRINT(ths[0]);
  PRINT(ths[1]);
  PRINT(ths[2]);
  os<<endl;
  
#undef PRINT
  
  return os;
}

/// Data stored in the 3pts file
struct data_t
{
  int nc;   ///< Conf id
  int size;
  vector<double> HA,HV;
};

/// Data stored in the 2pts file
struct data2_t
{
  int nc,size;
};

/// Structure to read Nazario data
struct nazarioReader
{
  const string path;
  
  int tmax;
  int x0;
  int stype;
  int phptype;
  int z0;
  int ninv;
  int neinv;
  int nsolv;
  int nhits;
  int ncomb;
  int ngsm;
  double epsgsm;
  int nqsml,nqsm0,nqsm;
  double epsqsm;
  flavour_t gflv;
  vector<combination_t> comb;
  vector<inv_t> inv;
  vector<einv_t> einv;
  
  const size_t L;
  
  index_t index3pts;
  vector<djvec_t> threePts;
  
  vector<djvec_t> twoPts;
  
  vector<double> eG;
  
  vector<double> kZ;
  
  /// Finds the opposite combination
  size_t iOppComb(const size_t icomb) const
  {
    size_t iout=0;
    
    const combination_t& c=comb[icomb];
    
    bool found=false;
    
    do
      {
	const combination_t& d=comb[iout];
	found=(c.mu1==d.mu2 and
	       c.mu2==d.mu1 and
	       c.th0[2]==d.th0[2] and
	       c.tht[2]==d.tht[2] and
	       c.ths[2]==d.ths[2]);
	
	if(not found) iout++;
      }
    while(iout<comb.size() and not found);
    
    return iout;
  }
  
  /// Reads the header
  void readHeader(raw_file_t& fin)
  {
    for(auto& i : {&tmax,&x0,&stype,&phptype,&z0,&ninv,&neinv,&nsolv,&nhits,&ncomb,&ngsm,&nqsml,&nqsm0,&nqsm})
      fin.bin_read(*i);
    
    cout<<ncomb<<endl;
    
    inv.resize(ninv);
    comb.resize(ncomb);
    einv.resize(neinv);
    
    for(auto& d : {&epsgsm,&epsqsm})
      fin.bin_read(*d);
    
    fin.bin_read(gflv.qhat);
    
    for(auto& d : {&gflv.kappa,&gflv.mu,&gflv.su3csw,&gflv.u1csw,&gflv.cF,&gflv.cF_prime,&gflv.th1,&gflv.th2,&gflv.th3})
      fin.bin_read(*d);
    
    for(int icomb=0;icomb<ncomb;++icomb)
      {
	auto& c=comb[icomb];
	
	for(auto& i : {&c.i0,&c.it,&c.is})
	  fin.bin_read(*i);
	
	for(auto& d : {&c.mu1,&c.mu2,&c.off,&c.th0[0],&c.th0[1],&c.th0[2],&c.tht[0],&c.tht[1],&c.tht[2],&c.ths[0],&c.ths[1],&c.ths[2]})
	  fin.bin_read(*d);
	
	cout<<c<<endl;
      }
    
    // prints opposite of all combo
    for(int icomb=0;icomb<ncomb;++icomb)
      cout<<"Opposite of comb "<<icomb<<": "<<iOppComb(icomb)<<endl;
    
    for(int i=0;i<ninv;i++)
      {
	auto& v=inv[i];
	
	for(auto& d : {&v.mu,&v.th[0],&v.th[1],&v.th[2]})
	  fin.bin_read(*d);
      }
  }
  
  /// Read three pts
  vector<djvec_t> readData()
  {
    raw_file_t fin(path+"/data/conf.virtualph.dat","r");
    readHeader(fin);
    
    const int ncorrs=2;
    
    /// Index to access three point functions in the Nazario ordering
    index3pts.set_ranges({{"corr",ncorrs},{"comb",ncomb},{"sl",nqsml},{"mu",ndim},{"alpha",ndim},{"ri",2}});
    
    vector<djvec_t> out(index3pts.max(),djvec_t(tmax));
    
    data_t data;
    data.size=2*ndim*ndim*tmax*ncomb*nqsml;
    
    int nTotConfs=(fin.size()-fin.get_pos())/(2*data.size*sizeof(double)+sizeof(int));
    cout<<"nTotConfs: "<<nTotConfs<<endl;
    
    const int clustSize=nTotConfs/njacks;
    nTotConfs=clustSize*njacks;
    cout<<"nTotConfs after rounding: "<<nTotConfs<<endl;
    
    vector<double> corr(data.size);
    
    for(int ic=0;ic<nTotConfs;ic++)
      {
	fin.bin_read(data.nc);
	//cout<<"nc: "<<data.nc<<endl;
	
	const int iclust=ic/clustSize;
	
	//HA,HV
	for(size_t icorr=0;icorr<ncorrs;icorr++)
	  {
	    fin.bin_read(corr);
	    
	    for(size_t icomb=0;icomb<(size_t)ncomb;icomb++)
	      for(size_t isl=0;isl<(size_t)nqsml;isl++)
		for(size_t mu=0;mu<ndim;mu++)
		  for(size_t alpha=0;alpha<ndim;alpha++)
		    for(size_t t=0;t<(size_t)tmax;t++)
		      for(size_t ire=0;ire<2;ire++)
			{
			  const int iin=ire+2*(t+tmax*(alpha+ndim*(mu+ndim*(isl+nqsml*icomb))));
			  const int iout=index3pts({icorr,icomb,isl,mu,alpha,ire});
			  
			  // if(icorr==0 and icomb==0 and isl==0 and mu==0 and alpha==0 and ire==0)
			  //   cout<<"t "<<t<<" "<<corr[iin]<<endl;
			  
			  out[iout][t][iclust]+=corr[iin];
			}
	  }
      }
    
    for(auto& j : out)
      j.clusterize(clustSize);
    
    return out;
  }
  
  /// Read two points
  vector<djvec_t> readData2()
  {
    raw_file_t fin(path+"/data/conf.virtualph.dat2","r");
    readHeader(fin);
    
    const int ncorrs=5;
    index_t index2pts({{"Corr",ncorrs},{"Inv2",ninv},{"Inv1",ninv},{"Sm",2},{"T",tmax},{"RI",2}});
    
    vector<djvec_t> out(2,djvec_t(tmax));
    
    data2_t data2;
    data2.size=2*tmax*ninv*ninv*nqsml;
    
    int nTotConfs=(fin.size()-fin.get_pos())/(ncorrs*data2.size*sizeof(double)+sizeof(int));
    cout<<"nTotConfs: "<<nTotConfs<<endl;
    
    const int clustSize=nTotConfs/njacks;
    nTotConfs=clustSize*njacks;
    cout<<"nTotConfs after rounding: "<<nTotConfs<<endl;
    
    vector<double> corr(data2.size*ncorrs);
    
    for(int ic=0;ic<nTotConfs;ic++)
      {
	fin.bin_read(data2.nc);
	
	const int iclust=ic/clustSize;
	
	//PP,PA0,PA1,PA2,PA3
	fin.bin_read(corr);
	
	constexpr size_t is=0,il=1; //hack
	
	for(size_t isl=0;isl<(size_t)nqsml;isl++)
	  for(size_t t=0;t<(size_t)tmax;t++)
	    out[isl][t][iclust]+=corr[index2pts({0,is,il,isl,t,0})];
      }
    
    for(auto& j : out)
      j.clusterize(clustSize);
    
    return out;
  }
  
  /// Gets the 3pts for a given insertion, combination, smearing level, mu and alpha
  djvec_t get3pts(const HAV hAV,const size_t icomb,const size_t isl,const size_t mu,const size_t alpha)
  {
    const size_t ire=(hAV==HAV::HA)?0:1;
    
    const size_t ic=index3pts({hAV,icomb,isl,mu,alpha,ire});
    //cout<<"ic: "<<ic<<" , "<<index3pts.descr(ic)<<endl;
    
    djvec_t out=threePts[ic];
    
    for(int it=0;it<tmax;it++)
      {
	const double dt=fabs(tmax/2.0-it);
	const double arg=-eG[icomb]*dt;
	out[it]/=exp(arg);
      }
    
    return out.symmetrized((hAV==HV)?-1:+1);
  }
  
  /// Momentum of the photon for a given combination of theta
  double Kz(const double th0,const double tht)
  {
    const double mom0=2*M_PI*th0/L;
    const double momt=2*M_PI*tht/L;
    
    //cout<<"Th0: "<<th0<<", tht: "<<tht<<endl;
    
    const double kDec=mom0-momt;
    
    return kDec;
  }
  
  /// momenum of the photon for a given combination
  double Kz(const int icomb)
  {
    return Kz(comb[icomb].th0[2],comb[icomb].tht[2]);
  }
  
  /// Energy of the photon for a given combination
  double Eg(const int icomb)
  {
    const double kDec=Kz(icomb);
    const double kHatDec=2*sin(kDec/2);
    
    return 2*asinh(fabs(kHatDec)/2);
  }
  
  nazarioReader(const string path,const size_t L) : path(path),L(L)
  {
    threePts=readData();
    twoPts=readData2();
    
    eG.resize(ncomb);
    kZ.resize(ncomb);
    for(int icomb=0;icomb<ncomb;++icomb)
      {
	kZ[icomb]=Kz(icomb);
	eG[icomb]=Eg(icomb);
	cout<<"Eg["<<icomb<<"]: "<<eG[icomb]<<endl;
      }
  }
};

#endif
