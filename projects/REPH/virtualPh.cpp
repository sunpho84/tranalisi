#include <tranalisi.hpp>

const int L=24;
const int ndim=4;

index_t index3pts;

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

struct inv_t
{
  int isolv;
  double mu,th[3];
};

struct einv_t
{
  int isolv;
  double mu,th0[3],tht[3],off;
};

struct combination_t
{
  int i0,it,is;
  double mu1,mu2,off;
  double th0[3],tht[3],ths[3];
};

struct file_head_t
{
  int tmax;
  int x0;
  int stype;
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
};

file_head_t f;

struct data_t
{
  int nc;   ///< Conf id
  int size;
  vector<double> HA,HV;
};

struct data2_t
{
  int nc,size;
};

double Eg(int icomb)
{
  const double th0=f.comb[icomb].th0[2];
  const double tht=f.comb[icomb].tht[2];
  const double mom0=2*M_PI*th0/L;
  const double momt=2*M_PI*tht/L;
  
  cout<<"Th0: "<<th0<<", tht: "<<tht<<endl;
  
  const double kDec=mom0-momt;
  const double kHatDec=2*sin(kDec/2);
  
  return 2*asinh(fabs(kHatDec)/2);
}

void readHeader(raw_file_t& fin)
{
  for(auto& i : {&f.tmax,&f.x0,&f.stype,&f.ninv,&f.neinv,&f.nsolv,&f.nhits,&f.ncomb,&f.ngsm,&f.nqsml,&f.nqsm0,&f.nqsm})
    fin.bin_read(*i);
  
  f.inv.resize(f.ninv);
  f.comb.resize(f.ncomb);
  f.einv.resize(f.neinv);
  
  for(auto& d : {&f.epsgsm,&f.epsqsm})
    fin.bin_read(*d);
   
  fin.bin_read(f.gflv.qhat);
  
  for(auto& d : {&f.gflv.kappa,&f.gflv.mu,&f.gflv.su3csw,&f.gflv.u1csw,&f.gflv.cF,&f.gflv.cF_prime,&f.gflv.th1,&f.gflv.th2,&f.gflv.th3})
    fin.bin_read(*d);
  
  for(int icomb=0;icomb<f.ncomb;++icomb)
    {
      auto& c=f.comb[icomb];
      
      for(auto& i : {&c.i0,&c.it,&c.is})
	fin.bin_read(*i);
      
      for(auto& d : {&c.mu1,&c.mu2,&c.off,&c.th0[0],&c.th0[1],&c.th0[2],&c.tht[0],&c.tht[1],&c.tht[2],&c.ths[0],&c.ths[1],&c.ths[2]})
	fin.bin_read(*d);
    }
  
  for(int inv=0;inv<f.ninv;++inv)
    {
      auto& v=f.inv[inv];
      
      for(auto& d : {&v.mu,&v.th[0],&v.th[1],&v.th[2]})
	fin.bin_read(*d);
    }
}

/// three pts
vector<djvec_t> readData()
{
  raw_file_t fin("data/conf.virtualph.dat","r");
  readHeader(fin);
  
  const int ncorrs=2;
  index3pts.set_ranges({{"corr",ncorrs},{"comb",f.ncomb},{"sl",f.nqsml},{"mu",ndim},{"alpha",ndim},{"ri",2}});
  
  vector<djvec_t> out(index3pts.max(),djvec_t(f.tmax));
  
  data_t data;
  data.size=2*ndim*ndim*f.tmax*f.ncomb*f.nqsml;
  
  int nTotConfs=(fin.size()-fin.get_pos())/(2*data.size*sizeof(double)+sizeof(int));
  cout<<"nTotConfs: "<<nTotConfs<<endl;
  
  const int clustSize=nTotConfs/njacks;
  nTotConfs=clustSize*njacks;
  cout<<"nTotConfs after rounding: "<<nTotConfs<<endl;
  
  vector<double> corr(data.size);
  
  for(int ic=0;ic<nTotConfs;ic++)
    {
      fin.bin_read(data.nc);
      
      const int iclust=ic/clustSize;
      
      //HA,HV
      for(size_t icorr=0;icorr<ncorrs;icorr++)
	{
	  fin.bin_read(corr);
	  
	  for(size_t icomb=0;icomb<(size_t)f.ncomb;icomb++)
	    for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
	      for(size_t mu=0;mu<ndim;mu++)
		for(size_t alpha=0;alpha<ndim;alpha++)
		  for(size_t t=0;t<(size_t)f.tmax;t++)
		    for(size_t ire=0;ire<2;ire++)
		      {
			const int iin=ire+2*(t+f.tmax*(alpha+ndim*(mu+ndim*(isl+f.nqsml*icomb))));
			const int iout=index3pts({icorr,icomb,isl,mu,alpha,ire});
			
			if(icorr==0 and icomb==0 and isl==0 and mu==0 and alpha==0 and ire==0)
			  cout<<"t "<<t<<" "<<corr[iin]<<endl;
			
			out[iout][t][iclust]+=corr[iin];
		      }
	}
    }
  
  for(auto& j : out)
    j.clusterize(clustSize);
  
  return out;
}

int index2pts(int icorr,int iinv2,int iinv1,int isl,int ire)
{
  return ire+2*(isl+f.nqsml*(iinv1+f.ninv*(iinv2+f.ninv*icorr)));
}

vector<djvec_t> readData2()
{
  raw_file_t fin("data/conf.virtualph.dat2","r");
  readHeader(fin);
  
  const int ncorrs=5;
  vector<djvec_t> out(index2pts(ncorrs,0,0,0,0),djvec_t(f.tmax));
  
  data2_t data2;
  data2.size=2*f.tmax*f.ninv*f.ninv*f.nqsml;
  
  int nTotConfs=(fin.size()-fin.get_pos())/(ncorrs*data2.size*sizeof(double)+sizeof(int));
  cout<<"nTotConfs: "<<nTotConfs<<endl;
  
  const int clustSize=nTotConfs/njacks;
  nTotConfs=clustSize*njacks;
  cout<<"nTotConfs after rounding: "<<nTotConfs<<endl;
  
  vector<double> corr(data2.size);
  
  for(int ic=0;ic<nTotConfs;ic++)
    {
      fin.bin_read(data2.nc);
      
      const int iclust=ic/clustSize;
      
      //PP,PA0,PA1,PA2,PA3
      for(int icorr=0;icorr<5;icorr++)
	{
	  fin.bin_read(corr);
	  
	  for(int iinv2=0;iinv2<f.ninv;iinv2++)
	    for(int iinv1=0;iinv1<f.ninv;iinv1++)
	      for(int isl=0;isl<f.nqsml;isl++)
	      for(int t=0;t<f.tmax;t++)
		for(int ire=0;ire<2;ire++)
		  {
		    const int iin=ire+2*(t+f.tmax*(isl+f.nqsml*(iinv1+f.ninv*iinv2)));
		    const int iout=index2pts(icorr,iinv2,iinv1,isl,ire);
		    
		    out[iout][t][iclust]+=corr[iin];
		  }
	}
    }
  
  for(auto& j : out)
    j.clusterize(clustSize);
  
  return out;
}

int main()
{
  set_njacks(10);
  
  const auto threePts=readData();
  const auto twoPts=readData2();
  
  cout<<"f.nqsml: "<<f.nqsml<<endl;
  
  enum{PP,PA0,PA1,PA2,PA3};
  const int is=0,il=1; //index to be fetched from inv list
  const int tmin=6,tmax=8;
  djvec_t mP(f.nqsml);
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    mP[isl]=constant_fit(effective_mass(twoPts[index2pts(PP,is,il,isl,0)].symmetrized()),tmin,tmax,combine("plots/PP_ll_sm%zu.xmg",isl));
  
  constexpr int eps[4][4]={{0,0,0,0},
			   {0,-1,-1,0},
			   {0,+1,-1,0},
			   {0,0,0,0}};
  
  
  enum{HA,HV};
  
  const double eG=Eg(1);
  cout<<"Eg: "<<eG<<endl;
  
  // djvec_t dt(f.tmax/2+1);
  // for(int t=0;t<=f.tmax/2;t++)
  //   dt[t]=exp(-eG*t);
    
  djvec_t dt(f.tmax);
  for(int t=0;t<f.tmax;t++)
    dt[t]=exp(-eG*((t<(f.tmax/2))?t:(f.tmax/2-t)));
    
  for(size_t isl=0;isl<(size_t)f.nqsml;isl++)
    {
      auto load3pts=[&isl,&eps,&threePts](size_t icomb) -> djvec_t
		    {
		      size_t ire=0;
		      
		      /// Charges a given polarization r, and weak current alpha
		      auto load3ptsPol=[&isl,&eps,&threePts,&icomb,&ire](size_t alpha,size_t r) -> djvec_t// r and alpha go from 1 to 2
				       {
					 auto t=[&isl,&threePts,&icomb,&alpha,&ire](size_t mu)
						{
						  size_t i=index3pts({HA,icomb,isl,mu,alpha,ire});
						  
						  return threePts[i];
						};
					 
					 return t(1)*eps[r][1]+t(2)*eps[r][2];
				       };
		      
		      djvec_t res(f.tmax);
		      for(int alpha=1;alpha<=2;alpha++)
			for(int r=1;r<=2;r++)
			  res+=load3ptsPol(alpha,r)/eps[r][alpha];
		      
		      return res;//.symmetrized();
		    };
      
      const djvec_t t0=load3pts(0);
      const djvec_t t1=load3pts(1);
      t1.ave_err().write("plots/t1.xmg");
      
      djvec_t r=t1/t0;
      r*=dt;
      
      // const djack_t xG=2*eG/mP;
      // cout<<"Xg: "<<smart_print(xG)<<endl;
      r.symmetrized().ave_err().write(combine("plots/threePts_sml%zu.xmg",isl));
    }
  
  for(size_t icombo=0;icombo<2;icombo++)
    for(size_t isl=0;isl<1;isl++)
      for(size_t mu=0;mu<ndim;mu++)
	for(size_t alpha=0;alpha<ndim;alpha++)
	for(size_t ri=0;ri<2;ri++)
	  {
	    const size_t i=index3pts({HA,icombo,isl,mu,alpha,ri});
	    auto b=threePts[i];
	    b.ave_err().write(combine("plots/naz/HA_%s.xmg",index3pts.escaped_descr(i).c_str()));
	  }
  dt.ave_err().write("plots/naz/dt.xmg");
  
  return 0;
}
