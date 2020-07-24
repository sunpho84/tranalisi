#include <tranalisi.hpp>

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

struct data_t
{
  int nc;   ///< Conf id
  int size;
  vector<double> HA,HV;
};

struct data2_t
{
  int nc,size;
  vector<double> PP,PA0,PA1,PA2,PA3;
};

file_head_t readHeader(raw_file_t& fin)
{
  file_head_t f;
  
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
  
  return f;
}

void readData()
{
  raw_file_t fin("data/conf.virtualph.dat","r");
  file_head_t f=readHeader(fin);
  
  data_t data;
  data.size=32*f.tmax*f.ncomb*f.nqsml;
  
  const int nTotConfs=(fin.size()-fin.get_pos())/(2*data.size*sizeof(double)+sizeof(int));
  cout<<"nTotConfs: "<<nTotConfs<<endl;
  
  for(int ic=0;ic<nTotConfs;ic++)
    {
      fin.bin_read(data.nc);
      for(auto& h : {&data.HA,&data.HV})
	{
	  h->resize(data.size);
	  fin.bin_read(*h);
	}
    }
}

void readData2()
{
  raw_file_t fin("data/conf.virtualph.dat2","r");
  file_head_t f=readHeader(fin);
  
  data2_t data2;
  data2.size=2*f.tmax*f.ninv*f.ninv*f.nqsml;
  const int nTotConfs=(fin.size()-fin.get_pos())/(5*data2.size*sizeof(double)+sizeof(int));
  cout<<"nTotConfs: "<<nTotConfs<<endl;
  
  for(int ic=0;ic<nTotConfs;ic++)
    {
      fin.bin_read(data2.nc);
      
      for(auto& h : {&data2.PP,&data2.PA0,&data2.PA1,&data2.PA2,&data2.PA3})
	{
	  h->resize(data2.size);
	  fin.bin_read(*h);
	}
    }
}

int main()
{
  readData();
  readData2();
  
  return 0;
}
