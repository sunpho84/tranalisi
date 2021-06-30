#include <tranalisi.hpp>

const size_t T=128;
const size_t nconfs=28;
const index_t id({{"r2",2},{"r1",2},{"corr",61},{"T",T},{"ri",2}});

std::pair<vector<double>,vector<double>> read(const char* path)
{
  ifstream input(path);
  
  if(not input.good())
    CRASH("unable to open %s",path);
  
  vector<double> data;
  
  string is;
  while(input>>is)
    {
      istringstream iss(is);
      double d;
      if(iss>>d)
	data.push_back(d);
    }
  
  if(data.size()!=id.max())
    CRASH("mismatch, %zu %zu",data.size(),id.max());
  
  vector<double> VV(T/2+1,0.0),PP(T/2+1,0.0);
  for(size_t r1=0;r1<2;r1++)
    //for(size_t r2=0;r2<2;r2++)
      for(size_t _t=0;_t<T;_t++)
	{
	  const size_t t=std::min(T-_t,_t);
	  const int c=(t==0)?2:1;
	  
	  const size_t r2=r1;
	  for(size_t corr : {35,36,37})
	    VV[t]+=data[id({r2,r1,corr,_t,0})]*c;
	  PP[t]+=data[id({r2,r1,5,_t,0})]*c;
	}
  
  for(auto& pp : PP)
    pp/=2*2;
  
  for(auto& vv : VV)
    vv/=2*2*3;
  
  return {PP,VV};
}

void an(const vector<double>& data,const char* tag)
{
  const size_t clust_size=nconfs/njacks;
  
  djvec_t c((T/2+1)*4),a((T/2+1)*2);
  index_t id({{"copy",2},{"conf",nconfs},{"T",T/2+1}});
  for(size_t t=0;t<=T/2;t++)
    for(size_t ijack=0;ijack<njacks+1;ijack++)
      for(size_t iconf=0;iconf<nconfs;iconf++)
	if(iconf<clust_size*ijack or iconf>=(clust_size+1)*ijack)
	  for(size_t icopy=0;icopy<2;icopy++)
	    {
	      const double& x=data[id({icopy,iconf,t})];
	      a[icopy+2*t][ijack]+=x;
	      
	      for(size_t jcopy=0;jcopy<2;jcopy++)
		{
		  const double& y=data[id({jcopy,iconf,t})];
		  c[jcopy+2*(icopy+2*t)][ijack]+=x*y;
		}
	    }
  
  djvec_t ave(T/2+1),err(T/2+1),corr(T/2+1);
  for(size_t t=0;t<=T/2;t++)
    {
      for(size_t icopy=0;icopy<2;icopy++)
	{
	  for(size_t ijack=0;ijack<njacks;ijack++)
	    a[icopy+2*t][ijack]/=nconfs-clust_size;
	  a[icopy+2*t][njacks]/=nconfs;
	}
      
      for(size_t icopy=0;icopy<2;icopy++)
	for(size_t jcopy=0;jcopy<2;jcopy++)
	  {
	    for(size_t ijack=0;ijack<njacks;ijack++)
	      c[jcopy+2*(icopy+2*t)][ijack]/=nconfs-clust_size;
	    c[jcopy+2*(icopy+2*t)][njacks]/=nconfs;
	  }
      
      for(size_t icopy=0;icopy<2;icopy++)
	for(size_t jcopy=0;jcopy<2;jcopy++)
	  c[jcopy+2*(icopy+2*t)]-=a[icopy+2*t]*a[jcopy+2*t];
	
      for(size_t icopy=0;icopy<2;icopy++)
	{
	  ave[t]+=a[icopy+2*t]/2;
	  err[t]+=sqr(c[icopy+2*(icopy+2*t)]);
	}
      err[t]=sqrt(err[t]);
      corr[t]=c[0+2*(1+2*t)]/sqrt(c[0+2*(0+2*t)]*c[1+2*(1+2*t)]);
    }
  
  ave.ave_err().write(combine("/tmp/%s_ave.xmg",tag));
  err.ave_err().write(combine("/tmp/%s_err.xmg",tag));
  corr.ave_err().write(combine("/tmp/%s_corr.xmg",tag));
}

int main(int narg,char** arg)
{
  // if(narg<2)
  //   CRASH("use %s file",arg[0]);
  
  vector<double> PP,VV;
  
  for(const char* p: {"out","run2/out"})
    for(size_t iconf=0;iconf<nconfs;iconf++)
    {
      vector<double> _PP,_VV;
      tie(_PP,_VV)=read(combine("%s/%04zu/mes_contr_2pts_ll",p,iconf+1).c_str());
      
      for(const double& pp : _PP)
	PP.push_back(pp);
      for(const double& vv : _VV)
	VV.push_back(vv);
    }
  
  set_njacks(28);
  
  an(PP,"PP");
  an(VV,"VV");
  
  return 0;
}
