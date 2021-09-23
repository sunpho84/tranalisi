#include <tranalisi.hpp>

const size_t T=128;
const size_t nconfs=196;
const size_t ncopies=3;
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
	  const int c=(t==0 or t==T/2)?2:1;
	  
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
  
  djvec_t c((T/2+1)*ncopies*ncopies),a((T/2+1)*ncopies),aa((T/2+1)),cc(T/2+1);
  index_t id({{"copy",ncopies},{"conf",nconfs},{"T",T/2+1}});
  for(size_t t=0;t<=T/2;t++)
    for(size_t ijack=0;ijack<njacks+1;ijack++)
      for(size_t iconf=0;iconf<nconfs;iconf++)
	if(iconf<clust_size*ijack or iconf>=clust_size*(ijack+1))
	  {
	    double copy_ave=0;
	    for(size_t icopy=0;icopy<ncopies;icopy++)
	      copy_ave+=data[id({icopy,iconf,t})];
	    copy_ave/=ncopies;
	    aa[t][ijack]+=copy_ave;
	    cc[t][ijack]+=copy_ave*copy_ave;
	    
	    for(size_t icopy=0;icopy<ncopies;icopy++)
	      {
		const double& x=data[id({icopy,iconf,t})];
		a[icopy+ncopies*t][ijack]+=x;
		
		for(size_t jcopy=0;jcopy<ncopies;jcopy++)
		  {
		    const double& y=data[id({jcopy,iconf,t})];
		    c[jcopy+ncopies*(icopy+ncopies*t)][ijack]+=x*y;
		  }
	      }
	  }
  
  djvec_t ave(T/2+1),err(T/2+1),corr(T/2+1);
  for(size_t t=0;t<=T/2;t++)
    {
      for(size_t ijack=0;ijack<njacks;ijack++)
	{
	  aa[t][ijack]/=nconfs-clust_size;
	  cc[t][ijack]/=nconfs-clust_size;
	}
      aa[t][njacks]/=nconfs;
      cc[t][njacks]/=nconfs;
      cc[t]-=aa[t]*aa[t];
      cc[t]=sqrt(cc[t]/(nconfs-1));
      
      for(size_t icopy=0;icopy<ncopies;icopy++)
	{
	  for(size_t ijack=0;ijack<njacks;ijack++)
	    a[icopy+ncopies*t][ijack]/=nconfs-clust_size;
	  a[icopy+ncopies*t][njacks]/=nconfs;
	}
      
      for(size_t icopy=0;icopy<ncopies;icopy++)
	for(size_t jcopy=0;jcopy<ncopies;jcopy++)
	  {
	    for(size_t ijack=0;ijack<njacks;ijack++)
	      c[jcopy+ncopies*(icopy+ncopies*t)][ijack]/=nconfs-clust_size;
	    c[jcopy+ncopies*(icopy+ncopies*t)][njacks]/=nconfs;
	  }
      
      for(size_t icopy=0;icopy<ncopies;icopy++)
	for(size_t jcopy=0;jcopy<ncopies;jcopy++)
	  c[jcopy+ncopies*(icopy+ncopies*t)]-=a[icopy+ncopies*t]*a[jcopy+ncopies*t];
	
      ave[t]+=a[0+ncopies*t];
      err[t]+=c[0+ncopies*(0+ncopies*t)];
      
      err[t]=sqrt(err[t]/(nconfs-1));
      corr[t]=0;
      for(size_t icopy=0;icopy<ncopies;icopy++)
	for(size_t jcopy=icopy+1;jcopy<ncopies;jcopy++)
	  corr[t]+=c[icopy+ncopies*(jcopy+ncopies*t)]/sqrt(c[icopy+ncopies*(icopy+ncopies*t)]*c[jcopy+ncopies*(jcopy+ncopies*t)]);
      corr[t]/=ncopies*(ncopies-1)/2.0;
    }
  
  // ave.ave_err().write(combine("/tmp/%s_ave.xmg",tag));
  // aa.ave_err().write(combine("/tmp/%s_combo_ave.xmg",tag));
  // cc.ave_err().write(combine("/tmp/%s_combo_err.xmg",tag));
  // err.ave_err().write(combine("/tmp/%s_err.xmg",tag));
  corr.ave_err().write(combine("plots/%s_stat_corr.xmg",tag));
}

int main(int narg,char** arg)
{
  // if(narg<2)
  //   CRASH("use %s file",arg[0]);
  
  vector<double> PP,VV;
  
  for(const char* p: {"out","run2/out","run3/out"})
    for(size_t iconf=0;iconf<nconfs;iconf++)
    {
      cout<<p<<" "<<iconf<<endl;
      
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
