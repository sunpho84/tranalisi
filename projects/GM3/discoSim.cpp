#include <tranalisi.hpp>

using namespace std;

const size_t nSources=2;
const size_t nConfs=401;
const size_t T=160;

int main()
{
  set_njacks(nConfs);
  
  index_t id({{"Source",nSources},{"Conf",nConfs},{"T",T}});
  
  vector<double> rawData(id.max());
  
  for(size_t iSource=0;iSource<nSources;iSource++)
    {
      //const string finPath="C80_gen_light_S"+std::to_string(iSource+1)+"_g5.txt";
      const string finPath="C80_std_light_S"+std::to_string(iSource+1)+"_g5g0.txt";
      ifstream fin(finPath);
      
      if(not fin.good())
	CRASH("Error opening %s",finPath.c_str());
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	{
	  string dum;
	  fin>>dum>>dum;
	  
	  for(size_t t=0;t<T;t++)
	    {
	      double dum;
	      //if(not (fin>>rawData[id({iSource,iConf,t})]>>dum))
	      if(not (fin>>dum>>rawData[id({iSource,iConf,t})]))
		CRASH("Reading conf %zu t %zu",iConf,t);
	    }
	}
    }
  
  vector<double> all(nConfs*T,0.0),unb(nConfs*T,0.0);
  
  for(size_t iConf=0;iConf<nConfs;iConf++)
    {
      for(size_t t1=0;t1<T;t1++)
	for(size_t dt=0;dt<T;dt++)
	  {
	    const size_t t2=(t1+dt)%T;
	    
	    double s12=0,s1=0,s2=0;
	    for(size_t iSource=0;iSource<nSources;iSource++)
	      {
		const double& x1=rawData[id({iSource,iConf,t1})];
		const double& x2=rawData[id({iSource,iConf,t2})];
		s12+=x1*x2;
		s1+=x1;
		s2+=x2;
	      }
	    
	    all[dt+T*iConf]+=s1*s2;
	    unb[dt+T*iConf]+=s1*s2-s12;
	  }
    }
  
  all/=T*sqr(nSources);
  unb/=T*nSources*(nSources-1);
  
  djvec_t a(T,0.0),u(T,0.0);
  
  jackknivesFill(nConfs,
		 [&a,
		  &u,
		  &all,
		  &unb](const size_t& iConf,
			const size_t& iClust,
			const double& w)
		 {
		   for(size_t t=0;t<T;t++)
		     {
		       a[t][iClust]+=all[t+T*iConf]*w;
		       u[t][iClust]+=unb[t+T*iConf]*w;
		     }
		 });
  
  a.clusterize((double)nConfs/njacks);
  u.clusterize((double)nConfs/njacks);
  
  a.ave_err().write("plots/all.xmg");
  u.ave_err().write("plots/unb.xmg");
  (a-u).ave_err().write("plots/bias.xmg");
  
  effective_mass(a.symmetrized()).ave_err().write("plots/alleff.xmg");
  effective_mass(u.symmetrized()).ave_err().write("plots/unbeff.xmg");
  
  return 0;
}
