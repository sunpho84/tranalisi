#include <tranalisi.hpp>

gen_t gen(235423465);

int main()
{
  set_njacks(50);

  const int T=32;
  djvec_t corr(T);
  
  ifstream l("/tmp/l");
  vector<double> d;
  double e;
  while(l>>e)
    d.push_back(e);
  
  const size_t nl=d.size();
  const size_t nc=nl/T/T;
  
  cout<<nl<<" "<<nc<<endl;

  size_t i=0;
  vector<double> corrs(T*nc,0);
  for(size_t ic=0;ic<nc;ic++)
    for(size_t ts=0;ts<T;ts++)
      for(size_t t=0;t<T;t++)
	{
	  size_t tp=(int(t)+T-int(ts))%T;
	  corrs[tp+T*ic]+=d[i];
	  i++;
	  }
  corrs/=T;
  
  for(size_t t=0;t<T;t++)
    jackknivesFill(nc,[&corr,&corrs,t](const size_t& iConf,const size_t& iClust,const double& weight)
    {
      corr[t][iClust]+=weight*corrs[t+T*iConf];
    });
  
  corr.clusterize(nc/50.0).symmetrize();
  constant_fit(effective_mass(corr),8,15,"/tmp/l.xmg");
  
  return 0;
}
