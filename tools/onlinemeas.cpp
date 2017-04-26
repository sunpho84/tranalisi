#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  raw_file_t input("input.txt","r");
  
  size_t T=input.read<size_t>("T");
  size_t TH=T/2;
  size_t nconfs=input.read<size_t>("NConfs");
  size_t ntherm=input.read<size_t>("NTherm");
  set_njacks(input.read<size_t>("NJacks"));
  size_t nposs_meas=nconfs-ntherm;
  cout<<"Nposs_meas: "<<nposs_meas<<endl;
  size_t clust_size=nposs_meas/njacks;
  cout<<"Clust_size: "<<clust_size<<endl;
  size_t nmeas=clust_size*njacks;
  cout<<"Nmeas: "<<nmeas<<endl;
  
  djvec_t corr(TH+1);
  
  for(size_t imeas=0;imeas<nmeas;imeas++)
    {
      size_t conf_id=ntherm+imeas;
      size_t ijack=imeas/clust_size;
      raw_file_t file(combine("onlinemeas.%06zu",conf_id),"r");
      for(size_t t=0;t<=TH;t++)
	{
	  double dum;
	  file.expect({"1","1",to_string(t)});
	  double temp;
	  file.read(temp);
	  corr[t][ijack]+=temp;
	  file.read(dum);
	}
      file.close();
    }
  
  corr.clusterize(clust_size);
  cout<<constant_fit(effective_mass(corr),12,24,"eff.xmg").ave_err()<<endl;
  
  return 0;
}
