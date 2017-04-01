#include <tranalisi.hpp>

size_t T=80,TH=T/2;

int main(int narg,char **arg)
{
  set_njacks(5);
  djvec_t corr(TH+1);
  
  for(size_t ijack=0;ijack<5;ijack++)
    {
      size_t iconf=ijack+110;
      raw_file_t file(combine("onlinemeas.%06zu",iconf),"r");
      for(size_t t=0;t<=TH;t++)
	{
	  double dum;
	  file.expect({"1","1",to_string(t)});
	  file.read(corr[t][ijack]);
	  file.read(dum);
	  cout<<t<<" "<<ijack<<" "<<corr[t][ijack]<<endl;
	}
      file.close();
      cout<<"ff"<<endl;
    }
  
  corr.clusterize();
  effective_mass(corr).ave_err().write("eff.xmg");
  
  return 0;
}
