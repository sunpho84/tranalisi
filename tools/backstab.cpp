#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=32,TH=T/2;

djvec_t read(const string &path)
{
  ifstream fin(path);
  if(not fin.good()) CRASH("Opening %s",path.c_str());
  
  djvec_t vec(T);
  vec=0.0;
  
  size_t n;
  if(not (fin>>n)) CRASH("Reading n");
  size_t clust_size=n/njacks;
  
  for(size_t i=0;i<n;i++)
    {
      size_t ijack=i/clust_size;
      for(size_t t=0;t<=TH;t++)
	{
	  double in;
	  if(not (fin>>in)) CRASH("Error reading %zu %zu",t,ijack);
	  vec[t][ijack]+=in;
	  vec[(T-t)%T][ijack]+=in;
	}
    }
  vec.clusterize(clust_size);
  
  return vec.symmetrized();
}

int main(int narg,char **arg)
{
  set_njacks(15);
  
  djvec_t vecA=read("1/sanfomeas");
  djvec_t vecB=read("1m/sanfomeas");
  
  grace_file_t gr("1.xmg");
  gr.write_vec_ave_err(effective_mass(vecA).ave_err());
  gr.write_vec_ave_err(effective_mass(vecB).ave_err());
  
  djvec_t test=(vecA/vecB)-1.0;
  test.ave_err().write("1diff.xmg");
  
  return 0;
}
