#include <tranalisi.hpp>

const size_t T=48,L=24,VSPAT=L*L*L,V=VSPAT*T;
const int nconfs=138;

djvec_t load(string path)
{
  djvec_t out(T);
  out=0.0;
  
  array<int,T> n;
  for(auto &nt : n) nt=0;
  
  ifstream in(path);
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    {
      array<array<complex<double>,T>,T> d;
      for(size_t i=0;i<T;i++)
	for(size_t j=0;j<T;j++)
	  {
	    double re,im;
	    if(!(in>>re>>im)) CRASH("Unable to read %s conf %zu, %zu %zu",path.c_str(),iconf,i,j);
	    d[i][j]=complex<double>(re*V/12,im*V/12);
	  }
      
      for(size_t i=0;i<T;i++)
	for(size_t j=0;j<T;j++)
	  {
	    size_t t=(T+i-j)%T;
	    
	    for(size_t iso1=0;iso1<T;iso1++)
	      for(size_t iso2=0;iso2<T;iso2++)
		if(iso1!=iso2)
		  {
		    out[t][iconf]+=-(d[iso1][i]*conj(d[iso2][j])).real();
		    if(iconf==0) n[t]++;
		  }
	  }
    }
  out.clusterize();
  
  for(size_t t=0;t<T;t++) out[t]/=n[t];
  
   double dum;
   if(in>>dum) CRASH("Should have reached EOF, obtained %lg",dum);
  
  return -out.symmetrized(1);
}

int main()
{
  set_njacks(nconfs);
  djvec_t data=load("data/manette2");
  
  data.ave_err().write("plots/manette2.xmg");
  
  return 0;
}
