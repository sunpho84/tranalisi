#include <tranalisi.hpp>

const size_t T=48;
const int nconfs=150;

djvec_t load(string path)
{
  djvec_t out(T);
  out=0.0;
  
  ifstream in(path);
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    for(size_t i=0;i<T;i++)
      for(size_t j=0;j<T;j++)
	{
	  double re,im;
	  if(!(in>>re>>im)) CRASH("Unable to read %s conf %zu, %zu %zu",path.c_str(),iconf,i,j);
	  size_t t=(T+i-j)%T;
	  out[t][iconf]+=re;
	}
  
  out.clusterize(1);
  
  out*=-1.0/(T*T);
  out[0]=out[1];
  
  return out.symmetrized(1);
}

djvec_t load_LO(string path)
{
  djvec_t out(T);
  out=0.0;
  
  ifstream in(path);
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    for(size_t t=0;t<T;t++)
      {
	double re,im;
	if(!(in>>re>>im)) CRASH("Unable to read %s conf=%zu, t=%zu",path.c_str(),iconf,t);
	out[t][iconf]=re;
      }
  
  out.clusterize();
  out[0]=out[1];
  
  return out.symmetrized(1);
}

int main()
{
  set_njacks(nconfs);
  cout<<"NConfs: "<<nconfs<<endl;
  
  djvec_t manette=load("data/manette");
  manette.ave_err().write("plots/manette.xmg");
  
  djvec_t mesone00=load_LO("data/mesone00");
  effective_mass(mesone00).ave_err().write("plots/mesone00.xmg");
  
  djvec_t mesone01=load_LO("data/mesone01");
  effective_mass(mesone01).ave_err().write("plots/mesone01.xmg");
  
  djvec_t ratio=manette/mesone01;
  ratio.ave_err().write("plots/ratio.xmg");
  forward_derivative(ratio).ave_err().write("plots/ratio_der.xmg");
  
  // plot_ratio<<forward_derivative(ratio).ave_err()<<endl;
  // plot_ratio.new_data_set();
  
  // djack_t M,A,SL;
  // two_pts_with_ins_ratio_fit(M,A,SL,mesone01,manette,T/2,8,15,"plots/pi_mass.xmg","plots/pi_slope.xmg");
  
  
  return 0;
}
