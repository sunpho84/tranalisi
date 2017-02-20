#include <tranalisi.hpp>

const size_t L=24;
const size_t T=48;
const int nconfs=150;

djvec_t load(string path)
{
  djvec_t out(T);
  out=0.0;
  size_t n=0;
  
  ifstream in(path);
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    for(size_t i=0;i<T;i++)
      for(size_t j=i;j<T;j++)
	{
	  double re,im;
	  if(!(in>>re>>im)) CRASH("Unable to read %s conf %zu, %zu %zu",path.c_str(),iconf,i,j);
	  size_t t=(T+i-j)%T;
	  out[t][iconf]+=re;
	  if(iconf==0) n++;
	}
  
  out.clusterize(1);
  
  out*=-1.0/n;
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

djvec_t load_bubble(const char *path)
{
  djvec_t out(T);
  djvec_t ave(2);
  ave=0.0;
  out=0.0;
  
  array<int,T> n;
  for(auto &nt : n) nt=0;
  
  ifstream in0(combine(path,0));
  ifstream in1(combine(path,1));
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    {
      array<complex<double>,T> d;
      for(size_t i=0;i<T;i++)
	{
	  double re0,im0;
	  double re1,im1;
	  if(!(in0>>re0>>im0)) CRASH("Unable to read %s conf0 %zu, %zu",combine(path,0).c_str(),iconf,i);
	  if(!(in1>>re1>>im1)) CRASH("Unable to read %s conf1 %zu, %zu",combine(path,1).c_str(),iconf,i);
	  double re=(re0-re1)/2;
	  double im=(im0-im1)/2;
	  d[i]=complex<double>(re,im);
	  ave[0][iconf]+=re;
	  ave[1][iconf]+=im;
	}
      
      for(size_t i=0;i<T;i++)
	for(size_t j=0;j<T;j++)
	  {
	    size_t t=(T+i-j)%T;
	    
	    out[t][iconf]+=-(d[i]*conj(d[j])).real();
	    if(iconf==0) n[t]++;
	  }
    }
  out.clusterize();
  ave.clusterize();
  ave/=T;
  cout<<ave.ave_err()<<endl;
  for(size_t t=0;t<T;t++) out[t]/=n[t];
  out+=sqr(ave[0])+sqr(ave[1]);
  
  out*=pow(L,3);
  
  double dum;
  if(in0>>dum) CRASH("Should have reached EOF, obtained %lg",dum);
  if(in1>>dum) CRASH("Should have reached EOF, obtained %lg",dum);
  
  return -out.symmetrized(1);
}

int main()
{
  set_njacks(nconfs);
  cout<<"NConfs: "<<nconfs<<endl;
  
  djvec_t manette=load("data/manette");
  manette.ave_err().write("plots/manette.xmg");
  
  djvec_t mesone00=load_LO("data/mesone00");
  mesone00.ave_err().write("plots/mesone00.xmg");
  effective_mass(mesone00).ave_err().write("plots/mesone_eff00.xmg");
  
  djvec_t mesone01=load_LO("data/mesone01");
  mesone01.ave_err().write("plots/mesone01.xmg");
  effective_mass(mesone01).ave_err().write("plots/mesone_eff01.xmg");
  
  djvec_t ratio=manette/mesone01;
  ratio.ave_err().write("plots/ratio.xmg");
  forward_derivative(ratio).ave_err().write("plots/ratio_der.xmg");
  
  djvec_t disc_PP=load_bubble("data/bubble_P5%d");
  disc_PP.ave_err().write("plots/disc_PP.xmg");

   djvec_t disc_VV=djvec_t(load_bubble("data/bubble_V1%d")+load_bubble("data/bubble_V2%d")+load_bubble("data/bubble_V3%d"))/3.0;
   disc_VV.ave_err().write("plots/disc_VV.xmg");

  
  // plot_ratio<<forward_derivative(ratio).ave_err()<<endl;
  // plot_ratio.new_data_set();
  
  // djack_t M,A,SL;
  // two_pts_with_ins_ratio_fit(M,A,SL,mesone01,manette,T/2,8,15,"plots/pi_mass.xmg","plots/pi_slope.xmg");
  
  
  return 0;
}
