#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  //pars
  int T=-1,ext_njacks=-1,iel=0,par=1,tmin=-1,tmax=-1;
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:n:m:M:i::p::"))!= -1)
    switch(c)
    {
      case 'T': T=to_int(optarg); break;
      case 'n': ext_njacks=to_int(optarg); break;
      case 'i': iel=to_int(optarg); break;
      case 'p': par=to_int(optarg); break;
      case 'm': tmin=to_int(optarg); break;
      case 'M': tmax=to_int(optarg); break;
      case '?': exit(0);break;
      default: CRASH("Unknown option -%c",optopt);
      }
  
  //check mandatory options
  if(T==-1) cerr<<"Missing argument T"<<endl;
  if(ext_njacks==-1) cerr<<"Missing argument n"<<endl;
  
  //parse paths
  string path_in="/dev/stdin",path_out="/dev/stdout";
  for(int i=optind;i<narg;i++)
    {
      if(i-optind==0) path_in=arg[i];
      if(i-optind==1) path_out=arg[i];
    }
  
  if(T==-1 or ext_njacks==-1 or path_in=="" or tmin==-1 or tmax==-1)
    close_with_mess("Use: %s -T=size -n=njacks -i=iel[0] -m=tmin -M=tmax -p=par[1] path_in[stdin] path_out[stdout]",arg[0]);
  
  //put a warning
  if(path_in=="/dev/stdin") cerr<<"Reading from stdin"<<endl;
  
  //set njacks and read
  set_njacks(ext_njacks);
  djvec_t data=read_djvec(path_in,T,iel);
  
  const size_t TH=T/2;
  
  //write average and error
  cout.precision(16);
  cout<<constant_fit(effective_mass(data.symmetrized(par),TH,par),tmin,tmax,path_out)<<endl;
  
  djack_t Z2,M;
  data=data.symmetrized(par);
  two_pts_fit(Z2,M,data,TH,tmin,tmax,"","",par);
  cout<<"Z: "<<djack_t(sqrt(sqrt(Z2*Z2))).ave_err()<<endl;

  djvec_t ecc(TH+1);
  for(size_t t=0;t<=TH;t++)
    ecc[t]=data[t]-two_pts_corr_fun(Z2,M,TH,t,par);
  
  effective_mass(ecc,TH,par).ave_err().write("/tmp/exc");
  
  return 0;
}
