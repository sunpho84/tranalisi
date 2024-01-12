#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  set_njacks(15);
  
  //pars
  int T=-1,par=1;
  int tmin=-1,tmax=-1;
  double mu=0;
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:m:M:n:i::u:p::"))!= -1)
    switch (c)
      {
      case 'T':
	T=to_int(optarg);
	cout<<"parsing T="<<T<<endl;
	break;
      case 'm':
	tmin=to_int(optarg);
	cout<<"parsing tmin="<<tmin<<endl;
	break;
      case 'M':
	tmax=to_int(optarg);
	cout<<"parsing tmax="<<tmax<<endl;
	break;
      case 'u':
	mu=strtod(optarg,nullptr);
	cout<<"Parsing mu= "<<mu<<endl;
	break;
      case 'p':
	par=to_int(optarg);
	cout<<"Parsing par= "<<mu<<endl;
	break;
      default: CRASH("Unknown option -%c",optopt);
      }
  
  //check mandatory options
  if(T==-1) cerr<<"Missing argument T"<<endl;
  if(tmin==-1) cerr<<"Missing argument tmin"<<endl;
  if(tmax==-1) cerr<<"Missing argument tmax"<<endl;
  
  //parse paths
  string path_in="/dev/stdin",path_out="/dev/stdout";
  for(int i=optind;i<narg;i++)
    {
      if(i-optind==0) path_in=arg[i];
      if(i-optind==1) path_out=arg[i];
    }
  
  if(T==-1 or path_in=="" or tmin==-1 or tmax==-1) close_with_mess("Use: %s -T=size -m=tmin -M=tmax -p=par[1] path_in[stdin] path_out[stdout]",arg[0]);
  
  //put a warning
  if(path_in=="/dev/stdin") cerr<<"Reading from stdin"<<endl;
  
  vector<double> ave(T/2+1),err(T/2+1);
  djvec_t data(T/2+1);
  ifstream file_in(path_in);
  if(!file_in.good()) CRASH("Opening %s to read",path_in.c_str());
  for(int t=0;t<T/2+1;t++)
    if(!(file_in>>ave[t]>>err[t]))
      CRASH("Reading t=%d",t);
    else
      data[t].fill_gauss({ave[t],err[t],3224+t});
  
  djack_t Z2,M;
  two_pts_fit(Z2,M,data,T/2,tmin,tmax,path_out,"",par);
  
  cout<<"M: "<<M.ave_err()<<endl;
  cout<<"Z: "<<sqrt(Z2).ave_err()<<endl;
  
  if(mu)
    {
      djack_t f=2*mu*sqrt(Z2)/(2*sqr(M));
      cout<<"f: "<<f.ave_err()<<endl;
    }
  
  // effective_mass(data,T/2,par).ave_err().write(path_out);
  
  return 0;
}
