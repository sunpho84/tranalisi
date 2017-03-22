#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  //pars
  int T=-1,par=1;
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:n:i::p::"))!= -1)
    switch (c)
      {
      case 'T': T=to_int(optarg); break;
      case 'p': par=to_int(optarg); break;
      case '?': exit(0);break;
      default: CRASH("Unknown option -%c",optopt);
      }
  
  //check mandatory options
  if(T==-1) cerr<<"Missing argument T"<<endl;
  
  //parse paths
  string path_in="/dev/stdin",path_out="/dev/stdout";
  for(int i=optind;i<narg;i++)
    {
      if(i-optind==0) path_in=arg[i];
      if(i-optind==1) path_out=arg[i];
    }
  
  if(T==-1 or path_in=="") close_with_mess("Use: %s -T=size -p=par[1] path_in[stdin] path_out[stdout]",arg[0]);
  
  //put a warning
  if(path_in=="/dev/stdin") cerr<<"Reading from stdin"<<endl;
  
  vector<double> data(T);
  ifstream file_in(path_in);
  if(!file_in.good()) CRASH("Opening %s to read",path_in.c_str());
  for(int t=0;t<T;t++) if(!(file_in>>data[t])) CRASH("Reading t=%zu",t);
  
  ofstream file_out(path_out);
  if(!file_out.good()) CRASH("Opening %s to write",path_out.c_str());
  for(int t=0;t<=T/2;t++) file_out<<t+1<<" "<<effective_mass(data[t],data[t+1],t,T/2,par)<<endl;
  
  return 0;
}
