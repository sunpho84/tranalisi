#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  //pars
  int T=-1,ext_njacks=-1,iel=0,par=1;
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:n:i::p::"))!= -1)
    switch (c)
      {
      case 'T': T=to_int(optarg); break;
      case 'n': ext_njacks=to_int(optarg); break;
      case 'i': iel=to_int(optarg); break;
      case 'p': par=to_int(optarg); break;
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
  
  if(T==-1 or ext_njacks==-1 or path_in=="") close_with_mess("Use: %s -T=size -n=njacks -i=iel[0] -p=par[1] path_in[stdin] path_out[stdout]",arg[0]);
  
  //put a warning
  if(path_in=="/dev/stdin") cerr<<"Reading from stdin"<<endl;
  
  //set njacks and read
  set_njacks(ext_njacks);
  djvec_t data=read_djvec(path_in,T,iel);
  
  //write average and error
  grace_file_t out(path_out);
  out.write_vec_ave_err(effective_mass(data.symmetrized(par),T/2,par).ave_err(),grace::RED,grace::SQUARE);
  
  return 0;
}
