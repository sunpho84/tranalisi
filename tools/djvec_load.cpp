#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  //pars
  int T=-1,ext_njacks=-1,iel=0,parity=-2;
  bool append=false;
  string bin_path="";
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:ab:n:p:i::"))!= -1)
    {
      printf("%c\n",c);
      switch (c)
      {
      case 'T': T=to_int(optarg); break;
      case 'a': append=true; break;
      case 'b': bin_path=optarg;break;
      case 'n': ext_njacks=to_int(optarg); break;
      case 'p': parity=to_int(optarg); break;
      case 'i': iel=to_int(optarg); break;
      case '?': exit(0);break;
      default: CRASH("Unknown option -%c",optopt);
      }
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
  
  if(T==-1 or ext_njacks==-1 or path_in=="") close_with_mess("Use: %s -T=size -n=njacks -i=iel[0] path_in[stdin] -a -p=parity -b=bin_path path_out[stdout]",arg[0]);
  
  //put a warning
  if(path_in=="/dev/stdin") cerr<<"Reading from stdin"<<endl;
  
  //set njacks and read
  set_njacks(ext_njacks);
  djvec_t data=read_djvec(path_in,T,iel);
  
  if(parity!=-2)
    data=data.symmetrized(parity);
  
  //write average and error
  grace_file_t out(path_out);
  out.write_vec_ave_err(data.ave_err(),grace::RED,grace::SQUARE);
  
  if(bin_path!="")
    {
      raw_file_t out(bin_path,append?"a":"w");
      out.bin_write(data);
    }
  
  return 0;
}
