#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  //pars
  int T=-1,ext_njacks=-1,iel1=0,iel2=0;
  string path_out="";
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:n:a:b:o::"))!= -1)
    switch (c)
      {
      case 'T': T=to_int(optarg); break;
      case 'n': ext_njacks=to_int(optarg); break;
      case 'a': iel1=to_int(optarg); break;
      case 'b': iel2=to_int(optarg); break;
      case 'o': path_out=optarg; break;
      case '?': exit(0);break;
      default: CRASH("Unknown option -%c",optopt);
      }
  
  //check mandatory options
  if(T==-1) cerr<<"Missing argument T"<<endl;
  if(ext_njacks==-1) cerr<<"Missing argument n"<<endl;
  
  //parse paths
  string path_in1="/dev/stdin",path_in2="/dev/stdin";
  for(int i=optind;i<narg;i++)
    {
      if(i-optind==0) path_in1=arg[i];
      if(i-optind==1) path_in2=arg[i];
    }
  
  if(T==-1 or ext_njacks==-1) close_with_mess("Use: %s -T=size -n=njacks -a=iel1 path_in1 -b=iel2 path_in2",arg[0]);
  
  //set njacks and read
  set_njacks(ext_njacks);
  const djvec_t data1=read_djvec(path_in1,T,iel1);
  const djvec_t data2=read_djvec(path_in2,T,iel2);
  const djvec_t rat=data1/data2;
  
  if(path_out=="")
    cout<<rat.ave_err()<<endl;
  else
    rat.ave_err().write(path_out);
  
  return 0;
}
