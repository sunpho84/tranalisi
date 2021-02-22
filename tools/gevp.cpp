#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  //pars
  int T=-1,ext_njacks=-1,iel1=0,iel2=0,iel3=0,t0=-1;
  string path_out="";
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"T:n:a:b:c:o:0::"))!= -1)
    switch (c)
      {
      case 'T': T=to_int(optarg); break;
      case 'n': ext_njacks=to_int(optarg); break;
      case 'a': iel1=to_int(optarg); break;
      case 'b': iel2=to_int(optarg); break;
      case 'c': iel3=to_int(optarg); break;
      case '0': t0=to_int(optarg); break;
      case 'o': path_out=optarg; break;
      case '?': exit(0);break;
      default: CRASH("Unknown option -%c",optopt);
      }
  
  //check mandatory options
  if(T==-1) CRASH("Missing argument T");
  if(t0==-1) CRASH("Missing argument 0");
  if(ext_njacks==-1) CRASH("Missing argument n");
  
  cout<<"T: "<<T<<endl;
  cout<<"t0: "<<t0<<endl;
  
  //parse paths
  vector<string> path_in(3);
  if(narg-optind!=3)
    CRASH("Please specify the paths");
  
  for(int i=optind;i<narg;i++)
    {
      const int j=i-optind;
      cout<<"Reading path "<<j<<" from: "<<arg[i]<<endl;
      path_in[j]=arg[i];
    }
  
 set_njacks(ext_njacks);
  
  djvec_t c00=read_djvec(path_in[0],T,iel1).symmetrized();
  djvec_t c01=read_djvec(path_in[1],T,iel2).symmetrized();
  djvec_t c11=read_djvec(path_in[2],T,iel3).symmetrized();
  
  c00.ave_err().write("plots/c00.xmg");
  c01.ave_err().write("plots/c01.xmg");
  c11.ave_err().write("plots/c11.xmg");
  
  effective_mass(c00).ave_err().write("plots/eff_c00.xmg");
  effective_mass(c01).ave_err().write("plots/eff_c01.xmg");
  effective_mass(c11).ave_err().write("plots/eff_c11.xmg");
  
  /////////////////////////////////////////////////////////////////
  
  vector<djvec_t> eig;
  vector<djvec_t> recastEigvec;
  vector<djvec_t> origEigvec;
  
  tie(eig,recastEigvec,ignore)=gevp({c00,c01,c01,c11},t0);
  
  eig[0].ave_err().write("plots/eig1.xmg");
  eig[1].ave_err().write("plots/eig2.xmg");
  
  effective_mass(eig[0]).ave_err().write("plots/eff_eig1.xmg");
  effective_mass(eig[1]).ave_err().write("plots/eff_eig2.xmg");
  
  djvec_t rat=effective_mass(eig[1])/effective_mass(eig[0]);
  rat.ave_err().write("plots/eff_rat.xmg");
  
  for(auto r : recastEigvec)
    cout<<r.ave_err()<<endl;
  
  return 0;
}
