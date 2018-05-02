#include <tranalisi.hpp>

const size_t T=48;

const string base="/home/francesco/QCD/LAVORI/GM2/romiti/A40.24_light/raw_data/jacks";
const size_t nm=1,nr=2;

index_t idx({{"m1",nm},{"m2",nm},{"r",nr}});

djvec_t get(const string base,const char si_TV,const char so_TV,const size_t im1,const size_t im2,const size_t par)
{
  djvec_t A(T);
  
  for(size_t si_123=1;si_123<=3;si_123++)
    for(size_t r=0;r<nr;r++)
      A+=read_djvec(combine("%s/corr_%c%zu%c%zu",base.c_str(),si_TV,si_123,so_TV,si_123),T,idx({im1,im2,r}));
  
  return A.symmetrized(par)/3/nr;
}

int main(int narg,char **arg)
{
  cout<<"Reading from: "<<base<<endl;
 
 set_njacks(15);
  
  const size_t im=0;
  djvec_t VKVK=get(base,'V','V',im,im,1);
  djvec_t VKTK=get(base,'V','T',im,im,-1);
  djvec_t TKVK=-get(base,'T','V',im,im,-1);
  djvec_t TKTK=-get(base,'T','T',im,im,1);
  
  VKVK.ave_err().write("/tmp/VKVK");
  VKTK.ave_err().write("/tmp/VKTK");
  TKVK.ave_err().write("/tmp/TKVK");
  TKTK.ave_err().write("/tmp/TKTK");
  
  /////////////////////////////////////////////////////////////////
  
  const size_t t0=3;
  const vector<djvec_t> eig=gevp({VKVK,VKTK,VKTK,TKTK},t0);
  
  eig[0].ave_err().write("/tmp/eig1");
  eig[1].ave_err().write("/tmp/eig2");
  
  effective_mass(eig[0]).ave_err().write("/tmp/eff_eig1");
  effective_mass(eig[1]).ave_err().write("/tmp/eff_eig2");
  
  djvec_t rat=effective_mass(eig[1])/effective_mass(eig[0]);
  rat.ave_err().write("/tmp/eff_rat");
  
  return 0;
}
