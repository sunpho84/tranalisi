#include <tranalisi.hpp>

const size_t T=48,TH=T/2;

const string base="/home/francesco/QCD/LAVORI/GM2/romiti/A40.24/raw_data/jacks";
const size_t nm=3,nr=2;

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
  
  typedef Matrix2d Matr;
  
  GeneralizedEigenSolver<Matr> ges;
  
  djvec_t eig1(TH+1),eig2(TH+1);
  
  const size_t t0=3;
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      Matr b;
      
      b(0,0)=VKVK[t0][ijack];
      b(0,1)=VKTK[t0][ijack];
      b(1,0)=b(0,1);
      b(1,1)=TKTK[t0][ijack];
      
      for(size_t t=0;t<=TH;t++)
	{
	  Matr a;
	  
	  a(0,0)=VKVK[t][ijack];
	  a(0,1)=VKTK[t][ijack];
	  a(1,0)=a(0,1);
	  a(1,1)=TKTK[t][ijack];
	  
	  ges.compute(a,b);
	  
	  eig1[t][ijack]=ges.eigenvalues()(0).real();
	  eig2[t][ijack]=ges.eigenvalues()(1).real();
	}
    }
  
  eig1.ave_err().write("/tmp/eig1");
  eig2.ave_err().write("/tmp/eig2");
  
  effective_mass(eig1).ave_err().write("/tmp/eff_eig1");
  effective_mass(eig2).ave_err().write("/tmp/eff_eig2");
  
  djvec_t rat=effective_mass(eig2)/effective_mass(eig1);
  rat.ave_err().write("/tmp/eff_rat");
  
  return 0;
}
