#include <tranalisi.hpp>

#include <loadCorrCov.hpp>

const size_t T=224;
  
const string basePath=
  "/home/francesco/QCD/LAVORI/GM3/C.06.112/data/mix_fixed_l_s1_TM_";

djvec_t read(const string& tag)
{
  const auto [readT,corr,cov]=loadCorrCov(basePath+tag,true);
  
  if(readT!=T)
    CRASH("wrong read T");
  
  effective_mass(corr).ave_err().write("plots/"+tag+".xmg");
  
  return corr;
}

int main()
{
  set_njacks(50);
  
#define READ(TAG) const djvec_t TAG=read(#TAG)
  READ(V0V0);
  READ(VKVK);
  READ(A0A0);
  READ(AKAK);
#undef READ
  
  return 0;
}
