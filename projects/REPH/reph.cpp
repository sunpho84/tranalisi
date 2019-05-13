#include <tranalisi.hpp>

size_t L,T;
size_t nMass;
vector<double> mass;
size_t nMoms;
vector<array<double,3>> moms;

int main(int narg,char **arg)
{
  set_njacks(15);
  
  raw_file_t fin("jacks/input.txt","r");
  
  L=fin.read<size_t>("L");
  T=fin.read<size_t>("T");
  
  nMass=fin.read<size_t>("NMass");
  mass.resize(nMass);
  for(size_t iMass=0;iMass<nMass;iMass++)
    mass[iMass]=fin.read<double>();
  
  nMoms=fin.read<size_t>("NMoms");
  moms.resize(nMoms);
  for(size_t iMom=0;iMom<nMoms;iMom++)
    for(size_t mu=0;mu<3;mu++)
      moms[iMom][mu]=fin.read<double>();
  
  index_t ind({{"iks",nMass},
	       {"ikt",nMass},
	       {"moms",nMoms},
	       {"momt",nMoms},
	       {"gamma",1},
	       {"reim",2}});
  
  grace_file_t plot("plots/2pts.xmg");
  
  size_t iMs=0;
  size_t iMt=0;
  
  for(size_t iMoms=0;iMoms<nMoms;iMoms++)
    for(size_t iMomt=iMoms;iMomt<nMoms;iMomt++)
      {
	
	const size_t i=ind({iMs,iMt,iMoms,iMomt,0,0});
	
	const djvec_t corr=read_djvec("jacks/oPPo-ss",T,i);
	
	plot.write_vec_ave_err(effective_mass(corr.symmetrized()).ave_err());
      }
  
  return 0;
}
