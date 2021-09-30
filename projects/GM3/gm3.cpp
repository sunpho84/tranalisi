#define EXTERN
#define INIT_EXTERN_TO(ARGS...) ARGS

#include <fit.hpp>

#include <GM3/aLaLuscherRepresentation.hpp>
#include <GM3/kernels.hpp>
#include <GM3/perens.hpp>

using namespace std;

void perens_t::loadDisco()
{
  ifstream file("disco.dat");
  if(not file.good())
    CRASH("Cannot open disco.dat");
  
  string conf;
  size_t t;
  double y;
  vector<double> l;
  while(file>>conf>>t>>y)
    l.push_back(y);
  
  if(l.size()%T!=0)
    CRASH("number of lines %zu is not multiple of T=%zu",l.size(),T);
  
  const size_t nDiscoConfs=
    l.size()/T;
  
  const double discoClustSize=
    nDiscoConfs/double(njacks);
  
  cout<<"nDiscoConfs: "<<nDiscoConfs<<endl;
  cout<<"discoClustSize: "<<discoClustSize<<endl;
  djvec_t d(T,0.0);
  for(size_t t=0;t<T;t++)
    jackknivesFill(nDiscoConfs,[this,&d,&l,t](const size_t& iConf,const size_t& iClust,const double& w)
    {
      d[t][iClust]+=l[t+T*iConf]*w;
    });
  d.clusterize(discoClustSize).symmetrize();
  
  d.ave_err().write("plots/disco.xmg");
  
  const djvec_t conn=
    getAveForRego(0,nSources,idVKVK,REGO_OS)*2// *5.0/9.0/4.0
    ;
  
  const djvec_t disc=
    d// /36.0
    ;
  
  const djvec_t tot=
    conn+disc;
  
  const djvec_t rat=
    disc/tot;
  
  effective_mass(conn).ave_err().write("/tmp/conn.xmg");
  effective_mass(disc).ave_err().write("/tmp/disc.xmg");
  effective_mass(tot).ave_err().write("/tmp/tot.xmg");
  rat.ave_err().write("/tmp/rat.xmg");
}

int main(int narg,char **arg)
{
  init(narg,arg);
  
  perens_t ens(".");
  
  ens.loadData(narg,arg);
  
  if(MPIrank==0)
    {
      ens.computeMPCAC();
      ens.Z=ens.determineRenoConst();
      
      if(ens.canUseRawData)
	{
	  ens.analyzeRawData();
	  ens.convertForSilvano();
	}
      
      if(file_exists("disco.dat"))
	ens.loadDisco();
      
      for(const RegoType& rego : {REGO_TM,REGO_OS})
	{
	  const auto [rep,repInfVol]=
	    ens.fitVKVK(rego,ens.nLevels);
	  
	ens.computeAmu(rego,rep,repInfVol);
	}
    }
  
  close();
  
  return 0;
}
