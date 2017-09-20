#include <tranalisi.hpp>

size_t conf_min=1000,conf_max=1140,conf_each=10;
size_t nconfs=(conf_max-conf_min)/conf_each+1;
range_t range({conf_min,conf_each,conf_max});

size_t T=48,TH=T/2;
size_t tmin=12,tmax=20;

const size_t nlevs=5;
const size_t levs[nlevs]={0,1,4,9,16};

djvec_t read(size_t ilev,size_t jlev,size_t klev,size_t llev)
{
  const size_t ncols_tot=2;
  const vector<size_t> cols={0};
  return read_conf_set_t("%04d/"+combine("mes_contr_%03zu_%03zu_LS_%03zu_%03zu",levs[ilev],levs[jlev],levs[klev],levs[llev]),range,ncols_tot,cols).symmetrized();
}

int main(int narg,char **arg)
{
  set_njacks(nconfs);
  
  for(size_t ilev=0;ilev<nlevs;ilev++)
    for(size_t jlev=0;jlev<nlevs;jlev++)
      for(size_t klev=0;klev<nlevs;klev++)
	for(size_t llev=0;llev<nlevs;llev++)
	  {
	    djvec_t corr=read(ilev,jlev,klev,llev);
	    (corr).ave_err().write(combine("plots/%03zu_%03zu_LS_%03zu_%03zu.xmg",levs[ilev],levs[jlev],levs[klev],levs[llev]));
	  }
  
  return 0;
}
