#include <tranalisi.hpp>

size_t conf_min=1000,conf_max=1140,conf_each=10;
size_t nconfs=(conf_max-conf_min)/conf_each+1;
range_t range({conf_min,conf_each,conf_max});

size_t T=48,TH=T/2;
size_t tmin=12,tmax=20;

int main(int narg,char **arg)
{
  set_njacks(nconfs);
  
  input_file_t input("input");
  size_t lev_min=0;
  size_t lev_max=input.read<size_t>("LevMax");
  size_t lev_each=input.read<size_t>("LevEach");
  size_t nlevs=(lev_max-lev_min)/lev_each+1;
  
  size_t ncols_tot=2;
  vector<size_t> cols={0};
  
  djvec_t Z2(nlevs),M(nlevs);
  djvec_t corr[nlevs];
  djvec_t exc(nlevs);
  
  vector<double> levs(nlevs);
  
  for(size_t ilev=0;ilev<nlevs;ilev++)
    {
      size_t lev=lev_min+ilev*lev_each;
      levs[ilev]=lev;
      
      //load
      corr[ilev]=read_conf_set_t("%04d/"+combine("mes_contr_K%03zu",lev),range,ncols_tot,cols).symmetrized();
      
      //fit
      two_pts_fit(Z2[ilev],M[ilev],corr[ilev],TH,tmin,tmax,combine("plots/effmass_%03d.xmg",lev));
      
      //subtract contribution from ground state
      exc[ilev]=0.0;
      for(size_t t=1;t<=tmin;t++)
	exc[ilev]+=abs(corr[ilev][t]/two_pts_corr_fun(Z2[ilev],M[ilev],TH,t,1)-1.0);
    }
  
  //plot Z2 and M as a function of nlevs
  Z2.ave_err().write(levs,"plots/Z2.xmg");
  M.ave_err().write(levs,"plots/M.xmg");
  exc.ave_err().write(levs,"plots/exc.xmg");
  
  return 0;
}
