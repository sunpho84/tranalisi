#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <boot.hpp>
#include <oper.hpp>

ostream& operator<<(ostream &out,const ave_err_t &ae)
{
  if(!isnan(ae.ave) && !isnan(ae.err)) out<<ae.ave<<" "<<ae.err;
  
  return out;
}

double effective_mass(double ct,double ct_p_dt,int t,int TH,double guess,int par,int dt)
{
  double target=ct_p_dt/ct;
  return Brent_solve([TH,t,dt,par,target](double x){return corr_fun_effmass_ratio(x,TH,t,dt,par)-target;},guess);
}
