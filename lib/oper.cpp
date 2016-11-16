#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <boot.hpp>
#include <oper.hpp>
#include <tools.hpp>

const double NaN=numeric_limits<double>::quiet_NaN();

ostream& operator<<(ostream &out,const ave_err_t &ae)
{
  if(ae.is_printable()) out<<ae.ave<<" "<<ae.err;
  
  return out;
}

double effective_mass(double ct,double ct_p_dt,int t,int TH,double guess,int par,int dt)
{
  //compute the target, c(t+dt)/c(t)
  double target=ct_p_dt/ct;
  if(target<=0) return NaN;
  
  //function to minimize
  auto fun=[TH,t,dt,par,target](double x){return corr_fun_effmass_ratio(x,TH,t,dt,par)-target;};
  
  //guess must be positive
  if(guess<=0 or !isfinite(guess))
    {
#ifdef DEBUG
      cout<<"Guess "<<guess<<" invalid, setting it to 1"<<endl;
#endif
      guess=1;
    }
  
  //define reasonbale limits
  double lims[2]={guess,guess};
  double funs[2]={fun(guess),fun(guess)};
  bool isdiff,bracketted;
  double incr_fact=1.001;
  do
    {
      double trial_lims[2]={lims[0]/incr_fact,lims[1]*incr_fact};
      double trial_funs[2]={fun(trial_lims[0]),fun(trial_lims[1])};
      incr_fact*=incr_fact;
      
      //check limits are valid and changed
      isdiff=false;
      for(int il=0;il<2;il++)
	if(isfinite(trial_lims[il]) and isfinite(trial_funs[il]) and trial_lims[il]!=0 and trial_lims[il]>epsilon)
	  {
	    isdiff=true;
	    lims[il]=trial_lims[il];
	    funs[il]=trial_funs[il];
	  }
#ifdef DEBUG
	else
	  cout<<"Limit "<<il<<" corresponding to "<<trial_lims[il]<<" not finte: "<<trial_funs[il]<<endl;
      #endif

      bracketted=!same_sign(funs[0],funs[1]);
#ifdef DEBUG
      cout<<"Isdiff: "<<isdiff<<", increasing range, searching target "<<target<<" in range: ["<<lims[0]<<", "<<lims[1]<<"] ,extreme values: "<<funs[0]<<", "<<funs[1]<<endl;
#endif
    }
  while(!bracketted and isdiff);
  
  //return nan if it was not possible to set the problem
  if(!bracketted) return NaN;
  
  return Brent_solve(fun,lims[0],lims[1]);
}
