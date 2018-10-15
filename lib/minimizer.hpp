#ifndef _MINIMIZER_HPP
#define _MINIMIZER_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_MINIMIZER
 #define EXTERN_MINIMIZER extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#include <string>
#include <vector>

using namespace std;

//////////////////////////////////////////////// Minuit1 implementation //////////////////////////

#if MINIMIZER_TYPE == MINUIT

#include <cify.hpp>
#include <TMinuit.h>

//! set the level of verbosity
EXTERN_MINIMIZER int glb_print_level;
inline void set_printlevel(int level)
{glb_print_level=level;}

//! host a single par
class par_t
{
 public:
  string name;
  double val;
  double err;
  bool is_fixed;
  
  par_t(const string &name,double val,double err) : name(name),val(val),err(err),is_fixed(false) {}
};

//! base class
class minimizer_fun_t
{
public:
  virtual ~minimizer_fun_t() {}
  virtual double operator()(const std::vector<double>& x) const=0;
  virtual double Up() const=0;
};

//! wrapper against minimization pars
class minimizer_pars_t
{
public:
  vector<par_t> pars;
  
  //! add a parameter
  void add(const string &name,double val,double err)
  {pars.push_back(par_t(name,val,err));}
  
  //! set a parameter
  void set(int ipar,double val)
  {pars[ipar].val=val;}
  
  //! fix a parameter
  void fix(int ipar)
  {pars[ipar].is_fixed=true;}
  
  //! fix a parameter to a given value
  void fix_to(int ipar,double val);
  
  //! return the number of parameters
  size_t size() const
  {return pars.size();}
  
  //! return the number of parameters
  size_t nfree_pars() const
  {
    size_t n=0;
    for(size_t i=0;i<pars.size();i++) if(not pars[i].is_fixed) n++;
    return n;
  }
};

//! pointer used by fcn
EXTERN_MINIMIZER const minimizer_fun_t *fun_ptr;
//! number of parameters
EXTERN_MINIMIZER int fun_npars;
//! wrapper to the function
EXTERN_MINIMIZER void fcn(int &npar,double *fuf,double &ch,double *p,int flag);

//! wrapper against minimization class
class minimizer_t
{
  //! function to be minimized
  const minimizer_fun_t &fun;
  
  //! root minimizer
  TMinuit minu;
  minimizer_t();
  
 public:
  //construct from migrad
  minimizer_t(const minimizer_fun_t &fun,const minimizer_pars_t &pars) : fun(fun)
  {
    //set verbosity
    minu.SetPrintLevel(glb_print_level);
    
    //set tolerance
    //double tol=1e-16;
    //int iflag;
    minu.SetMaxIterations(10000000);
    //minu.mnexcm("SET ERR",&tol,1,iflag);
    
    //set the external pointer
    fun_ptr=&fun;
    minu.SetFCN(fcn);
    
    //set the pars
    fun_npars=pars.size();
    for(size_t ipar=0;ipar<pars.size();ipar++)
      {
	minu.DefineParameter(ipar,pars.pars[ipar].name.c_str(),pars.pars[ipar].val,pars.pars[ipar].err,0.0,0.0);
	if(pars.pars[ipar].is_fixed) minu.FixParameter(ipar);
      }
  }
  
  //! minimizer
  vector<double> minimize()
  {
    int npars=minu.GetNumPars();
    vector<double> pars(npars);
    minu.Migrad();
    double dum;
    for(int it=0;it<npars;it++)
      minu.GetParameter(it,pars[it],dum);
    return pars;
  }
  
  //! evaluate the function
  double eval(const vector<double> &pars)
  {
    double res;
    double in_pars[pars.size()];
    for(size_t it=0;it<pars.size();it++) in_pars[it]=pars[it];
    minu.Eval(pars.size(),NULL,res,in_pars,0);
    return res;
  }
};

#endif

//////////////////////////////////////////////// Minuit2 implementation ///////////////////////////////////////////

#if MINIMIZER_TYPE == MINUIT2

// #include <Math/MinimizerOptions.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

using minimizer_fun_t=ROOT::Minuit2::FCNBase;

//! wrapper against minimization pars
class minimizer_pars_t
{
public:
  ROOT::Minuit2::MnUserParameters pars;
  
  //! add a parameter
  void add(const string &name,double val,double err)
  {pars.Add(name,val,err);}
  
  void setlimits(const string &name,double min,double max)
  {pars.SetLimits(name,min,max);}
  
  //! set a parameter
  void set(int ipar,double val)
  {pars.SetValue(ipar,val);}
  
  //! fix a parameter
  void fix(int ipar)
  {pars.Fix(ipar);}
  
  //! fix a parameter to a given value
  void fix_to(int ipar,double val);
  
  //! return the number of parameters
  size_t size() const
  {return pars.Parameters().size();}
  
  //! return the number of parameters
  size_t nfree_pars() const
  {
    size_t n=0;
    for(size_t i=0;i<size();i++) if(not pars.Parameter(i).IsFixed()) n++;
    return n;
  }
};

//! wrapper against minimization class
class minimizer_t
{
  ROOT::Minuit2::MnMigrad migrad;
  size_t npars;
  minimizer_t();
  
 public:
  //construct from migrad
  minimizer_t(const minimizer_fun_t &fun,const minimizer_pars_t &pars) : migrad(fun,pars.pars),npars(pars.size()) {}
  
  //! minimizer
  vector<double> minimize()
  {
    //call minimizer
    ROOT::Minuit2::FunctionMinimum min=migrad();
    ROOT::Minuit2::MnUserParameterState par_state=min.UserState();
    if(!min.IsValid())
      {
	cerr<<"WARNING! Minimization failed"<<endl;
	cerr<<"Has accurate cov: "<<min.HasAccurateCovar()<<endl;
	cerr<<"Has positive definite cov: "<<min.HasPosDefCovar()<<endl;
	cerr<<"Has valid covariance"<<min.HasValidCovariance()<<endl;
	cerr<<"Has valid parameters: "<<min.HasValidParameters()<<endl;
	cerr<<"Has reached call limit: "<<min.HasReachedCallLimit()<<endl;
	cerr<<"Hesse failed: "<<min.HesseFailed()<<endl;
	cerr<<"Has cov: "<<min.HasCovariance()<<endl;
	cerr<<"Is above max edm: "<<min.IsAboveMaxEdm()<<endl;
	cout<<"FunctionCalls: "<<migrad.NumOfCalls()<<endl;
	cerr<<par_state<<endl;
      }
    
    //copy the result
    vector<double> pars(npars);
    for(size_t ipar=0;ipar<npars;ipar++) pars[ipar]=migrad.Value(ipar);
    
    return pars;
  }
  
  //! evaluate the function
  double eval(const vector<double> &pars)
  {return migrad().Fval();}
};

//! set the level of verbosity
inline void set_printlevel(int lev)
{ROOT::Minuit2::MnPrint::SetLevel(lev);}

#endif

inline void minimizer_pars_t::fix_to(int ipar,double val)
{
  set(ipar,val);
  fix(ipar);
}

#undef INIT_TO
#undef EXTERN_MINIMIZER

#endif
