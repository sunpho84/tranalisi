#ifndef _MINIMIZER_HPP
#define _MINIMIZER_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string>
#include <vector>

using namespace std;

//////////////////////////////////////////////// Minuit2 implementation //////////////////////////

#if MINIMIZER_TYPE == MINUIT2

// #include <Math/MinimizerOptions.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

using minimizer_base_t=ROOT::Minuit2::FCNBase;

//! wrapper against minimization pars
class minimizer_pars_t
{
public:
  ROOT::Minuit2::MnUserParameters pars;
  
  //! add a parameter
  void add(const string &name,double val,double err)
  {pars.Add(name,val,err);}
  
  //! set a parameter
  void set(int ipar,double val)
  {pars.SetValue(ipar,val);}
  
  //! fix a parameter
  void fix(int ipar)
  {pars.Fix(ipar);}
  
  //! fix a parameter to a given value
  void fix_to(int ipar,double val)
  {
    set(ipar,val);
    fix(ipar);
  }
  
  //! return the number of parameters
  size_t size()
  {return pars.Parameters().size();}
};

//! wrapper against minimization class
class minimizer_t
{
  ROOT::Minuit2::MnMigrad migrad;
  minimizer_t();
  
 public:
  //construct from migrad
  template <class function_t> minimizer_t(const function_t &fun,const minimizer_pars_t &pars) : migrad(fun,pars.pars) {}
  
  //! minimizer
  vector<double> minimize()
  {
    //call minimizer
    ROOT::Minuit2::FunctionMinimum min=migrad();
    ROOT::Minuit2::MinimumParameters par_min=min.Parameters();
    
    //copy the result
    size_t size=par_min.Vec().size();
    vector<double> pars(size);
    for(size_t it=0;it<size;it++) pars[it]=par_min.Vec()[it];
    
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

#endif
