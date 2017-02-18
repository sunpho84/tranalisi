#ifndef _MINIMIZER_HPP
#define _MINIMIZER_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <vector>

//////////////////////////////////////////////// Minuit2 implementation //////////////////////////

#if MINIMIZER_TYPE == MINUIT2

#include <Math/MinimizerOptions.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

using namespace ROOT::Minuit2;
using namespace ROOT::Math;

//! wrapper against minimization class
class minimizer_t
{
 public:
  
};

//! set the level of verbosity
inline void set_printlevel(int lev)
{MnPrint::SetLevel(lev);}

#endif

#endif
