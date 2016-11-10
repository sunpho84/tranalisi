#ifndef _FIT_HPP
#define _FIT_HPP

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>
#include <Math/MinimizerOptions.h>

using namespace ROOT::Minuit2;

//! set the level of verbosity
void set_print_level(int lev);

#endif
