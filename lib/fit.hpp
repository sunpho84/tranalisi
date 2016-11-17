#ifndef _FIT_HPP
#define _FIT_HPP

#include <grace.hpp>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>
#include <Math/MinimizerOptions.h>
#include <vector>

using namespace std;
using namespace ROOT::Minuit2;

//! set the level of verbosity
inline void set_printlevel(int lev)
{ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(lev);}

//! perform a fit to constant (uncorrelated)
template <class TV> auto constant_fit(const TV &in,int xmin,int xmax,string path="") -> typename TV::base_type
{
  using T=typename TV::base_type;
  T out(init_nel(in[0]));
  
  double norm=0;
  
  //take weighted average
  for(int iel=max(xmin,0);iel<=min(xmax,(int)in.size()-1);iel++)
    {
      auto ele=in[iel];
      double err=in[iel].err();
      double weight=1/(err*err);
      if(!std::isnan(err)&&err!=0)
        {
          out+=ele*weight;
          norm+=weight;
        }
    }
  
  //take simply average if error was zero
  if(norm==0)
    for(int iel=max(xmin,0);iel<=min(xmax,(int)in.size()-1);iel++)
      {
        norm++;
        out+=in[iel];
      }
  
  //normalize
  out/=norm;
  
  if(path!="")
    {
      grace_file_t plot(path);
      plot.write_constant_band(xmin,xmax,out,grace::RED);
      
      plot.no_line();
      plot.set_colors(grace::BLACK);
      plot<<in.ave_err();
    }
  
  return out;
}

#endif
