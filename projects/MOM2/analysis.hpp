#ifndef _ANALYSIS_HPP
#define _ANALYSIS_HPP

#ifndef EXTERN_ANALYSIS
 #define EXTERN_ANALYSIS extern
 #define INIT_ANALYSIS_TO(...)
#else
 #define INIT_ANALYSIS_TO(...) __VA_ARGS__
#endif

#include <MOM2/perens.hpp>

EXTERN_ANALYSIS map<string,perens_t> data;

inline void compute_or_load_all(const string &grp_name)
{
  for(auto &path : pars::ens)
    {
      string name=path+"_"+grp_name;
      data[name]
	.read_pars(path)
	.set_indices()
	.allocate()
	.read_or_compute();
    }
}

#endif
