#ifndef _ANALYSIS_HPP
#define _ANALYSIS_HPP

#ifndef EXTERN_ANALYSIS
 #define EXTERN_ANALYSIS extern
 #define INIT_ANALYSIS_TO(...)
#else
 #define INIT_ANALYSIS_TO(...) __VA_ARGS__
#endif

#include <perens.hpp>

EXTERN_ANALYSIS map<string,perens_t> data;

inline void create_all(const string &prefix)
{
  for(auto &name : ens)
    data[prefix+name].create_from_scratch();
}

#endif
