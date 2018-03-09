#ifndef _ANALYSIS_HPP
#define _ANALYSIS_HPP

#ifndef EXTERN_ANALYSIS
 #define EXTERN_ANALYSIS extern
 #define INIT_ANALYSIS_TO(...)
#else
 #define INIT_ANALYSIS_TO(...) __VA_ARGS__
#endif

#include <MOM2/perens.hpp>

EXTERN_ANALYSIS map<string,perens_t> _data;

//! returns access to the data
perens_t& data(const string &key);

//! remove an entry
void data_erase(const string &key);

inline void compute_or_load_all(const string &grp_name)
{
  for(auto &path : pars::ens)
    {
      const string name=path+"_"+grp_name;
      data(name)
	.read_pars(path)
	.set_pars_for_scratch()
	.get_deltam_cr()
	.set_indices()
	.allocate()
	.read_or_compute();
    }
}

inline void plot_all_Z(const string &grp_name,const string &suffix)
{
  for(auto &path : pars::ens)
    {
      const string name=path+"_"+grp_name;
      data(name)
	.plot_Z(suffix);
    }
}

#define DEFINE_SINGLE_COMMAND_ALL(ALL_COMMAND,SINGLE_COMMAND)		\
  inline void ALL_COMMAND(const string &grp_name_out,const string &grp_name_in)	\
  {									\
    for(auto &path : pars::ens)						\
      {									\
	const string name_out=path+"_"+grp_name_out;			\
	string name_in=path+"_"+grp_name_in;				\
	data(name_out)=data(name_in).SINGLE_COMMAND();			\
      }									\
  }

DEFINE_SINGLE_COMMAND_ALL(average_all_r,average_r)
DEFINE_SINGLE_COMMAND_ALL(average_all_equiv_momenta,average_equiv_momenta)
DEFINE_SINGLE_COMMAND_ALL(val_chir_extrap_all,val_chir_extrap)

/////////////////////////////////////////////////////////////////

void average(const string out,const string in1,const string in2);

#endif
