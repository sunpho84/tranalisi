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
#define ASSERT_PRESENT true
#define PRESENCE_NOT_NEEDED false

//! returns access to the data
perens_t& data(const string &key,const bool assert_present_flag);

//! remove an entry
void data_erase(const string &key);

inline void compute_or_load_all()
{
  cout<<"Going to rotate propagators: "<<(pars::twisted_run and pars::phys_basis)<<endl;
  
  for(auto &name : pars::ens)
    data(name,PRESENCE_NOT_NEEDED)
      .read_pars(name)
      .set_pars_for_scratch()
      .set_indices()
      .get_deltam()
      .get_mPCAC()
      .get_meson_mass()
      .allocate()
      .read_or_compute();
}

inline void plot_all_Z(const string &suffix)
{
  for(auto &path : pars::ens)
    {
      const string name=path;
      data(name,ASSERT_PRESENT)
	.plot_Z(suffix);
    }
}

#define DEFINE_SINGLE_COMMAND_ALL(ALL_COMMAND,SINGLE_COMMAND)		\
  inline void ALL_COMMAND()						\
  {									\
    for(auto &path : pars::ens)						\
      {									\
	const string name_out=path;					\
	string name_in=path;						\
	data(name_out,PRESENCE_NOT_NEEDED)=				\
	  data(name_in,ASSERT_PRESENT).SINGLE_COMMAND();		\
      }									\
  }

DEFINE_SINGLE_COMMAND_ALL(average_all_r,average_r)
DEFINE_SINGLE_COMMAND_ALL(average_all_equiv_momenta,average_equiv_momenta)
DEFINE_SINGLE_COMMAND_ALL(val_chir_extrap_all,val_chir_extrap)

/////////////////////////////////////////////////////////////////

//! average in1 and in2 to form out, removing the in
void average(const string out,const string in1,const string in2);

//! print all ensembles available
void list_ensembles();

//! extrapolate w.r.t sea mass
void sea_chir_extrap(const string out,const vector<string> &ens_list);

//! print momenta with discretization
void print_discr();

#endif
