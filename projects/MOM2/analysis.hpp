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

void compute_or_load_all_ingredients();

void plot_all_Z(const string &suffix);

#define DEFINE_SINGLE_SELF_COMMAND_ALL(ALL_COMMAND,SINGLE_COMMAND)	\
  inline void ALL_COMMAND()						\
  {									\
    for(auto &path : pars::ens)						\
      data(path,ASSERT_PRESENT).SINGLE_COMMAND();			\
  }

DEFINE_SINGLE_SELF_COMMAND_ALL(compute_deltam_from_prop_all,compute_deltam_from_prop)
// DEFINE_SINGLE_SELF_COMMAND_ALL(compute_Z_all,compute_Z)

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

// DEFINE_SINGLE_COMMAND_ALL(average_all_r,average_r)
DEFINE_SINGLE_COMMAND_ALL(average_all_equiv_momenta,average_equiv_momenta)
// DEFINE_SINGLE_COMMAND_ALL(val_chir_extrap_all,val_chir_extrap)

// /////////////////////////////////////////////////////////////////

// //! average in1 and in2 to form out, removing the in
// void average(const string out,const string in1,const string in2);

// //! in1/in2-1 to form out, removing the in
// void ratio_minus_one(const string out,const string in1,const string in2);

// //! print all ensembles available
// void list_ensembles();

// //! extrapolate w.r.t sea mass
// void sea_chir_extrap(const string out,const vector<string> &ens_list);

// //! print momenta with discretization
// void print_discr();

#endif
