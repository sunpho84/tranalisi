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

//! store whether ingredients are computed
EXTERN_ANALYSIS bool ingredients_computed INIT_ANALYSIS_TO({false});

//! marks that ingredients are computed
inline void validate_ingredients()
{
  ingredients_computed=true;
}

//! ask permission to read ingredients
inline void needs_to_read_ingredients()
{
  if(not ingredients_computed) CRASH("ingredients not valid");
}

//! marks that ingredients are not computed
inline void invalidate_ingredients()
{
  ingredients_computed=false;
}

//! store whether Z are computed
EXTERN_ANALYSIS bool Z_computed INIT_ANALYSIS_TO({false});

//! marks that Z are computed
inline void validate_Z()
{
  Z_computed=true;
}

//! ask permission to read Z
inline void needs_to_read_Z()
{
  if(not Z_computed) CRASH("Z not valid");
}

//! marks that Z are not computed
inline void invalidate_Z()
{
  Z_computed=false;
}

/////////////////////////////////////////////////////////////////

//! returns access to the data
perens_t& data(const string &key,const bool assert_present_flag);

//! add an entry
void add_ens(const string &name);
//! remove an entry
void data_erase(const string &key);

void compute_or_load_all_ingredients();

void plot_all_Z(const string &suffix);

#define DEFINE_SINGLE_SELF_COMMAND_ALL(ALL_COMMAND,SINGLE_COMMAND,CHECK,CLAUSE) \
  inline void ALL_COMMAND()						\
  {									\
    CHECK();								\
    for(auto &path : pars::ens)						\
      data(path,ASSERT_PRESENT).SINGLE_COMMAND();			\
    CLAUSE();								\
  }

DEFINE_SINGLE_SELF_COMMAND_ALL(recompute_deltam_all,recompute_deltam,needs_to_read_ingredients,invalidate_Z)
DEFINE_SINGLE_SELF_COMMAND_ALL(compute_Z_all,compute_Z,needs_to_read_ingredients,validate_Z)

#define DEFINE_SINGLE_COMMAND_ALL(ALL_COMMAND,SINGLE_COMMAND,CHECK,CLAUSE) \
  inline void ALL_COMMAND()						\
  {									\
    CHECK();								\
    for(auto &path : pars::ens)						\
      {									\
	const string name_out=path;					\
 	string name_in=path;						\
 	data(name_out,PRESENCE_NOT_NEEDED)=				\
 	  data(name_in,ASSERT_PRESENT).SINGLE_COMMAND();		\
      }									\
    CLAUSE();								\
  }

DEFINE_SINGLE_COMMAND_ALL(average_all_r,average_r,needs_to_read_ingredients,invalidate_Z)
DEFINE_SINGLE_COMMAND_ALL(average_all_equiv_momenta,average_equiv_momenta,needs_to_read_ingredients,invalidate_Z)
DEFINE_SINGLE_COMMAND_ALL(val_chir_extrap_all,val_chir_extrap,needs_to_read_ingredients,invalidate_Z)
DEFINE_SINGLE_COMMAND_ALL(interpolate_to_p2ref_all,interpolate_to_p2ref,needs_to_read_Z,invalidate_ingredients)

// /////////////////////////////////////////////////////////////////

//! average Z in1 and in2 to form out, removing the in
void average_Z(const string out,const string in1,const string in2);

//! average ingredients in1 and in2 to form out, removing the in
void average_ingredients(const string out,const string in1,const string in2);

//! in1/in2-1 to form out, removing the in
void ratio_Z_minus_one(const string out,const string in1,const string in2);

// //! print all ensembles available
void list_ensembles();

//! extrapolate w.r.t sea mass
void sea_chir_extrap(const string out,const vector<string> &ens_list);

//! print momenta with discretization
void print_discr();

//! type used to specify a combined extrapolation
typedef tuple<string,double,vector<string>> comb_extr_t;

//! perform a combined sea chiral extrapolation
void combined_sea_chir_extrap(const vector<comb_extr_t> &list);

#endif
