#ifndef _READ_HPP
#define _READ_HPP

#include <corrections.hpp>
#include <evolutions.hpp>

#include <string>

using namespace std;

#ifndef EXTERN_READ
#define EXTERN_READ extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

//! number of momenta between each print
EXTERN_READ size_t print_each_mom;
//! gauge action
EXTERN_READ gaz::type_t act;
//! beta
EXTERN_READ double beta;
//! plaquette
EXTERN_READ double plaq;
//! maximal quark mass
EXTERN_READ double am_max;
//! minimal quark mass
EXTERN_READ double am_min;
//! index of sea mass
EXTERN_READ size_t im_sea;
//! number of hits
EXTERN_READ size_t nhits;
//! number of hits to be used
EXTERN_READ size_t nhits_to_use;
//! minimal time for fit
EXTERN_READ size_t tmin;
//! maximal time for fit
EXTERN_READ size_t tmax;
//! inverse lattice spacing
EXTERN_READ double ainv;
//! conf range
EXTERN_READ range_t conf_range;
//! number of flavour
EXTERN_READ ev::Nf_t Nf;
//! order of alpha to use
EXTERN_READ size_t ord;
//! coupling
EXTERN_READ double g2;
//! boosted coupling
EXTERN_READ double g2tilde;
//! list of existing confs
EXTERN_READ vector<size_t> conf_list;
//! cluster size
EXTERN_READ size_t clust_size;

//! read input file and setup basic things
void read_input(const string &input_path);

//! prepare the list of confs
void prepare_list_of_confs();

//! prepare a list of reading task, to be executed in parallel
vector<task_list_t> prepare_read_prop_taks(vector<m_r_mom_conf_props_t> &props,const vector<size_t> &conf_list,const size_t imom);

#undef EXTERN_READ
#undef INIT_TO

#endif
