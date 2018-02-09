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
//! hadron path
EXTERN_READ string prop_hadr_path;
//! lepton path
EXTERN_READ string prop_lep_path;
//! number of flavour
EXTERN_READ ev::Nf_t Nf;
//! order of alpha to use
EXTERN_READ size_t evo_ord;
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

//! open all qprops
vector<raw_file_t> setup_read_all_qprops_mom(const vector<size_t> &conf_list);

//! open all lprops
vector<raw_file_t> setup_read_all_lprops_mom(const vector<size_t> &conf_list);

//! real all propagators for a given i_in_clust_ihit and imom, for quarks
vector<m_r_mom_conf_qprops_t> read_all_qprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);

//! real all propagators for a given i_in_clust_ihit and imom, for leptons
vector<mom_conf_lprops_t> read_all_lprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);

#undef EXTERN_READ
#undef INIT_TO

#endif
