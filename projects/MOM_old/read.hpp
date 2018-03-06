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

//! read input file and setup basic things
void read_input(const string &input_path);

//! prepare the list of confs
void prepare_list_of_confs();

//! open all qprops
vector<raw_file_t> setup_read_all_qprops_mom(const vector<size_t> &conf_list);

//! open all lprops
vector<raw_file_t> setup_read_all_lprops_mom(const vector<size_t> &conf_list);

//! real all propagators for a given i_in_clust_ihit and imom, for quarks
vector<m_r_mom_conf_qprops_t> read_all_qprops_mom(vector<raw_file_t> &files,const index_t &im_r_ijack_ind,const size_t i_in_clust_ihit,const size_t imom);

//! real all propagators for a given i_in_clust_ihit and imom, for leptons
vector<mom_conf_lprops_t> read_all_lprops_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom);

#undef EXTERN_READ
#undef INIT_TO

#endif
