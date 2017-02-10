#ifndef _COMMON_HPP
#define _COMMON_HPP

#include <cstdlib>
#include <tranalisi.hpp>

#ifndef EXTERN_COMMON
 #define EXTERN_COMMON extern
 #define INIT_TO(A)
 #define INIT_ARGS(A)
#else
 #define INIT_TO(A) =A
 #define INIT_ARGS(...) (__VA_ARGS__)
#endif

//totoal number of possible ensemble
const size_t nens_total=15;

//number of input analysis
const size_t ninput_an=8;

//! number of beta
const size_t nbeta=3;

//! number of boostrap
const size_t nboots=100;

//! charge of up quark
const double eu=2.0/3;
//! charge of down quark
const double ed=-1.0/3;
//! charge of strange quark
const double es=ed;
//! charge of charm
const double ec=eu;

//! hold the jacknife index for the given bootstrap
EXTERN_COMMON boot_init_t jack_index[ninput_an][nens_total];

//! symbols for plots
const vector<size_t> symbol={grace::SQUARE,grace::CIRCLE,grace::DIAMOND};
//! colors for plots
const vector<size_t> color={grace::GREEN4,grace::RED,grace::BLUE};

//! class to hold results from mass determination
class lat_par_t
{
public:
  dboot_t ml,ms,mc,r0,f0,B0;
  dbvec_t ainv,Z;
 
  lat_par_t() : ainv(nbeta),Z(nbeta) {}
};

//! results from mass determinations
EXTERN_COMMON vector<lat_par_t> lat_par INIT_ARGS(ninput_an);

//! initialize the common part of the IB analsys
void init_common_IB(string ens_pars);

//! read an ascii bootstrap
dboot_t read_boot(const raw_file_t &file);

//! perform the analysis according to eq.28
ave_err_t eq_28_analysis(const dbvec_t &v);

#undef INIT_TO
#undef EXTERN_COMMON

#endif
