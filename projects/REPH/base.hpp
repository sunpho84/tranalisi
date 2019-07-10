#ifndef _BASE_HPP
#define _BASE_HPP

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <tranalisi.hpp>
#include <compatibility.hpp>
#include <phys_point.hpp>

const char VA_tag[2][2]={"V","A"};

const char ST_tag[2]={'S','T'};

//! Holds all info for a given physical quark
using quark_t=tuple<double,size_t,vector<size_t>>;

//! List of physical quarks
//!
//! Each entry in the map is connected to a physical quark (Up, Down, etc)
//! Each value is a tuple: the first component is the charge, the second is a vector, listing all masses in the data
map<string,quark_t> quarkList;

//! Holds all info for a given meson: the first entry is the name, the second is the bw quark, the third is the fw one
using meson_t=tuple<string,string,string,double>;

//! List of all known mesons
//!
//! Each entry is connected to a physical meson (Pi+, D, etc)
//! Each value is a pair of strings, pointing to the physical quark
vector<meson_t> mesonList;

//! renormalization constants taken from 1403.4504, table 20 pag. 56, first all M1, then M2
const vector<ave_err_t> Zp_ae({{{0.529,0.007},{0.509,0.004},{0.516,0.002},{0.574,0.004},{0.546,0.002},{0.545,0.002}}});
const vector<ave_err_t> Za_ae({{{0.731,0.008},{0.737,0.005},{0.762,0.004},{0.703,0.002},{0.714,0.002},{0.752,0.002}}});
const vector<ave_err_t> Zv_ae({{{0.587,0.004},{0.603,0.003},{0.655,0.003},{0.608,0.003},{0.614,0.002},{0.657,0.002}}});
const vector<ave_err_t> Zt_ae({{{0.711,0.005},{0.724,0.004},{0.774,0.004},{0.700,0.003},{0.711,0.002},{0.767,0.002}}});

dbvec_t Za(2*nbeta),Zp(2*nbeta),Zv(2*nbeta),Zt(2*nbeta);

#endif
