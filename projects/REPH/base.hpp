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
using meson_t=tuple<string,string,string>;

//! List of all known mesons
//!
//! Each entry is connected to a physical meson (Pi+, D, etc)
//! Each value is a pair of strings, pointing to the physical quark
vector<meson_t> mesonList;

#endif
