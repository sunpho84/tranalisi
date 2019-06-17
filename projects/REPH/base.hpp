#ifndef _BASE_HPP
#define _BASE_HPP

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <raw_file.hpp>

const char VA_tag[2][2]={"V","A"};

//! Holds all info for a given physical quark
using quark_t=tuple<double,vector<size_t>>;

//! List of physical quarks
//!
//! Each entry in the map is connected to a physical quark (Up, Down, etc)
//! Each value is a tuple: the first component is the charge, the second is a vector, listing all masses in the data
map<string,quark_t> quarks;

//! Holds all info for a given meson
using meson_t=pair<string,string>;

//! List of all known mesons
//!
//! Each entry in the map is connected to a physical meson (Pi+, D, etc)
//! Each value is a pair of strings, pointing to the physical quark
map<string,meson_t> mesons;

#endif
