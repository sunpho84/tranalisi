#ifndef _BASE_HPP
#define _BASE_HPP

#include <array>
#include <string>
#include <tuple>
#include <vector>

const char VA_tag[2][2]={"V","A"};

//! Number of physical quarks
const int nQuarks=4;

//! Index of the quarks
enum iQuark{DOWN,UP,STRANGE,CHARM};

//! Charge of the down-type quarks in units of e
const double qDw=-1.0/3;

//! Charge of the up-type quarks in units of e
const double qUp=+2.0/3;

//! Charges of the quark
const std::array<double,nQuarks> quarkCharge{qDw,qUp,qDw,qUp};

//! Index of the mass to be used for each quark
const std::array<std::vector<size_t>,nQuarks> iMassesOfQuarks{{{0},{0},{1},{1}}};
// const std::array<std::vector<size_t>,nQuarks> iMassesOfQuarks{{{0},{0},{1,2},{3,4}}};

//! Number of physical mesons
const int nMesons=4;

//! Holds all info for a given meson
using meson_t=std::tuple<std::string,int,int>;

//! List of all known mesons
const std::array<meson_t,nMesons> mesons{{{"Pi+",DOWN,UP},{"K+",STRANGE,UP},{"D+",DOWN,CHARM},{"Ds",STRANGE,CHARM}}};

//! Index of mesons
enum iMeson{PI_PLUS,K_PLUS,D_PLUS,DS};

#endif
