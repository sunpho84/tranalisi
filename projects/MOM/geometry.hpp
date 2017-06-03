#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <array>
#include <map>
#include <vector>

#ifndef EXTERN_GEOMETRY
 #define EXTERN_GEOMETRY extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

using namespace std;

//! number of dimensions
const int NDIM=4;

//! holds a coordinate
using coords_t=array<int,NDIM>;

//! holds a momentum - to be extended
class mom_t : public coords_t
{
public:
};

//! lattice sizes
EXTERN_GEOMETRY coords_t L;

//! list of momenta
EXTERN_GEOMETRY vector<mom_t> moms;

//! reads the list of momenta from file
void get_list_of_moms(const string &path);

//! list of class of equivalence of momenta
EXTERN_GEOMETRY map<mom_t,vector<size_t>> equiv_moms;

//! fills the list of class of equivalence
void get_class_of_equiv_moms();

#undef EXTERN_GEOMETRY
#undef INIT_TO

#endif
