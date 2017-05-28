#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <array>

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

//! holds a momentum
class mom_t : public coords_t
{
public:
};

//! lattice sizes
EXTERN_GEOMETRY coords_t L;

//! return the volume of a given box
inline size_t vol_of_lx(coords_t size)
{
  size_t vol=1;
  for(size_t mu=0;mu<NDIM;mu++) vol*=size[mu];
  return vol;
}

//! returns the coordinate inside a box
inline coords_t coord_of_lx(size_t ilx,const coords_t &s)
{
  coords_t x;
  for(int mu=NDIM-1;mu>=0;mu--)
    {
      x[mu]=ilx%s[mu];
      ilx/=s[mu];
    }
  
  return x;
}

//! sum two coords_t
inline coords_t operator+(const coords_t &a1,const coords_t &a2)
{
  coords_t s;
  for(int mu=0;mu<NDIM;mu++) s[mu]=(a1[mu]+a2[mu])%L[mu];
  return s;
}

//! summassign two coords_t
inline coords_t operator+=(coords_t &a1,const coords_t &a2)
{return a1=a1+a2;}

//! get mirrorized coords according to a bit decomposition of imir
inline coords_t get_mirrorized_site_coords(const coords_t &c,int imir)
{
  coords_t cmir;
  for(int mu=0;mu<NDIM;mu++)
    cmir[mu]=(L[mu]+(1-2*((imir>>mu)&1))*c[mu])%L[mu];
  return cmir;
}

#undef EXTERN_GEOMETRY
#undef INIT_TO

#endif
