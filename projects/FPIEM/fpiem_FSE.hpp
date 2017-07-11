#ifndef _FPIEM_FSE_HPP
#define _FPIEM_FSE_HPP

#include <jack.hpp>

using namespace std;

double FSE_V(double mpi,double L,double fpi,double B);
djack_t FSE_V(const djack_t &mpi,double L,const djack_t &fpi,double B)
{
  djack_t out;
  for(size_t ijack=0;ijack<=njacks;ijack++) out[ijack]=FSE_V(mpi[ijack],L,fpi[ijack],B);
  return out;
};

#endif
