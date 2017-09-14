#ifndef _FPIEM_FSE_HPP
#define _FPIEM_FSE_HPP

#include <jack.hpp>

using namespace std;

double FSE_V(double mpi,double L,double fpi,double B);

template <class T>
T FSE_V(const T &mpi,double L,const T &fpi,double B)
{
  T out;
  for(size_t iel=0;iel<mpi.size();iel++) out[iel]=FSE_V(mpi[iel],L,fpi[iel],B);
  return out;
};

#endif
