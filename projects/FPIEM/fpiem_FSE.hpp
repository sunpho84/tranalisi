#ifndef _FPIEM_FSE_HPP
#define _FPIEM_FSE_HPP

#include <array>

using namespace std;

double FSE(double mpi,double L,double fpi,const array<double,3> &B);
inline double FSE(double mpi,double L,double fpi,const double B)
{return FSE(mpi,L,fpi,{B,B,B});}
double elltheta(double tau,double theta);

#endif
