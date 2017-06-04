#ifndef _COMPL_HPP
#define _COMPL_HPP

#include <complex>

using namespace std;

using dcompl_t=complex<double>;

//! get the real part
template<typename T>
T& get_real(complex<T> &z)
{
  struct Z
  {
    T re,im;
  };
  static_assert(sizeof(Z)==sizeof(std::complex<T>),"!!");
  static_assert(alignof(Z)==alignof(std::complex<T>),"!!");
  return static_cast<Z&>(z).re;
}

//! get the imag part
template<typename T>
T& get_imag(complex<T> &z)
{
  struct Z
  {
    T re,im;
  };
  static_assert(sizeof(Z)==sizeof(std::complex<T>),"!!");
  static_assert(alignof(Z)==alignof(std::complex<T>),"!!");
  return static_cast<Z&>(z).im;
}

#endif
