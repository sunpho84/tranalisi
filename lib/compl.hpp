#ifndef _COMPL_HPP
#define _COMPL_HPP

#include <complex>

using namespace std;

using dcompl_t=complex<double>;

//! get the real or imag part
template<typename T>
T& get_re_or_im(complex<T> &z,size_t ri)
{return ((T*)(&z))[ri];}

//! get the real or imag part
template<typename T>
const T get_re_or_im(const complex<T> &z,size_t ri)
{return ((T*)(&z))[ri];}

//! get the real part
template<typename T>
T& get_real(complex<T> &z)
{return get_re_or_im(z,1);}

//! get the imag part
template<typename T>
T& get_imag(complex<T> &z)
{return get_re_or_im(z,1);}

#endif
