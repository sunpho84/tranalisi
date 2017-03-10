#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <iostream>
#include <random>
#include <vector>

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#include <ave_err.hpp>

using namespace std;

//! random generator
class gen_t: public mt19937_64
{
public:
  //! init with a seed
  gen_t(int seed) : mt19937_64(seed) {}
  
  //! return a real random number in the range [min,max)
  double get_double(double min,double max)
  {return uniform_real_distribution<double>(min,max)(*this);}
  
  //! return an integer random number in the range [min,max)
  double get_int(int min,int max)
  {return uniform_int_distribution<>(min,max-1)(*this);}
  
  //! return an integer random number in the range [min,max)
  double get_gauss(double ave,double sig)
  {return normal_distribution<>(ave,sig)(*this);}
  
private:
  //! init without a seed
  gen_t() : mt19937_64() {}
};

///////////////////////////////////////////////////////////// gauss_filler_t /////////////////////////////////////////////////////

//! allows to fill from gauss
class gauss_filler_t : pair<ave_err_t,int>
{
public:
  //! fill from ave_err and seed
  gauss_filler_t(const ave_err_t &ext_ae,int ext_seed)
  {
    ae=ext_ae;
    seed=ext_seed;
  }
  
  //! fill from ave, err and seed
  gauss_filler_t(double ave,double err,int ext_seed) : gauss_filler_t(ave_err_t(ave,err),ext_seed) {}
  
  //! rebind ave_err
  ave_err_t &ae=first;
  
  //! rebind seed
  int &seed=second;
};

//! if we ever needed a global generator...
//EXTERN_RANDOM gen_t glb_gen;

#undef EXTERN_RANDOM
#undef INIT_TO

#endif
