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

//! if we ever needed a global generator...
//EXTERN_RANDOM gen_t glb_gen;

#undef EXTERN_RANDOM
#undef INIT_TO

#endif
