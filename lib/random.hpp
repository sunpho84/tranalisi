#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <boot.hpp>
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
  //! init without a seed
  gen_t() : mt19937_64() {}
  
  //! init with a seed
  gen_t(int seed) : mt19937_64(seed) {}
  
  //! return a real random number in the range [min,max)
  double get_rand_double(double min,double max)
  {return uniform_real_distribution<double>(min,max)(*this);}
  
  //! return an integer random number in the range [min,max)
  double get_rand_int(int min,int max)
  {return uniform_int_distribution<>(min,max-1)(*this);}
};

EXTERN_RANDOM gen_t glb_gen;

void boot_init_t::fill(int seed)
{
  gen_t gen(seed);
  
  for(auto &it : *this) it=gen.get_rand_int(0,njacks);
}

#undef EXTERN_RANDOM
#undef INIT_TO

#endif
