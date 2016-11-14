#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_JACK
#include <jack.hpp>
#include <meas_vec.hpp>
#ifdef USE_OMP
 #include <omp.h>
#endif
#include <oper.hpp>
#include <tools.hpp>

void set_njacks(int ext_njacks)
{
  if(njacks==UNDEF_NJACKS) njacks=ext_njacks;
  else CRASH("Unbale to set njacks twice");
}

void check_njacks_init()
{if(njacks==UNDEF_NJACKS) CRASH("Set njacks before");}
