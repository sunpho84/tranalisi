#ifndef _JACK_HPP
#define _JACK_HPP

#ifndef EXTERN_JACK
 #define EXTERN_JACK extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

//! number of jackknife
#define UNDEF_NJACKS 0
EXTERN_JACK int njacks INIT_TO(UNDEF_NJACKS);

//! set the number of jackknives
void set_njacks(int ext_njacks);

//! crash if number of jackknives is not initialized
void check_njacks_init();

#undef EXTERN_JACK
#undef INIT_TO

#endif
