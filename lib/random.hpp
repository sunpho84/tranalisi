#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#undef EXTERN_BOOT
#undef INIT_TO

#endif
