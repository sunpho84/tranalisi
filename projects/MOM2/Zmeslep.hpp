#ifndef _ZMESLEP_HPP
#define _ZMESLEP_HPP

#ifndef EXTERN_ZMESLEP
 #define EXTERN_ZMESLEP extern
 #define INIT_ZMESLEP_TO(...)
#else
 #define INIT_ZMESLEP_TO(...) __VA_ARGS__
#endif

namespace meslep
{
  EXTERN_ZMESLEP double ql   INIT_ZMESLEP_TO(=-1.0);      //!< the program simulates muon *antiparticle*
  EXTERN_ZMESLEP double q_in INIT_ZMESLEP_TO(=+2.0/3.0);  //!< charge of the quark which comes into the vertex
  EXTERN_ZMESLEP double q_ou INIT_ZMESLEP_TO(=-1.0/3.0);  //!< charge of the quark which comes out the vertex
}

#undef INIT_ZMESLEP_TO
#undef EXTERN_ZMESLEP

#endif
