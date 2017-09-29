#ifndef _MACROS_HPP
#define _MACROS_HPP

#include <array>
#include <cmath>

using namespace std;

#define CRASH(...) internal_crash(__LINE__,__FILE__,__PRETTY_FUNCTION__,__VA_ARGS__)
#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)
#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)

const double gammaeul=0.5772156649015328;
const double alpha_em=1/137.04;
const double e2=4*M_PI*alpha_em;
const double mu_MS=2.0;

enum{RE,IM};

enum{SILENT,VERBOSE};

//! check if has method
#define DEFINE_HAS_METHOD(NAME)						\
  template<class C, typename T = void>					\
  struct has_method_##NAME						\
  {									\
    typedef char (&yes)[2];						\
    template<unsigned long> struct exists;				\
    template<typename V> static yes check(exists<sizeof(static_cast<T>(&V::NAME))>*); \
    template<typename> static char Check (...);				\
    static const bool value=(sizeof(Check<C>(0))==sizeof(yes));		\
  };									\
  template<class C>							\
  struct has_method_##NAME<C,void>					\
  {									\
    typedef char (&yes)[2];						\
    template<unsigned long> struct exists;				\
    template<typename V> static yes check(exists<sizeof(&V::NAME)>*);	\
    template<typename> static char check(...);				\
    static const bool value=(sizeof(check<C>(0))==sizeof(yes));		\
  }

#endif
