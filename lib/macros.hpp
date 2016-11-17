#ifndef _MACROS_HPP
#define _MACROS_HPP

#include <type_traits>

using namespace std;

#if __cplusplus<201402L
 template <bool B,class T=void> using enable_if_t=typename enable_if<B,T>::type;
#endif

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)
#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)

enum{RE,IM};

//! check if has method
#define DEFINE_HAS_METHOD(METHOD)					\
  template <typename T> class NAME2(has_method,METHOD)			\
  {									\
  private:								\
    typedef char yes;							\
    typedef yes no[2];							\
    									\
    template<typename C> static auto test(void*)->decltype(size_t{std::declval<C const>().METHOD()},yes{}); \
    template<typename> static no& test(...);				\
    									\
  public:								\
    static bool const value=sizeof(test<T>(0))==sizeof(yes);		\
  }

#endif
