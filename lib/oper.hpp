#ifndef _OPER_HPP
#define _OPER_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <tools.hpp>
#include <type_traits>
#include <vector>

using namespace std;
using namespace placeholders;

#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)
#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)

//! check agreement of sizes of two vectors
template <class T1,class T2> void check_match_size(const vector<T1> &first,const vector<T2> &second)
{if(first.size()!=second.size()) CRASH("Vectors do not agree in size, %d vs %d",first.size(),second.size());}

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

DEFINE_HAS_METHOD(size);

//////////////////////////////////////////////// sum /////////////////////////////////////////////////

//! operation between two vectors
#define DEFINE_BIN_OPERATOR(OP_NAME,SELF_OP_NAME,OP)			\
  template <class T,class = typename enable_if<has_method_size<T>::value>::type> T OP_NAME(const T &first,const T &second)	\
  {									\
    check_match_size(first,second);					\
    									\
    T out(first.size());						\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second[it]; \
    return out;								\
  }									\
  /* operation between vector and scalar */				\
  template <class TV,class TS,class = typename enable_if<!is_base_of<vector<TV>,TS>::value>::type> vector<TV> OP_NAME(const vector<TV> &first,const TS &second) \
  {									\
    vector<TV> out(first.size());					\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second;	\
    return out;								\
  }									\
  /* opposite */							\
  template <class TV,class TS,class = typename enable_if<!is_base_of<vector<TV>,TS>::value>::type> vector<TV> OP_NAME(const TS &first,const vector<TV> &second) \
  {									\
    vector<TV> out(first.size());					\
    for(size_t it=0;it<first.size();it++) out[it]=first OP second[it];	\
    return out;								\
  }									\
  /* self-version of a given operator */				\
  template <class T1,class T2> auto SELF_OP_NAME(T1 &first,const T2 &second) -> decltype(first OP second) \
  {return first=first OP second;}

DEFINE_BIN_OPERATOR(operator+,operator+=,+)
DEFINE_BIN_OPERATOR(operator-,operator-=,-)
DEFINE_BIN_OPERATOR(operator*,operator*=,*)
DEFINE_BIN_OPERATOR(operator/,operator/=,/)

//! function of a vector
#define DEFINE_NAMED_FUNCTION(OP_NAME,OP)				\
  template <class T,class ...Args,class = typename enable_if<has_method_size<T>::value>::type> T OP_NAME(const T &first,Args... args) \
  {									\
    T out(first.size());						\
    for(size_t it=0;it<first.size();it++) out[it]=OP(first[it],args...); \
    return out;								\
  }
#define DEFINE_FUNCTION(OP_NAME) DEFINE_NAMED_FUNCTION(OP_NAME,OP_NAME)

//helper for unary minus
template <class T> T uminus(const T &in)
{return -in;}

DEFINE_NAMED_FUNCTION(operator-,uminus)
DEFINE_FUNCTION(cos)
DEFINE_FUNCTION(cosh)
DEFINE_FUNCTION(exp)
DEFINE_FUNCTION(log)
DEFINE_FUNCTION(pow)
DEFINE_FUNCTION(sin)
DEFINE_FUNCTION(sinh)
DEFINE_FUNCTION(sqr)
DEFINE_FUNCTION(sqrt)
DEFINE_FUNCTION(tan)
DEFINE_FUNCTION(tanh)

#endif
