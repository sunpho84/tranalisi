#ifndef _OPER_HPP
#define _OPER_HPP

#include <algorithm>
#include <cmath>
#include <functions.hpp>
#include <functional>
#include <iostream>
#include <jack.hpp>
#include <macros.hpp>
#include <solve.hpp>
#include <tools.hpp>
#include <type_traits>
#include <valarray>
#include <vector>

using namespace std;
using namespace placeholders;

//! get the size to init an object, avoiding size() for boot and jack
template <class T> size_t init_nel(const valarray<T> &obj)
{return obj.size();}

//////////////////////////////////////////////// operations //////////////////////////////////////////////

//! operation between two vectors
// template <class T,class=enable_if_t<is_vector<T>::value>> T OP_NAME(const T &first,const T &second)
// {
//   check_match_size(first,second);
//
//   T out(init_nel(first));
//   for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second[it];
//   return out;
// }


#define DEFINE_BIN_OPERATOR(OP)						\
  /* operation between vector and scalar */				\
  template <class TV,class TS,class=enable_if_t<is_vector<TV>::value and is_arithmetic<TS>::value and !(is_same<typename TV::base_type,TS>::value)>> \
    TV operator OP(const TV &first,const TS &second)			\
  {									\
    TV out(init_nel(first));						\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second;	\
    return out;								\
  }									\
  /* opposite */							\
    template <class TV,class TS,class=enable_if_t<is_vector<TV>::value and is_arithmetic<TS>::value and !(is_same<typename TV::base_type,TS>::value)>> \
    TV operator OP(const TS &first,const TV &second)			\
  {									\
    TV out(init_nel(second));						\
    for(size_t it=0;it<second.size();it++) out[it]=first OP second[it];	\
    return out;								\
  }									\
  /* self-version of a given operator */				\
    template <class TV,class TS,class=enable_if_t<is_vector<TV>::value and is_arithmetic<TS>::value and !(is_same<typename TV::base_type,TS>::value)>> \
    auto operator OP##=(TV &first,const TS &second) -> decltype(first OP second) \
    {return first=first OP second;}

DEFINE_BIN_OPERATOR(+)
DEFINE_BIN_OPERATOR(-)
DEFINE_BIN_OPERATOR(*)
DEFINE_BIN_OPERATOR(/)

//! function of a vector
#define DEFINE_NAMED_FUNCTION(OP_NAME,OP)				\
  template <class T,class ...Args,class=enable_if<is_vector<T>::value>> T OP_NAME(const T &first,Args... args) \
  {									\
    T out(init_nel(first));						\
    for(size_t it=0;it<first.size();it++) out[it]=OP(first[it],args...); \
    return out;								\
  }
#define DEFINE_FUNCTION(OP_NAME) DEFINE_NAMED_FUNCTION(OP_NAME,OP_NAME)

// DEFINE_NAMED_FUNCTION(operator-,uminus)
// DEFINE_FUNCTION(cos)
// DEFINE_FUNCTION(cosh)
// DEFINE_FUNCTION(exp)
// DEFINE_FUNCTION(log)
// DEFINE_FUNCTION(pow)
// DEFINE_FUNCTION(sin)
// DEFINE_FUNCTION(sinh)
// DEFINE_FUNCTION(sqr)
// DEFINE_FUNCTION(sqrt)
// DEFINE_FUNCTION(tan)
// DEFINE_FUNCTION(tanh)

//! check that the object is printable
template <class T> enable_if_t<!has_method_is_printable<T>::value,bool> is_printable(const T &o)
{return true;}

//! check that the object is printable
template <class T> enable_if_t<has_method_is_printable<T>::value,bool> is_printable(const T &o)
{return o.is_printable();}

//! specify hot to print a vector
template <class T> ostream& operator<<(ostream &out,const vector<T> &v)
{
  for(size_t it=0;it<v.size();it++) if(is_printable(v[it])) out<<it<<" "<<v[it]<<endl;
  return out;
}

//! transpose a matrix
template <class T> vector<vector<T>> transpose(const vector<vector<T>> &in)
{
  auto start=take_time();
  
  vector<vector<T>> out(in[0].size(),vector<T>(in.size()));
  for(size_t ir=0;ir<in.size();ir++)
    for(size_t ic=0;ic<in[0].size();ic++)
      out[ic][ir]=in[ir][ic];
  
  cout<<elapsed_time(start)<<" to transpose"<<endl;
  
  return out;
}

//! take the forward derivative
template <class T> T forward_derivative(const T &v)
{
  T out(v.size()-1);
  for(size_t it=0;it<out.size();it++) out[it]=v[it+1]-v[it];
  return out;
}

//! extract a subset of a vector
template <class T> T subset(const T &v,size_t beg,size_t end)
{
  check_ordered({beg,end,v.size()});
  return T(&v[beg],&v[end]);
}

//! compute the covariance
template <class T> double cov(const T &x,const T &y)
{return (T(x*y).ave()-x.ave()*y.ave())*(njacks-1);}

//! variance
template <class T> double var(const T &x)
{return cov(x,x);}


#endif
