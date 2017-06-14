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

//! get the size to init a valarray
template <class T> size_t init_nel(const valarray<T> &obj)
{return obj.size();}

//////////////////////////////////////////////// operations //////////////////////////////////////////////

#define DEFINE_BIN_OPERATOR(OP)						\
  /* operation between vectors */					\
  template <class TS>							\
  vector<TS> operator OP(const vector<TS> &first,const vector<TS> &second) \
  {									\
    vector<TS> out(first.size());					\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second[it]; \
    return out;								\
  }									\
  /* self operation between vectors */					\
  template <class TS>							\
  vector<TS> operator OP##=(const vector<TS> &first,const vector<TS> &second) \
  {									\
    for(size_t it=0;it<first.size();it++) first[it] OP##= second[it];	\
    return first;							\
  }									\
  /* operation between arrays */					\
  template <class TS,int N>						\
  array<TS,N> operator OP(const array<TS,N> &first,const array<TS,N> &second) \
  {									\
    array<TS,N> out(first.size());					\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second[it]; \
    return out;								\
  }									\
  /* self operation between arrays */					\
  template <class TS,int N>						\
  array<TS,N> operator OP##=(const array<TS,N> &first,const array<TS,N> &second) \
  {									\
    for(size_t it=0;it<first.size();it++) first[it] OP##= second[it];	\
    return first;							\
  }									\
  /* operation between vector and scalar */				\
    template <class TV,class TS,class=enable_if_t<is_vector<TV>::value and is_arithmetic<TS>::value and !(is_same<base_type_t<TV>,TS>::value)>> \
  TV operator OP(const TV &first,const TS &second)			\
  {									\
    TV out(init_nel(first));						\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second;	\
    return out;								\
  }									\
  /* opposite */							\
    template <class TV,class TS,class=enable_if_t<is_vector<TV>::value and is_arithmetic<TS>::value and !(is_same<base_type_t<TV>,TS>::value)>> \
    TV operator OP(const TS &first,const TV &second)			\
  {									\
    TV out(init_nel(second));						\
    for(size_t it=0;it<second.size();it++) out[it]=first OP second[it];	\
    return out;								\
  }									\
  /* self-version of a given operator */				\
    template <class TV,class TS,class=enable_if_t<is_vector<TV>::value and is_arithmetic<TS>::value and !(is_same<base_type_t<TV>,TS>::value)>> \
    auto operator OP##=(TV &first,const TS &second) -> decltype(first OP second) \
    {return first=first OP second;}

DEFINE_BIN_OPERATOR(+)
DEFINE_BIN_OPERATOR(-)
DEFINE_BIN_OPERATOR(*)
DEFINE_BIN_OPERATOR(/)

#define DEFINE_BIN_BOOL_OPERATOR(OP)					\
  template <class TV,class=enable_if_t<is_vector<TV>::value>>		\
  bool operator OP(const TV &first,const TV &second)			\
  {									\
    bool out=true;							\
    for(size_t it=0;it<first.size();it++) out&=(first[it] OP second[it]); \
    return out;								\
  }
DEFINE_BIN_BOOL_OPERATOR(!=)
DEFINE_BIN_BOOL_OPERATOR(==)
DEFINE_BIN_BOOL_OPERATOR(>)
DEFINE_BIN_BOOL_OPERATOR(<)

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
DEFINE_FUNCTION(abs)
DEFINE_FUNCTION(asinh)
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

//! specify how to print a vector
template <class TV,class=enable_if_t<is_vector<TV>::value and !is_same<TV,string>::value>>
ostream& operator<<(ostream &out,const TV &v)
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

//! take the backward derivative
template <class T> T backward_derivative(const T &v)
{
  T out(v.size());
  out[0]=0.0;
  for(size_t it=1;it<out.size();it++) out[it]=v[it]-v[it-1];
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

//! correlation
template <class T> double corr(const T &x,const T &y)
{return cov(x,y)/sqrt(var(x)*var(y));}

#endif
