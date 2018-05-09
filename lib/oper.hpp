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
#include <stopwatch.hpp>
#include <tools.hpp>
#include <utility>
#include <valarray>
#include <vector>

#include <traits.hpp>

using namespace std;
using namespace placeholders;

//! get the size to init a valarray
template <class T> size_t
init_nel(const valarray<T> &obj)
{return obj.size();}

//! get the size to init a vector
template <class T> size_t
init_nel(const vector<T> &obj)
{return obj.size();}


////////////////////////////////////////////// a scalar wrapper //////////////////////////////////////////

template <class T,size_t deg>
class SC_t
{
 public:
  const SC_t<T,deg-1> val;
  SC_t(const T &val) : val(val) {}
};

//! 1-level wrapper
template <class T>
class SC_t<T,1>
{
 public:
  const T &val;
  SC_t(const T &val) : val(val) {}
};

//! return a scalar wrapper around T
template <class T,const size_t deg=1>
auto SC(const T &val) -> SC_t<T,deg>
{return SC_t<T,deg>(val);}

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
    {return first=first OP second;}					\
    /* operation between a generic vector and a scalar of type different from base_type */ \
    template<class TV,class TT=vector_traits<TV>,class TS,size_t deg,class=enable_if_t<TT::is_vector and !is_same<typename TT::base_type,TS>::value>> \
    TV operator OP(const TV &first,const SC_t<TS,deg> &second)		\
    {									\
      TV out(init_nel(first));						\
      for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second.val; \
      return out;							\
    }									\
    /* operation between a generic vector and a scalar of type different from base_type */ \
    template<class TV,class TT=vector_traits<TV>,class TS,size_t deg,class=enable_if_t<TT::is_vector and !is_same<typename TT::base_type,TS>::value>> \
    TV operator OP(const SC_t<TS,deg> &first,const TV &second)		\
    {									\
      TV out(init_nel(second));						\
      for(size_t it=0;it<first.size();it++) out[it]=first.val OP second[it]; \
      return out;							\
    }

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
  const size_t n=v.size();
  T out(n);
  for(size_t it=0;it<n;it++) out[it]=v[(it+1)%n]-v[it];
  return out;
}

//! take the backward derivative
template <class T> T backward_derivative(const T &v)
{
  const size_t n=v.size();
  T out(n);
  for(size_t it=0;it<n;it++) out[it]=v[it]-v[(it-1+n)%n];
  return out;
}

//! take the symmetric derivative
template <class T> T symmetric_derivative(const T &v)
{
  const size_t n=v.size();
  T out(n);
  for(size_t it=0;it<n;it++) out[it]=(v[(it+1)%n]-v[(it-1+n)%n])/2.0;
  return out;
}

//! extract a subset of a vector
template <class T> T subset(const T &v,size_t beg,size_t end)
{
  check_ordered({beg,end,v.size()});
  return T(&v[beg],&v[end]);
}

//! compute the covariance error
template <class T>
T distr_cov(const T &x,const T &y)
{
  T xy=x*y;
  
  size_t nel=x.size();
  T a;
  double xa=x.ave();
  double ya=y.ave();
  double xya=xy.ave();
  
  for(size_t iel=0;iel<nel;iel++)
    {
      double xi=x[iel];
      double yi=y[iel];
      double xyi=xy[iel];
      
      double xb=(xa*nel-xi)/(nel-1);
      double yb=(ya*nel-yi)/(nel-1);
      double xyb=(xya*nel-xyi)/(nel-1);
      
      a[iel]=(xyb-xb*yb)*(njacks-1);
    }
  
  // //renormalize
  // double t=a.ave();
  // for(size_t iel=0;iel<nel;iel++)
  //   a[iel]=t+(a[iel]-t)*(nel-1)/(njacks-1);
  
  return a;
}

//! compute the covariance
template <class T>
double cov(const T &x,const T &y)
{return (T(x*y).ave()-x.ave()*y.ave())*(njacks-1);}

//! variance
template <class T>
double var(const T &x)
{return cov(x,x);}

//! correlation
template <class T>
double corr(const T &x,const T &y)
{return cov(x,y)/sqrt(var(x)*var(y));}

//! filter a valarray
template <class T,class=enable_if_t<has_method_size<T>::value>>
T vec_filter(const T &v,const gslice &slice)
{return (T)(v[slice]);}

//! concatenate an arbitrary number of vectors
template <class T>
vector<T> concat_internal(const vector<const vector<T>*> &in)
{
  vector<T> out;
  
  //compute eize and reserve
  size_t size=0;
  for(auto &v : in) size+=v->size();
  out.reserve(size);
  
  //push back all
  for(auto &v : in)
    for(auto &vi : *v)
      out.push_back(vi);
  
  return out;
}

//! concatenate one vector to another
template <class T,typename ...Rest>
vector<T> concat_internal(vector<const vector<T>*> &vlist,const vector<T> &v,Rest&&... rest)
{
  vlist.push_back(&v);
  return concat_internal(vlist,forward<Rest>(rest)...);
}

//! concatenate one vector and a const
template <class T,typename ...Rest>
vector<T> concat_internal(vector<const vector<T>*> &vlist,const T &v,Rest&&... rest)
{
  vector<T> tv={v};
  vlist.push_back(&tv);
  return concat_internal(vlist,forward<Rest>(rest)...);
}

//! concatenate one vector and an initializer list
template <class T,typename ...Rest>
vector<T> concat_internal(vector<const vector<T>*> &vlist,const initializer_list<T> &v,Rest&&... rest)
{
  vector<T> tv={v};
  vlist.push_back(&tv);
  return concat_internal(vlist,forward<Rest>(rest)...);
}
//! concatenate at least two vectors
template <class T,typename... Rest>
vector<T> concat(const vector<T> &v,Rest&&... rest)
{
  vector<const vector<T>*> vlist({&v});
  return concat_internal(vlist,forward<Rest>(rest)...);
}

//! concatenate at least an element and a vectors
template <class T,typename... Rest>
vector<T> concat(const T&c,const vector<T> &v,Rest&&... rest)
{
  const vector<T> vc{c};
  vector<const vector<T>*> vlist({&vc,&v});
  return concat_internal(vlist,forward<Rest>(rest)...);
}

#endif
