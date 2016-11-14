#ifndef _OPER_HPP
#define _OPER_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <macros.hpp>
#include <tools.hpp>
#include <type_traits>
#include <vector>

using namespace std;
using namespace placeholders;

//! check agreement of sizes of two vectors
template <class T1,class T2> void check_match_size(const vector<T1> &first,const vector<T2> &second)
{if(first.size()!=second.size()) CRASH("Vectors do not agree in size, %d vs %d",first.size(),second.size());}

//! get the size to init an object, avoiding size() for boot and jack
template <class T> const size_t init_nel(const vector<T> &obj)
{return obj.size();}

//////////////////////////////////////////////// operations //////////////////////////////////////////////

//! operation between two vectors
#define DEFINE_BIN_OPERATOR(OP_NAME,SELF_OP_NAME,OP)			\
  template <class T,class=enable_if_t<is_vector<T>::value>> T OP_NAME(const T &first,const T &second)	\
  {									\
    check_match_size(first,second);					\
    									\
    T out(init_nel(first));						\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second[it]; \
    return out;								\
  }									\
  /* operation between vector and scalar */				\
  template <class TV,class TS,class=enable_if_t<is_vector<TV>::value&&(is_base_of<vector<TS>,TV>::value||is_arithmetic<TS>::value)>> \
    TV OP_NAME(const TV &first,const TS &second)			\
  {									\
    TV out(init_nel(first));						\
    for(size_t it=0;it<first.size();it++) out[it]=first[it] OP second;	\
    return out;								\
  }									\
  /* opposite */							\
  template <class TV,class TS,class=enable_if_t<is_vector<TV>::value&&(is_base_of<vector<TS>,TV>::value||is_arithmetic<TS>::value)>> \
    TV OP_NAME(const TS &first,const TV &second)			\
  {									\
    TV out(init_nel(first));						\
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
  template <class T,class ...Args,class=enable_if<is_vector<T>::value>> T OP_NAME(const T &first,Args... args) \
  {									\
    T out(init_nel(first));						\
    for(size_t it=0;it<first.size();it++) out[it]=OP(first[it],args...); \
    return out;								\
  }
#define DEFINE_FUNCTION(OP_NAME) DEFINE_NAMED_FUNCTION(OP_NAME,OP_NAME)

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

//! specify hot to print a vector
template <class T> ostream& operator<<(ostream &out,const vector<T> &v)
{
  for(size_t it=0;it<v.size();it++) out<<it<<" "<<v[it]<<endl;
  return out;
}

//! transpose a matrix
template <class T> vector<vector<T>> transpose(const vector<vector<T>> &in)
{
  START_TIME();
  
  vector<vector<T>> out(in[0].size(),vector<T>(in.size()));
  for(size_t ir=0;ir<in.size();ir++)
    for(size_t ic=0;ic<in[0].size();ic++)
      out[ic][ir]=in[ir][ic];
  
  cout<<ELAPSED_TIME()<<" to transpose"<<endl;
  
  return out;
}

//! hold a triplet to set a filter
class filter_t
{
public:
  size_t each;
  size_t offset;
  size_t how_many;
  
  //! init from three pars
  filter_t(size_t each,size_t offset=0,size_t how_many=1) : each(each),offset(offset),how_many(how_many) {}
  
  //! default constructor
  filter_t()=default;
  
  //! filter vector
  template <class T> T operator()(const T &v)
  {
    START_TIME();
    
    T out;
    for(size_t it=offset;it<v.size();it+=each)
      if(it+how_many<v.size()) //check that we store a whole bunch
	for(size_t sh=0;sh<how_many;sh++)
	  out.push_back(v[it+sh]);
    
    cout<<ELAPSED_TIME()<<" to filter data"<<endl;
    
    return out;
  }
};

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
  if(beg>=v.size()) CRASH("Asked to extract from %zu, beyond the end %zu",beg,v.size());
  if(end>=v.size()) CRASH("Asked to extract up to %zu, beyond the end %zu",end,v.size());
  
  return T(&v[beg],&v[end]);
}

#endif
