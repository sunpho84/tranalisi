#ifndef _OPER_HPP
#define _OPER_HPP

#include <algorithm>
#include <functional>
#include <iostream>
#include <tools.hpp>
#include <vector>

using namespace std;
using namespace placeholders;

//! check agreement of sizes of two vectors
// void check_match_size(const vector<T1> &first,const vector<T2> &second)
// {if(first.size()!=second.size()) CRASH("Vectors do not agree in size, %d vs %d",first.size(),second.size());}

//////////////////////////////////////////////// sum /////////////////////////////////////////////////

//! sum between two vectors
template <class T> vector<T> operator+(const vector<T> &first,const vector<T> &second)
{
  check_match_size(first,second);
  
  vector<T> out(first.size());
  for(size_t it=0;it<first.size();it++) out[it]=first[it]+second[it];
  return out;
}

//! sum between vector and a scalar
template <class TV,class TS,class = typename enable_if<!is_base_of<vector<TV>,TS>::value>::type> vector<TV> operator+(const vector<TV> &first,const TS &second)
{
  vector<TV> out(first.size());
  for(size_t it=0;it<first.size();it++) out[it]=first[it]+second;
  return out;
}

//! sum between scalar and vector
template <class TV,class TS,class = typename enable_if<!is_base_of<vector<TV>,TS>::value>::type> vector<TV> operator+(const TS &first,const vector<TV> &second)
{return second+first;}

//! summassign
template <class T1,class T2> auto operator+=(T1 &first,const T2 &second) -> decltype(first+second)
{return first=first+second;}

////////////////////////////////////////////// prod //////////////////////////////////////////////////

//! prod between two vectors
template <class T> vector<T> operator*(const vector<T> &first,const vector<T> &second)
{
  check_match_size(first,second);
  
  vector<T> out(first.size());
  for(size_t it=0;it<first.size();it++) out[it]=first[it]*second[it];
  return out;
}

//! prod between vector and a scalar
template <class TV,class TS,class = typename enable_if<!is_base_of<vector<TV>,TS>::value>::type> vector<TV> operator*(const vector<TV> &first,const TS &second)
{
  vector<TV> out(first.size());
  for(size_t it=0;it<first.size();it++) out[it]=first[it]*second;
  return out;
}

//! prod between scalar and vector
template <class TV,class TS,class = typename enable_if<!is_base_of<vector<TV>,TS>::value>::type> vector<TV> operator*(const TS &first,const vector<TV> &second)
{return second*first;}



#endif
