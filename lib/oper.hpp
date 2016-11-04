#ifndef _OPER_HPP
#define _OPER_HPP

#include <algorithm>
#include <functional>
#include <vector>

using namespace std;
using namespace placeholders;

//! sum between two vectors
template <class T> vector<T> operator+(const vector<T> &first,const vector<T> &second)
{
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

#endif
