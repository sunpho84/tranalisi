#ifndef _OPER_HPP
#define _OPER_HPP

#include <algorithm>
#include <functional>
#include <vector>

using namespace std;

//! sum between two vectors
template <class T> vector<T> operator+(const vector<T> &first,const vector<T> &second)
{
  vector<T> out(first.size());
  transform(first.begin(),first.end(),second.begin(),out.begin(),plus<T>());
  return out;
}

#endif
