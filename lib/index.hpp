#ifndef _INDEX_HPP
#define _INDEX_HPP

#include <array>
#include <iostream>
#include <string>
#include <tools.hpp>

using namespace std;

//! build index
class index_t : vector<size_t>
{
public:
  //! return the rank
  const size_t rank() const
  {return size();}
  
  //! set the range
  void set_ranges(const vector<size_t> &list)
  {
    this->resize(list.size());
    copy(list.begin(),list.end(),this->begin());
  }
  
  //! get an element
  size_t operator()(const vector<size_t> &oth) const
  {
    if(oth.size()!=rank()) CRASH("Calling index operator of rank %zu with the wrong size %zu",rank(),oth.size());
    size_t out=0;
    for(size_t i=0;i<rank();i++)
      {
	if((*this)[i]<=oth.begin()[i]) CRASH("Calling rank %zu with value %zu greater or equal to max %zu",i,oth.begin()[i],(*this)[i]);
	out=out*(*this)[i]+oth.begin()[i];
      }
    return out;
  }
  
  //! decompose
  vector<size_t> operator()(const size_t &oth) const
  {
    if(oth>=max()) CRASH("Calling decompose %zu which is beyond maximum value %zu",oth,max());
    vector<size_t> out(rank());
    
    size_t temp=oth;
    for(size_t invi=0;invi<rank();invi++)
      {
	int i=rank()-1-invi;
	out[i]=temp%(*this)[i];
	temp/=(*this)[i];
      }
    return out;
  }
  
  //! return the maximal index+1
  size_t max() const
  {
    size_t out=1;
    for(size_t i=0;i<rank();i++) out*=((*this)[i]);
    return out;
  }
  
  index_t() : vector<size_t>() {}
  
  index_t(const vector<size_t> &list) : vector<size_t>(list) {}
};

#endif
