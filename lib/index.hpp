#ifndef _INDEX_HPP
#define _INDEX_HPP

#include <array>
#include <iostream>
#include <string>
#include <tools.hpp>

using namespace std;

//! build index
template <size_t rank> class index_t : array<size_t,rank>
{
public:
  //! set the range
  void set_ranges(const initializer_list<size_t> &list)
  {
    if(list.size()!=rank) CRASH("Initializing index of rank %zu with the wrong size %zu",rank,list.size());
    copy(list.begin(),list.end(),this->begin());
  }
  
  //! get an element
  size_t operator()(const initializer_list<size_t> &oth)
  {
    if(oth.size()!=rank) CRASH("Calling index operator of rank %zu with the wrong size %zu",rank,oth.size());
    size_t out=0;
    for(size_t i=0;i<rank;i++)
      {
	if((*this)[i]<=oth.begin()[i]) CRASH("Calling rank %zu with value %zu greater or equal to max %zu",i,oth.begin()[i],(*this)[i]);
	out=out*(*this)[i]+oth.begin()[i];
      }
    return out;
  }
  
  size_t size()
  {
    size_t out=1;
    for(size_t i=0;i<rank;i++) out*=((*this)[i]);
    return out;
  }
};

#endif
