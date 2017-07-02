#ifndef _INDEX_HPP
#define _INDEX_HPP

#include <array>
#include <iostream>
#include <string>
#include <tools.hpp>

using namespace std;

//! build index
class index_t : vector<pair<string,size_t>>
{
public:
  //! return a reference to the base type
  vector<pair<string,size_t>> &base()
  {return *this;}
  
  //! return the rank
  size_t rank() const
  {return size();}
  
  //! return the name of the  index
  string name(size_t comp) const
  {
    if(comp>=rank()) CRASH("Aksing for the name of component %zu greater or equal to rank %zu",comp,rank());
    return (*this)[comp].first;
  }
  
  //! get the full description
  string descr(size_t i)
  {
    vector<size_t> comp=(*this)(i);
    string out="";
    for(size_t i=0;i<rank();i++)
      {
	out+=name(i)+"="+to_string(comp[i]);
	if(i!=rank()-1) out+=", ";
      }
    return out;
  }
  
  //! set the range
  void set_ranges(const vector<pair<string,size_t>> &list)
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
	if((*this)[i].second<=oth.begin()[i]) CRASH("Calling rank %zu with value %zu greater or equal to max %zu",i,oth.begin()[i],(*this)[i].second);
	out=out*(*this)[i].second+oth.begin()[i];
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
	out[i]=temp%(*this)[i].second;
	temp/=(*this)[i].second;
      }
    return out;
  }
  
  //! return the maximal index+1
  size_t max() const
  {
    size_t out=1;
    for(size_t i=0;i<rank();i++) out*=((*this)[i]).second;
    return out;
  }
  
  //! return the maximal value of a component
  size_t max(size_t icomp) const
  {return ((*this)[icomp]).second;}
  
  index_t() : vector<pair<string,size_t>>() {}
  
  index_t(const vector<pair<string,size_t>> &list) : vector<pair<string,size_t>>(list) {}
};

#endif
