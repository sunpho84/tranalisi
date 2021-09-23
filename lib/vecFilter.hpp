#ifndef _VECFILTER_HPP
#define _VECFILTER_HPP

#include <index.hpp>
#include <tools.hpp>

//! filter a vector
template <class T,class=enable_if_t<has_method_size<T>::value>>
T vec_filter(const T &v,
	     const size_t offset,
	     const vector<size_t>& lengths,
	     const vector<size_t>& strides)
{
  vector<pair<string,size_t>> idIniter;
  for(size_t i=0;i<lengths.size();i++)
    idIniter.push_back({"id"+to_string(i),lengths[i]});
  
  index_t id;
  id.set_ranges(idIniter);
  
  T res(id.max());
  
  for(size_t o=0;o<id.max();o++)
    {
      const vector<size_t> c=
	id(o);
      
      size_t i=
	offset;
      
      for(size_t d=0;d<strides.size();d++)
	i=c[d]*strides[d];
      
      res[o]=v[i];
    }
  
  return
    res;
}

#endif
