#ifndef _JACK_HPP
#define _JACK_HPP

#ifndef EXTERN_JACK
 #define EXTERN_JACK extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#include <file.hpp>
#include <fstream>
#include <iostream>
#include <tools.hpp>
#include <vector>

using namespace std;

//! number of jackknife
#define UNDEF_NJACKS 0
EXTERN_JACK size_t njacks INIT_TO(UNDEF_NJACKS);

//! set the number of jackknives
void set_njacks(int ext_njacks);

//! crash if number of jackknives is not initialized
void check_njacks_init();

template <class T> class jack_t : public vector<T>
{
public:
  //! creator
  jack_t() : vector<T>(njacks+1) {check_njacks_init();}
  
  //! creator from data
  jack_t(const vector<T> &data) : jack_t() {init_from_data(data);}
  
  void init_from_data(const vector<T> &data)
  {
    check_njacks_init();
    
    //compute cluster size
    size_t clust_size=data.size()/njacks;
    if(clust_size*njacks!=data.size()) CRASH("Data size=%d, njacks=%d are incommensurable",data.size(),njacks);
    
    //hold clusters
    vector<T> clust(njacks+1,0);
    
    //fill clusters and compute avarages
    (*this)[njacks]=0;
    for(size_t it=0;it<data.size();it++)
      {
	clust[it/clust_size]+=data[it];
	(*this)[njacks]+=data[it];
      }
    
    //clusterize
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[ijack]=((*this)[njacks]-clust[ijack])/((njacks-1)*clust_size);
    (*this)[njacks]/=data.size();
  }
};

//! typically we use jackknives of double
using djack_t=jack_t<double>;

/////////////////////////////////////////////////////////// vector of jackknives ////////////////////////////////////////////////

//! hold a triplet to set a filter
class filter_t
{
public:
  size_t each;
  size_t offset;
  size_t howmany;
  
  filter_t(size_t each,size_t offset=0,size_t howmany=1) : each(each),offset(offset),howmany(howmany) {}
  filter_t() {}
};

template <class T> class jvec_t : public vector<jack_t<T>>
{
 public:
  jvec_t()=default;
  
  jvec_t(const vector<vector<T>> o) : vector<jack_t<T>>(o.size())
    {for(size_t it=0;it<o.size();it++) (*this)[it]=o[it];}
  
  //! filter vector
  jvec_t filter(filter_t filter)
  {
    jvec_t out;
    for(size_t it=filter.offset;it<=this->size()-filter.howmany;it+=filter.each)
      for(size_t sh=0;sh<filter.howmany;sh++)
	out.push_back((*this)[it+sh]);
    
    return out;
  }
  //! wrapper
  jvec_t filter(size_t each,size_t offset=0,size_t howmany=1)
  {return filter(filter_t(each,offset,howmany));}
  
  //! write to a stream
  void bin_write(const raw_file_t &out)
  {
    out.bin_write<size_t>(njacks);
    out.bin_write(this->size());
    out.bin_write(*this);
  }
  
  //! wrapper with name
  void bin_write(const char *path)
  {bin_write(raw_file_t(path,"w"));}
  
  //! wrapper with name
  void bin_write(const string &path)
  {bin_write(path.c_str());}
};

using djvec_t=jvec_t<double>;

//! read from a set of confs
djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,int nlines);

#undef EXTERN_JACK
#undef INIT_TO

#endif
