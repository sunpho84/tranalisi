#ifndef _JACK_HPP
#define _JACK_HPP

#include <ave_err.hpp>
#include <file.hpp>
#include <fstream>
#include <iostream>
#include <tools.hpp>
#include <vector>

#ifndef EXTERN_JACK
 #define EXTERN_JACK extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

using namespace std;

//! number of jackknife
#define UNDEF_NJACKS 0
EXTERN_JACK size_t njacks INIT_TO(UNDEF_NJACKS);

//! set the number of jackknives
inline void set_njacks(int ext_njacks)
{
  if(njacks==UNDEF_NJACKS) njacks=ext_njacks;
  else CRASH("Unbale to set njacks twice");
}

//! crash if number of jackknives is not initialized
inline void check_njacks_init()
{if(njacks==UNDEF_NJACKS) CRASH("Set njacks before");}

template <class T> class jack_t : public vector<T>
{
public:
  //! base type of the jack
  typedef T base_type;
  
  //! creator
  jack_t() : vector<T>(njacks+1) {check_njacks_init();}
  
  //! create with size (only njacks is accepted)
  jack_t(size_t ext_njacks) : jack_t() {if(njacks!=ext_njacks) CRASH("NJacks %zu different from global value %zu",njacks,ext_njacks);}
  
  //! creator from data
  jack_t(const vector<T> &data) : jack_t() {init_from_data(data);}
  
  //! compute average and error
  ave_err_t ave_err() const
  {
    ave_err_t ae=range_ave_stddev(*this,njacks);
    ae.err*=sqrt(njacks-1);
    
#if MEAN_TYPE==DISTR_MEAN
#elif MEAN_TYPE==PROP_MEAN
    ae.ave=(*this)[njacks];
#else
     #error Unknown mean propagation
#endif
    
    return ae;
  }
  
  //! return only the average
  double ave() const {return ave_err().ave;}
  
  //! return only the error
  double err() const {return ave_err().err;}
  
  //! fill the clusters
  size_t fill_clusters(const vector<T> &data)
  {
    //compute cluster size
    size_t clust_size=data.size()/njacks;
    if(clust_size*njacks!=data.size()) CRASH("Data size %zu is not multiple of njacks %zu",data.size(),njacks);
    for(size_t it=0;it<data.size();it++) (*this)[it/clust_size]+=data[it];
    
    return clust_size;
  }
  
  //! clusterize
  void clusterize(size_t clust_size)
  {
    //fill clusters and compute avarages
    (*this)[njacks]=0;
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[njacks]+=(*this)[ijack];
    
    //clusterize
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[ijack]=((*this)[njacks]-(*this)[ijack])/((njacks-1)*clust_size);
    (*this)[njacks]/=clust_size*njacks;
  }
  
  //! initialize from vector of double, so to create jackknives
  void init_from_data(const vector<T> &data)
  {
    check_njacks_init();
    clusterize(fill_clusters(data));
  }
};

//! typically we use jackknives of double
using djack_t=jack_t<double>;

//! get the size needed to init a jack_t
template <class T> const size_t init_nel(const jack_t<T> &obj)
{return njacks;}

#undef EXTERN_JACK
#undef INIT_TO

#endif
