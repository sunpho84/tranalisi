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
void set_njacks(int ext_njacks);

//! crash if number of jackknives is not initialized
void check_njacks_init();

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
  
  //! initialize from vector of double, so to create jackknives
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

//! get the size needed to init a jack_t
template <class T> const size_t init_nel(const jack_t<T> &obj)
{return njacks;}

#undef EXTERN_JACK
#undef INIT_TO

#endif
