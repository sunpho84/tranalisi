#ifndef _JACK_HPP
#define _JACK_HPP

#ifndef EXTERN_JACK
 #define EXTERN_JACK extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

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
  jack_t() : vector<T>(njacks) {check_njacks_init();}
  
  //! creator from data
  jack_t(const vector<T> &data) : jack_t() {init_from_data(data);}
  
  void init_from_data(const vector<T> &data)
  {
    check_njacks_init();
    
    //compute cluster size
    size_t clust_size=data.size()/njacks;
    if(clust_size*njacks!=data.size()) CRASH("Data size=%d, njacks=%d are incommensurable",data.size(),njacks);
    
    //hold clusters
    vector<T> clust(njacks,0);
    
    //fill clusters and compute avarages
    T tot=0;
    for(size_t it=0;it<data.size();it++)
      {
	clust[it/clust_size]+=data[it];
	tot+=data[it];
      }
    
    //clusterize
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[ijack]=(tot-clust[ijack])/((njacks-1)*clust_size);
  }
};

//! typically we use jackknives of double
using djack_t=jack_t<double>;

/////////////////////////////////////////////////////////// vector of jackknives ////////////////////////////////////////////////

// template <class T> class jvec_t : public vector<jack_t<T>>
// {
// public:
//   //! creator from data
//   jvec_t(const vector<vector<T>> &data) : vector<jack_t<T>>(data) {init_from_data(data);}
  
//   //! initialize from a vector of vectors
//   void init_from_data(const vector<vector<T>> &data)
//   {for(size_t it=0;it<data.size();it++) (*this)[it]=data
    
//   }
// };

template <class T> class jvec_t : public vector<jack_t<T>>
{
 public:
  jvec_t(const vector<vector<T>> o) : vector<jack_t<T>>(o.size())
    {for(size_t it=0;it<o.size();it++) (*this)[it]=o[it];}
};

#undef EXTERN_JACK
#undef INIT_TO

#endif
