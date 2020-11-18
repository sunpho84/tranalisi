#ifndef _JACK_HPP
#define _JACK_HPP

#include <ave_err.hpp>
#include <complex>
#include <raw_file.hpp>
#include <fstream>
#include <iostream>
#include <random.hpp>
#include <sstream>
#include <tools.hpp>
#include <valarray>

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

template <class T> class jack_t : public valarray<T>
{
public:
  //! base type of the jack
  typedef T base_type;
  
  //! creator
  jack_t() : valarray<T>(njacks+1) {check_njacks_init();}
  
  //! create with size (only njacks is accepted)
  //explicit jack_t(size_t ext_njacks) : jack_t() {if(njacks!=ext_njacks) CRASH("NJacks %zu different from global value %zu",ext_njacks,njacks);}
  
  //! create from double
  jack_t(double ext) : jack_t() {*this=ext;}
  
  //! create from sliced array
  jack_t(const slice_array<jack_t> &slice) : valarray<T>(slice) {}
  
  //! creator from data
  jack_t(const valarray<T> &data) : jack_t() {init_from_data(data);}
  
  //! move constructor
  jack_t(jack_t&& oth) : valarray<T>(forward<valarray<T>>(oth)) {// cout<<"move const"<<endl;
  }
  
  //! construct from expr
  template<class _Dom>
  jack_t(const _Expr<_Dom,T> &oth) : valarray<T>(oth) {}
  
  //! constrcutor specifying gauss_filler
  explicit jack_t(const gauss_filler_t &gf) : jack_t() {fill_gauss(gf);}
  
  //! copy constructor
  jack_t(const jack_t &oth) : valarray<T>(oth) {// cout<<"copy const"<<endl;
  }
  
  //! move assignement
  jack_t &operator=(jack_t &&oth) noexcept =default;
  
  //! copy assignement
  jack_t &operator=(const jack_t &oth) =default;
  
  //! assignement
  template<class oth_t>
  jack_t &operator=(const oth_t &oth)
  {valarray<T>::operator=(oth);return *this;}
  
  //! fill the central with the average
  void fill_ave_with_components_ave()
  {
    (*this)[njacks]=0;
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[njacks]+=(*this)[ijack];
    (*this)[njacks]/=njacks;
  }
  
  //! compute average and error
  ave_err_t ave_err() const
  {
    ave_err_t ae=range_ave_stddev(*this,njacks);
    ae.err()*=sqrt(njacks-1);
    
#if MEAN_TYPE==DISTR_MEAN
    if(njacks==1)
      {
	ae.ave()=(*this)[njacks];
	ae.err()=0.0;
      }
    //do nothing, the previously computed is already correct
#elif MEAN_TYPE==PROP_MEAN
    ae.ave()=(*this)[njacks];
#else
     #error Unknown mean propagation
#endif
    
    return ae;
  }
  
  //! Return the error with the error
  ave_err_t err_with_err() const
  {
    return ::err_with_err(*this,njacks);
  }
  
  //! Return the skewness
  ave_err_t skewness() const
  {
    return ::skewness(*this,njacks);
  }
  
  //! return only the average
  T ave() const {return ave_err().ave();}
  
  //! return only the error
  T err() const {return ave_err().err();}
  
  //! significativity (number of sigma of difference from 0)
  T significativity() const
  {
    auto ae=ave_err();
    return fabs(ae.ave()/ae.err());
  }
  
  //! initialize from aver_err_t and a seed
  void fill_gauss(const gauss_filler_t &gf)
  {
    check_njacks_init();
    gen_t gen(gf.seed);
    for(size_t ijack=0;ijack<njacks;ijack++)
      (*this)[ijack]=gen.get_gauss(gf.ae.ave(),gf.ae.err()/sqrt(njacks-1));
    (*this)[njacks]=gf.ae.ave();
  }
  
  //! initialize from ave and err
  void fill_gauss(T ave,T err,int seed)
  {fill_gauss(gauss_filler_t(ave,err,seed));}
  
  //! intialize froma ave_err_t and seed
  void fill_gauss(const ave_err_t &ae,int seed)
  {fill_gauss(gauss_filler_t(ae,seed));}
  
  //! fill the clusters
  size_t fill_clusters(const valarray<T> &data)
  {
    //compute cluster size
    size_t clust_size=data.size()/njacks;
    if(clust_size*njacks!=data.size()) CRASH("Data size %zu is not multiple of njacks %zu",data.size(),njacks);
    for(size_t it=0;it<data.size();it++) (*this)[it/clust_size]+=data[it];
    
    return clust_size;
  }
  
  //! clusterize
  void clusterize(size_t clust_size=1)
  {
    if(clust_size==0) CRASH("clust_size is zero");
    
    //fill clusters and compute avarages
    set_to_zero((*this)[njacks]);
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[njacks]+=(*this)[ijack];
    
    //clusterize
    for(size_t ijack=0;ijack<njacks;ijack++) (*this)[ijack]=((*this)[njacks]-(*this)[ijack])/double((njacks-1)*clust_size);
    (*this)[njacks]/=clust_size*njacks;
  }
  
  //! initialize from vector of T, so to create jackknives
  void init_from_data(const valarray<T> &data)
  {
    check_njacks_init();
    clusterize(fill_clusters(data));
  }
  
  //! write to a stream
  void bin_write(const raw_file_t &out) const
  {out.bin_write(*this);}
  
  //! wrapper with name
  void bin_write(const char *path) const
  {bin_write(raw_file_t(path,"w"));}
  
  //! wrapper with name
  void bin_write(const string &path) const
  {bin_write(path.c_str());}
  
  //! read from a stream
  void bin_read(const raw_file_t &in)
  {in.bin_read(*this);}
  
  //! wrapper with name
  void bin_read(const char *path)
  {bin_read(raw_file_t(path,"r"));}
  
  //! wrapper with name
  void bin_read(const string &path)
  {bin_read(path.c_str());}
  
  //! init (as for STL containers)
  T* begin() {return &((*this)[0]);}
  const T* begin() const {return &((*this)[0]);}
  
  //! end (as for STL containers)
  T* end() {return &((*this)[0])+this->size();}
  const T* end() const {return &((*this)[0])+this->size();}
};

//! traits for jack_t
template <class TS> class vector_traits<jack_t<TS>> : public true_vector_traits<TS> {};

//! typically we use jackknives of double
using djack_t=jack_t<double>;

//! complex jackkinves
using cdjack_t=complex<jack_t<double>>;

//! return a string
template <class T>
string to_string(const jack_t<T> &obj)
{
  ave_err_t ae=obj.ave_err();
  ostringstream os;
  os<<ae.ave()<<" "<<ae.err();
  return os.str();
}

//! get the size needed to init a jack_t
template <class T>
size_t init_nel(const jack_t<T> &obj)
{return njacks;}

//! specify hot to print a jack_t
template <class T>
ostream& operator<<(ostream &out,const jack_t<T> &v)
{return out<<v.ave_err();}

//! trim a vector in such a way that its size is multiple of njacks, and return clust_size
template <class T>
size_t trim_to_njacks_multiple(vector<T> &v,bool verbosity=false)
{
  //compute max n compatible to njacks
  size_t clust_size=v.size()/njacks;
  size_t n=clust_size*njacks;
  bool to_trim=(v.size()!=n);
  
  //output
  if(verbosity)
    {
      if(to_trim) cout<<"Trimmed from "<<v.size()<<" to "<<n<<", clust_size="<<clust_size<<endl;
      else        cout<<"No need to trim, keeping nconfs="<<n<<endl;
    }
  
  //trim if needed
  if(to_trim) v.resize(n);
  
  return clust_size;
}

//! clusterize a generic vector
template <class T>
void clusterize(vector<T> &v,size_t clust_size=1)
{
  //compute avarages
  v[njacks]=0.0;
  for(size_t ijack=0;ijack<njacks;ijack++) v[njacks]+=v[ijack];
  
  //clusterize
  for(size_t ijack=0;ijack<njacks;ijack++) v[ijack]=(v[njacks]-v[ijack])/double((njacks-1)*clust_size);
  v[njacks]/=clust_size*njacks;
}

#undef EXTERN_JACK
#undef INIT_TO

#endif
