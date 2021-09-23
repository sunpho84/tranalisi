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
  if(njacks==UNDEF_NJACKS)
    njacks=ext_njacks;
  else
    CRASH("Unbale to set njacks twice");
}

template <typename T>
struct jack_t;

template <typename U>
decltype(auto) getJack(U&& u,
		       const size_t& ijack)
{
  return
    u;
}

template <typename U>
U& getJack(jack_t<U>& u,
	   const size_t& ijack)
{
  return
    u[ijack];
}

template <typename U>
const U& getJack(const jack_t<U>& u,
		 const size_t& ijack)
{
  return
    u[ijack];
}

//! crash if number of jackknives is not initialized
inline void check_njacks_init()
{
  if(njacks==UNDEF_NJACKS)
    CRASH("Set njacks before");
}

template <typename T>
struct jack_t
{
  vector<T> data;
  
  //! base type of the jack
  typedef T base_type;
  
  //! creator
  jack_t() :
    data(njacks+1)
  {
    check_njacks_init();
  }
  
  //! create with size (only njacks is accepted)
  //explicit jack_t(size_t ext_njacks) : jack_t() {if(njacks!=ext_njacks) CRASH("NJacks %zu different from global value %zu",ext_njacks,njacks);}
  
  
  // //! create from sliced array
  // jack_t(const slice_array<jack_t> &slice) :
  //   valarray<T>(slice)
  // {
  // }
  
  // //! creator from data
  // jack_t(const valarray<T> &data) :
  //   jack_t()
  // {
  //   init_from_data(data);
  // }
  
  //! move constructor
  jack_t(jack_t&& oth) :
    data(std::move(oth.data))
  {
  }
  
  // //! construct from expr
  // template<class _Dom>
  // jack_t(const _Expr<_Dom,T> &oth) :
  //   valarray<T>(oth)
  // {
  // }
  
  //! constrcutor specifying gauss_filler
  explicit jack_t(const gauss_filler_t &gf) :
    jack_t()
  {
    fill_gauss(gf);
  }
  
  //! copy constructor
  jack_t(const jack_t &oth) :
    data(oth.data)
  {
  }
  
  template <typename...Args>
  explicit jack_t(const Args&...args)
  {
    data.reserve(njacks+1);
    
    for(size_t ijack=0;ijack<=njacks;ijack++)
      data.emplace_back(getJack(args,ijack)...);
  }
  
  size_t size() const
  {
    return
      data.size();
  }
  
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)	\
  CONST T& operator[](const size_t& i) CONST	\
  {						\
    return					\
      data[i];					\
  }
  
  PROVIDE_SUBSCRIBE_OPERATOR(const);
  
  PROVIDE_SUBSCRIBE_OPERATOR(/*non const*/);
  
#undef PROVIDE_SUBSCRIBE_OPERATOR
  
  //! move assignment
  jack_t &operator=(jack_t&& oth) noexcept =default;
  
  //! copy assignment
  jack_t &operator=(const jack_t& oth) =default;
  
  //! assignment
  template <typename U>
  jack_t &operator=(const U& oth)
  {
    for(size_t ijack=0;ijack<=njacks;ijack++)
      data[ijack]=getJack(oth,ijack);
    
    return
      *this;
  }
  
#define PROVIDE_CALLABLE(CONST)			\
  template <typename...Args>			\
  auto operator()(Args&&...args) CONST		\
  {						\
    using R=							\
      decltype((*this)[0](getJack(args,/*ijack*/0)...));	\
								\
    jack_t<R> res;						\
								\
    for(size_t ijack=0;ijack<=njacks;ijack++)			\
      res[ijack]=(*this)[ijack](getJack(args,ijack)...);	\
								\
    return							\
      res;							\
  }
  
  PROVIDE_CALLABLE(const);
  
  PROVIDE_CALLABLE(/* non const*/);
  
#undef PROVIDE_CALLABLE
  
  //! fill the central with the average
  void fill_ave_with_components_ave()
  {
    (*this)[njacks]=0;
    
    for(size_t ijack=0;ijack<njacks;ijack++)
      (*this)[njacks]+=(*this)[ijack];
    
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
  T ave() const
  {
    return
      ave_err().ave();
  }
  
  //! return only the error
  T err() const
  {
    return
      ave_err().err();
  }
  
  //! significativity (number of sigma of difference from 0)
  T significativity() const
  {
    auto ae=
      ave_err();
    
    return
      fabs(ae.ave()/ae.err());
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
  {
    fill_gauss(gauss_filler_t(ave,err,seed));
  }
  
  //! intialize froma ave_err_t and seed
  void fill_gauss(const ave_err_t &ae,int seed)
  {
    fill_gauss(gauss_filler_t(ae,seed));
  }
  
  //! fill the clusters
  size_t fill_clusters(const vector<T> &data)
  {
    //compute cluster size
    const size_t clust_size=
      data.size()/njacks;
    
    if(clust_size*njacks!=data.size())
      CRASH("Data size %zu is not multiple of njacks %zu",data.size(),njacks);
    
    for(size_t it=0;it<data.size();it++)
      (*this)[it/clust_size]+=data[it];
    
    return
      clust_size;
  }
  
  //! clusterize
  void clusterize(double clust_size=1.0)
  {
    if(clust_size==0)
      CRASH("clust_size is zero");
    
    //fill clusters and compute avarages
    set_to_zero((*this)[njacks]);
    for(size_t ijack=0;ijack<njacks;ijack++)
      (*this)[njacks]+=(*this)[ijack];
    
    //clusterize
    for(size_t ijack=0;ijack<njacks;ijack++)
      (*this)[ijack]=((*this)[njacks]-(*this)[ijack])/double((njacks-1)*clust_size);
    (*this)[njacks]/=clust_size*njacks;
  }
  
  //! initialize from vector of T, so to create jackknives
  void init_from_data(const vector<T>& data)
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
void clusterize(vector<T> &v,const double clust_size=1.0)
{
  //compute avarages
  v[njacks]=0.0;
  for(size_t ijack=0;ijack<njacks;ijack++) v[njacks]+=v[ijack];
  
  //clusterize
  for(size_t ijack=0;ijack<njacks;ijack++) v[ijack]=(v[njacks]-v[ijack])/double((njacks-1)*clust_size);
  v[njacks]/=clust_size*njacks;
}

/// Use to fill a jackknife-like structure
///
/// Loops automatically on the configurations and populate the
/// cluster, assigning the correct weight to each contribution, even
/// when nNconfs cannot be divided exactly by njacks
void jackknivesFill(const size_t& nConfs,                                                                    ///< Number of confs
		    const function<void(const size_t& iConf,const size_t& iClust,const double& weight)>& f); ///< Function used to fill
		    

#define PROVIDE_BINARY_OPERATOR(SYMBOL)				\
  template <typename T>						\
  jack_t<T> operator SYMBOL(const jack_t<T>& a,			\
			    const jack_t<T>& b)			\
  {								\
    jack_t<T> c;						\
								\
    for(size_t ijack=0;ijack<=njacks;ijack++)			\
      c[ijack]=a[ijack] SYMBOL b[ijack];			\
								\
    return							\
      c;							\
  }								\
								\
  template <typename T>						\
  jack_t<T> operator SYMBOL(const jack_t<T>& a,			\
			    const T& b)				\
  {								\
    jack_t<T> c;						\
								\
    for(size_t ijack=0;ijack<=njacks;ijack++)			\
      c[ijack]=a[ijack] SYMBOL b;				\
								\
    return							\
      c;							\
  }								\
								\
  template <typename T>						\
  jack_t<T> operator SYMBOL(const T& a,				\
			    const jack_t<T>& b)			\
  {								\
    jack_t<T> c;						\
								\
    for(size_t ijack=0;ijack<=njacks;ijack++)			\
      c[ijack]=a SYMBOL b[ijack];				\
								\
    return							\
      c;							\
  }								\
								\
  template <typename T>						\
  jack_t<T>& operator SYMBOL ##=(jack_t<T>& a,			\
				 const jack_t<T>& b)		\
    {								\
      for(size_t ijack=0;ijack<=njacks;ijack++)			\
	a[ijack] SYMBOL ## = b[ijack];				\
								\
      return							\
      a;							\
    }								\
								\
  template <typename T>						\
  jack_t<T>& operator SYMBOL ##=(jack_t<T>& a,			\
				 const T& b)			\
    {								\
      for(size_t ijack=0;ijack<=njacks;ijack++)			\
	a[ijack] SYMBOL ## = b;					\
								\
      return							\
      a;							\
    }


PROVIDE_BINARY_OPERATOR(+)
PROVIDE_BINARY_OPERATOR(-)
PROVIDE_BINARY_OPERATOR(*)
PROVIDE_BINARY_OPERATOR(/)

#undef PROVIDE_BINARY_OPERATOR

#define PROVIDE_UNARY_OPERATOR(NAME)				\
  template <typename T>						\
  jack_t<T> NAME(const jack_t<T>& a)				\
  {								\
    jack_t<T> b;						\
								\
    for(size_t ijack=0;ijack<=njacks;ijack++)			\
      b[ijack]=NAME(a[ijack]);					\
    								\
    return							\
      b;							\
  }

PROVIDE_UNARY_OPERATOR(abs)
PROVIDE_UNARY_OPERATOR(exp)
PROVIDE_UNARY_OPERATOR(log)
PROVIDE_UNARY_OPERATOR(sqrt)
PROVIDE_UNARY_OPERATOR(sin)
PROVIDE_UNARY_OPERATOR(sinh)
PROVIDE_UNARY_OPERATOR(cos)
PROVIDE_UNARY_OPERATOR(cosh)
PROVIDE_UNARY_OPERATOR(tan)
PROVIDE_UNARY_OPERATOR(tanh)

#undef PROVIDE_UNARY_OPERATOR

#define PROVIDE_BINARY_FUNCTION(NAME)					\
  template <typename T>							\
  jack_t<T> NAME(const jack_t<T>& a,					\
		 const jack_t<T>& b)					\
  {									\
    jack_t<T> c;							\
    									\
    for(size_t ijack=0;ijack<=njacks;ijack++)				\
      c[ijack]=NAME(a[ijack],b[ijack]);					\
    									\
    return								\
      c;								\
  }									\
									\
  template <typename T>							\
  jack_t<T> NAME(const jack_t<T>& a,					\
		 const T& b)						\
  {									\
    jack_t<T> c;							\
    									\
    for(size_t ijack=0;ijack<=njacks;ijack++)				\
      c[ijack]=NAME(a[ijack],b);					\
    									\
    return								\
      c;								\
  }									\
									\
  template <typename T>							\
  jack_t<T> NAME(const T& a,						\
		 const jack_t<T>& b)					\
  {									\
    jack_t<T> c;							\
    									\
    for(size_t ijack=0;ijack<=njacks;ijack++)				\
      c[ijack]=NAME(a,b[ijack]);					\
    									\
    return								\
      c;								\
  }

PROVIDE_BINARY_FUNCTION(pow)
#undef PROVIDE_BINARY_FUNCTION

template <typename T>
jack_t<T> operator-(const jack_t<T>& a)
{
  jack_t<T> b;
  
    for(size_t ijack=0;ijack<=njacks;ijack++)
      b[ijack]=-a[ijack];
    
    return
      b;
}

#undef EXTERN_JACK
#undef INIT_TO

#endif
