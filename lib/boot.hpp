#ifndef _BOOT_HPP
#define _BOOT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <ave_err.hpp>
#include <cmath>
#include <iostream>
#include <jack.hpp>
#include <utility>
#include <random.hpp>
#include <vector>

using namespace std;

#ifndef EXTERN_BOOT
 #define EXTERN_BOOT extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

//! standard number of bootstrap sample, if not specified
EXTERN_BOOT size_t def_nboots INIT_TO(DEF_NBOOTS);

////////////////////////////////////////////////////// type to initialize a boot_t //////////////////////////////////////////

//! class to initialize a boot_t
class boot_init_t : public vector<size_t>
{
public:
  //! initialize with a given number of bootstrap
  explicit boot_init_t(size_t nboots=def_nboots) : vector<size_t>(nboots) {}
  
  //! fill with a seed
  void fill(int seed)
  {
    check_njacks_init();
    
    //init the fenerator
    gen_t gen(seed);
    
    //fill all elements
    for(auto &it : *this) it=gen.get_int(0,njacks);
  }
};

///////////////////////////////////////////////////////////// gauss_filler_t /////////////////////////////////////////////////////

//! allows to fill from gauss
class gauss_filler_t : pair<ave_err_t,int>
{
public:
  //! fill from ave_err and seed
  gauss_filler_t(const ave_err_t &ext_ae,int ext_seed)
  {
    ae=ext_ae;
    seed=ext_seed;
  }
  
  //! fill from ave, err and seed
  gauss_filler_t(double ave,double err,int ext_seed) : gauss_filler_t(ave_err_t(ave,err),ext_seed) {}
  
  //! rebind ave_err
  ave_err_t &ae=first;
  
  //! rebind seed
  int &seed=second;
};

//////////////////////////////////////////////////////////////// boot_t /////////////////////////////////////////////////////

//! type defining boot
template <class T> class boot_t : public valarray<T>
{
public:
  //! base type of the jack
  typedef T base_type;
  
  //! return the number of bootstrap
  size_t nboots() const {return this->size()-1;}
  
  //! constrcutor specifying nboots
  explicit boot_t(size_t nboots) : valarray<T>((T)0.0,nboots+1) {}
  
  //! constrcutor specifying gauss_filler
  explicit boot_t(const gauss_filler_t &gf) : boot_t() {fill_gauss(gf);}
  
  //! constrcutor specifying nboots and gauss_filler
  explicit boot_t(size_t nboots,const gauss_filler_t &gf) : boot_t(nboots+1) {fill_gauss(gf);}
  
  //! constrcutor specifying iboot_ind and a jack
  explicit boot_t(const boot_init_t &boot_init,const jack_t<T> &jack) : boot_t(boot_init.size()) {fill_from_jack(boot_init,jack);}
  
  //! default value of nboots used
  boot_t() : boot_t(def_nboots) {}
  
  //! move constructor
  boot_t(boot_t&& oth) : valarray<T>(forward<valarray<T>>(oth)) {cout<<"move const"<<endl;}
  
  //! copy constructor
  boot_t(const boot_t &oth) : valarray<T>(oth) {cout<<"copy const"<<endl;}
  
  //! move assignement
  boot_t &operator=(boot_t &&oth)=default;
  
  //! copy assignement
  boot_t &operator=(const boot_t &oth)// =default;
  {valarray<T>::operator=(oth);cout<<"copy"<<endl;return *this;}
  
  //! assign from a scalar
  boot_t& operator=(const T &oth) {for(auto &it : *this) it=oth;return *this;}
  
  //! compute average and error
  ave_err_t ave_err() const
  {
    ave_err_t ae=range_ave_stddev(*this,nboots());
    ae.err*=sqrt(njacks-1);
    
#if MEAN_TYPE==DISTR_MEAN
    //do nothing, the previously computed is already correct
#elif MEAN_TYPE==PROP_MEAN
    ae.ave=(*this)[nboots()];
#else
     #error Unknown mean propagation
#endif
    
    return ae;
  }
  
  //! return only the average
  double ave() const {return ave_err().ave;}
  
  //! return only the error
  double err() const {return ave_err().err;}
  
  //! initialize from aver_err_t and a seed
  void fill_gauss(const gauss_filler_t &gf)
  {
    check_njacks_init();
    gen_t gen(gf.seed);
    for(size_t iboot=0;iboot<nboots();iboot++)
      (*this)[iboot]=gen.get_gauss(gf.ae.ave,gf.ae.err/sqrt(njacks-1));
    (*this)[nboots()]=gf.ae.ave;
  }
  
  //! initialize from ave and err
  void fill_gauss(double ave,double err,int seed)
  {fill_gauss(ave_err_t(ave,err),seed);}
  
  //! initialize from a jackknife
  void fill_from_jack(const boot_init_t &iboot_ind,const jack_t<T> &jack)
  {
    for(size_t iboot=0;iboot<nboots();iboot++)
      {
	size_t ind=iboot_ind[iboot];
	if(ind>=njacks) CRASH("Index %d not in the interval [0,%d]",ind,njacks-1);
	(*this)[iboot]=jack[ind];
      }
    (*this)[nboots()]=jack[njacks];
  }
};

//! typically we will use double numbers
using dboot_t=boot_t<double>;

//! return the size needed to init a boot_t
template <class T> size_t init_nel(const boot_t<T> &obj)
{return obj.nboots();}

#undef EXTERN_BOOT
#undef INIT_TO

#endif
