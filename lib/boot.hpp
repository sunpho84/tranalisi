#ifndef _BOOT_HPP
#define _BOOT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <iostream>
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
EXTERN_BOOT int def_nboots INIT_TO(DEF_NBOOTS);

//! define the square of a double, float or integer
template <class T,class=typename enable_if<is_arithmetic<T>::value>::type> T sqr(T x)
{return x*x;}

/////////////////////////////////////////////////////////////// average and error /////////////////////////////////////////////////

//! average and error
class ave_err_t : pair<double,double>
{
public:
  //! rebind base constructor
  ave_err_t(double a=0,double b=0) : pair<double,double>(a,b) {};
  
  //! move constructor
  ave_err_t(ave_err_t&& oth)=default;
  
  //! copy constructor
  ave_err_t(const ave_err_t &oth)=default;
  
  //! move assignement
  ave_err_t &operator=(ave_err_t &&)=default;
  
  //! copy assignement
  ave_err_t &operator=(const ave_err_t &oth)
  {pair<double,double>::operator=(oth);return *this;}
  
  //! rebind average
  double &ave=first;
  
  //! rebind error
  double &err=second;
};

////////////////////////////////////////////////////// type to initialize a boot_t //////////////////////////////////////////

//! class to initialize a boot_t
class boot_init_t : public vector<int>
{
public:
  //! initialize with a given number of bootstrap
  explicit boot_init_t(int nboots=def_nboots) : vector<int>(nboots) {}
  
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

//////////////////////////////////////////////////////////////// boot_t /////////////////////////////////////////////////////

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

//! type defining boot
template <class T> class boot_t : public vector<T>
{
public:
  //! constrcutor specifying nboots
  explicit boot_t(int nboots) : vector<T>(nboots,0) {}
  
  //! constrcutor specifying gauss_filler
  explicit boot_t(const gauss_filler_t &gf) : boot_t() {fill_gauss(gf);}
  
  //! constrcutor specifying nboots and gauss_filler
  explicit boot_t(int nboots,const gauss_filler_t &gf) : boot_t(nboots) {fill_gauss(gf);}
  
  //! default value of nboots used
  boot_t() : boot_t(def_nboots) {}
  
  //! move constructor
  boot_t(boot_t&& oth) : vector<T>(forward<vector<T>>(oth)) {cout<<"move const"<<endl;}
  
  //! copy constructor
  boot_t(const boot_t &oth) : vector<T>(oth) {cout<<"copy const"<<endl;}
  
  //! move assignement
  boot_t &operator=(boot_t &&oth)=default;
  
  //! copy assignement
  boot_t &operator=(const boot_t &oth)// =default;
  {vector<T>::operator=(oth);cout<<"copy"<<endl;return *this;}
  
  //! assign from a scalar
  boot_t& operator=(const T &oth) {for(auto &it : *this) it=oth;return *this;}
  
  //! compute average and error
  ave_err_t ave_err()
  {
    ave_err_t ae;
    double &ave=ae.ave;
    double &err=ae.err;
    
    for(auto & x : *this)
    {
      ave+=x;
      err+=sqr(x);
    }
    ave/=this->size();
    err/=this->size();
    err-=sqr(ave);
    err=sqrt(fabs(err)*(njacks-1));
    
    return ae;
  }
  
  //! initialize from aver_err_t and a seed
  void fill_gauss(const gauss_filler_t &gf)
  {
    check_njacks_init();
    gen_t gen(gf.seed);
    for(auto &it : *this) it=gen.get_gauss(gf.ae.ave,gf.ae.err/sqrt(njacks-1));
  }
  
  //! initialize from ave and err
  void fill_gauss(double ave,double err,int seed)
  {fill_gauss(ave_err_t(ave,err),seed);}
  
  //! initialize from a set of jackknife
  void init_from_jackknives(const boot_init_t &iboot,jack_t &jack)
  {
    for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
    out.data[nboot]=in.data[njack];
  }
};

//! typically we will use double numbers
using dboot_t=boot_t<double>;

////////////////////////////////////////////////////////// vector of boot_t /////////////////////////////////////////////////

//! type defining boot vec
template <class T> class bvec_t : public vector<boot_t<T>>
{
public:
  //! constrcutor specifying nel and nboots
  explicit bvec_t(int nel,int nboots) : vector<boot_t<T>>(nel,boot_t<T>(nboots)) {}
  
  //! constrcutor specifying nel only (avoid copy constructor)
  explicit bvec_t(int nel=0) : vector<boot_t<T>>(nel) {}
  
  //! move constructor
  bvec_t(bvec_t&& oth) : vector<boot_t<T>>(forward<vector<boot_t<T>>>(oth)) {cout<<"vec move const"<<endl;}
  
  //! copy constructor
  bvec_t(const bvec_t &oth) : vector<boot_t<T>>(oth) {cout<<"vec copy const"<<endl;}
  
  //! move assignement
  bvec_t &operator=(bvec_t &&oth)=default;
  
  //! copy assignement
  bvec_t &operator=(const bvec_t &oth)// =default;
  {vector<boot_t<T>>::operator=(oth);cout<<"vec copy"<<endl;return *this;}
  
  //! assign from a scalar
  bvec_t& operator=(const T &oth) {for(auto &it : *this) it=oth;return *this;}
};

//! typically we use double boot
using dbvec_t=bvec_t<double>;

#undef EXTERN_BOOT
#undef INIT_TO

#endif
