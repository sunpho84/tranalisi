#ifndef _BOOT_HPP
#define _BOOT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <iostream>
#include <utility>
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

//! number of jacknife
EXTERN_BOOT int njacks INIT_TO(1);

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
  void fill(int seed);
};

//////////////////////////////////////////////////////////////// boot_t /////////////////////////////////////////////////////

//! type defining boot
template <class T> class boot_t : public vector<T>
{
public:
  //! constrcutor specifying nboots
  explicit boot_t(int nboots) : vector<T>(nboots,0) {}
  
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
