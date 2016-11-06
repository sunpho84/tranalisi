#ifndef _BOOT_HPP
#define _BOOT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>
#include <vector>

using namespace std;

//! standard number of bootstrap sample, if not specified
int def_nboots=DEF_NBOOTS;

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

#endif
