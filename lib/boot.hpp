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

template <class T> class boot_t : public vector<T>
{
public:
  //! constrcutor specifying nboots
  explicit boot_t(int nboots) : vector<T>(nboots,0) {}
  
  //! default value of nboots used
  boot_t() : boot_t(def_nboots) {}
  
  //! init from vector<T>
  boot_t(const vector<T> oth) : vector<T>(oth) {cout<<"Conversion from vector"<<endl;}
  
  //! move constructor
  boot_t(const boot_t&& oth) : vector<T>(move(oth)) {};
  
  //! copy constructor
  boot_t(const boot_t &oth) : vector<T>(oth) {cout<<"copia"<<endl;}
  
  //! move assignement
  boot_t &operator=(const boot_t &&oth) {vector<T>(move(oth));cout<<"move"<<endl;return *this;}
  
  //! assign from a scalar
  boot_t& operator=(const T &oth) {for(auto &it : *this) it=oth;return *this;}
};

#endif
