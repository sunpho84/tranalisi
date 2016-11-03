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
  boot_t(int nboots) : vector<T>(nboots,0) {}
  
  //! default value of nboots used
  boot_t() : boot_t(def_nboots) {}
  
  //! move constructor
  boot_t(boot_t&&) = default;
  
  //! copy constructor
  boot_t(boot_t &oth) {cout<<"copia"<<endl;}
};

#endif
