#ifndef _BOOT_HPP
#define _BOOT_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <vector>

using namespace std;

int nboots=DEF_NBOOTS; //<standard number of bootstrap sample, if not specified

template <class T> class boot_t : vector<T>
{
};

#endif
