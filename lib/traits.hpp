#ifndef _TRAITS_HPP
#define _TRAITS_HPP

#include <array>
#include <type_traits>
#include <valarray>
#include <vector>

using namespace std;

#if __cplusplus<201402L
 template <bool B,class T=void> using enable_if_t=typename enable_if<B,T>::type;
#endif

//! return a dummy vector traits
template <class T>
class false_vector_traits
{
public:
  typedef void base_type;
  static const bool is_vector=false;
};

//! true vector traits
template <class T>
class true_vector_traits
{
public:
  typedef T base_type;
  static const bool is_vector=true;
};

//! return a dummy vector traits
template <class T> class vector_traits : public false_vector_traits<T> {};

//! traits for array
template <class TS,size_t N> class vector_traits<array<TS,N>> : public true_vector_traits<TS> {};

//! traits for valarray
template <class TS> class vector_traits<valarray<TS>> : public true_vector_traits<TS> {};

//! traits for vector
template <class TS> class vector_traits<vector<TS>> : public true_vector_traits<TS> {};

//! helper for is_vector
template <class T> class is_vector
{
public:
  static const bool value=vector_traits<T>::is_vector;
};

//! helper for base_type
template <class T> using base_type_t=typename vector_traits<T>::base_type;

//! true only if vector
template <class TV> using enable_if_vector=enable_if_t<is_vector<TV>::value>;

//! true only if vector of double
template <class TV> using enable_if_vector_of_double=enable_if_t<is_vector<TV>::value and is_same<typename vector_traits<TV>::base_type,double>::value>;

#endif
