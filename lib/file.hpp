#ifndef _FILE_HPP
#define _FILE_HPP

#include <vector>
#include <tools.hpp>

using namespace std;

const size_t max_word_length=128;
const size_t max_line_length=1024;
using word_t=char[max_word_length];
using line_t=char[max_line_length];

//! class to format data
template <class T> class format_str {public: static const enable_if<is_void<T>::value,char> *value(){return "";}};
template <> class format_str<bool> {public: static const char *value(){return "%d";}};
template <> class format_str<int> {public: static const char *value(){return "%d";}};
template <> class format_str<size_t> {public: static const char *value(){return "%zu";}};
template <> class format_str<double> {public: static const char *value(){return "%lg";}};
template <> class format_str<char> {public: static const char *value(){return "%c";}};
template <> class format_str<char*> {public: static const char *value(){return "%s";}};
template <> class format_str<string> {public: static const char *value(){return "%s";}};

//! class to handle return type from reader
template<class T> struct return_type {typedef T type;};
template<size_t N> struct return_type<char[N]> {typedef string type;};

//! class to handle working type from reader
template<class T> struct working_type {typedef T type;};
template<> struct working_type<string> {typedef word_t type;};
template<> struct working_type<bool> {typedef int type;};

//! class to handle working type from reader
template<class T> T *get_ptr(T &in) {return &in;}
template<class T> T *get_ptr(T *in) {return in;}
template<class T,size_t n> T *get_ptr(T in[n]) {return in;}

//! get a vector of all id pointing to existing paths, on the basis of a certain template and range
vector<size_t> get_existing_paths_in_range(const string &template_path,const range_t &range,bool verbosity=false);

#endif
