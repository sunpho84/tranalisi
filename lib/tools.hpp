#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <array>
#include <macros.hpp>
#include <string>
#include <vector>

using namespace std;

//! crashes emitting the message
void internal_crash(int line,const char *file,const char *temp,...);

//! combine arguments in a single string
string combine(const char *format,...);

//!check if a file exists
int file_exists(string path);

//! check if a directoy exists
int dir_exists(string path);

//! return a range of int
class range_t : private array<size_t,3>
{
public:
  //bind
  size_t &start=(*this)[0];
  size_t &each=(*this)[1];
  size_t &end=(*this)[2];
  
  //! default constructor
  range_t() {}
  
  //! init from triplet
  range_t(initializer_list<size_t> list)
  {
    if(list.size()!=3) CRASH("list size %d while expecting 3",list.size());
    copy(list.begin(),list.end(),this->begin());
  }
};

DEFINE_HAS_METHOD(size);
#define is_vector has_method_size

#endif
