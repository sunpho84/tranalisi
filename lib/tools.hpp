#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <array>
#include <string>
#include <vector>

using namespace std;

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)
#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)
#define NAME4(s1,s2,s3,s4) NAME3(CONCAT3(s1,_,s2),s3,s4)

//! crashes emitting the message
void internal_crash(int line,const char *file,const char *temp,...);

//! combine arguments in a single string
string combine(const char *format,...);

//!check if a file exists
int file_exists(string path);

//! check if a directoy exists
int dir_exists(string path);

//! return a range of int
class range_t : public array<size_t,3>
{
public:
  //bind
  size_t &start=(*this)[0];
  size_t &each=(*this)[1];
  size_t &end=(*this)[2];
  
  //init from triplet
  range_t(initializer_list<size_t> list)
  {copy(list.begin(),list.end(),this->begin());}
};

//! check if has method
#define DEFINE_HAS_METHOD(METHOD)					\
  template <typename T> class NAME2(has_method,METHOD)			\
  {									\
  private:								\
    typedef char yes;							\
    typedef yes no[2];							\
    									\
    template<typename C> static auto test(void*)->decltype(size_t{std::declval<C const>().METHOD()},yes{}); \
    template<typename> static no& test(...);				\
    									\
  public:								\
    static bool const value=sizeof(test<T>(0))==sizeof(yes);		\
  }

DEFINE_HAS_METHOD(size);
#define is_vector has_method_size

#endif
