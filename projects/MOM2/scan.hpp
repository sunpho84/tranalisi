#ifndef _SCAN_HPP
#define _SCAN_HPP

#ifndef EXTERN_SCAN
 #define EXTERN_SCAN extern
 #define INIT_SCAN_TO(...)
#else
 #define INIT_SCAN_TO(...) __VA_ARGS__
#endif

#include <fstream>

#include <tranalisi.hpp>

using namespace std;

int parser_parse(void*);
  
namespace scan
{
  //! opaque flex object
  EXTERN_SCAN void *scanner;
  
  //! stream from which to read
  EXTERN_SCAN ifstream fin;
  
  void init_scanner();
  void destroy_scanner();
  
  inline void scan(const string &path)
  {
    cout<<"Opening \""<<<<path<<"\""<<endl;
    
    fin.open(path);
    if(not fin.good()) CRASH("Opening %s",path.c_str());
    
    init_scanner();
    parser_parse(nullptr);
    destroy_scanner();
  }
}

#undef INIT_SCAN_TO
#undef EXTERN_SCAN

#endif
