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
  EXTERN_SCAN vector<pair<ifstream,string>> fin;
  
  void init_scanner();
  void destroy_scanner();
  
  //! init the scanner
  void scan(const string& path);
  
  //! open a subfile
  void include(const string& path);
  
  //! print the list of open file
  void print_input_files_stack();
  
  //! read a character
  void lex(char* buf,int &result,int max_size);
}

#undef INIT_SCAN_TO
#undef EXTERN_SCAN

#endif
