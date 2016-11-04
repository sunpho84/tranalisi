#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <string>

using namespace std;

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

//! crashes emitting the message
void internal_crash(int line,const char *file,const char *temp,...);

//! combine arguments in a single string
string combine(const char *format,...);

//!check if a file exists
int file_exists(string path);

//! check if a directoy exists
int dir_exists(string path);


#endif
