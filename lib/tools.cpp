#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>
#include <cstdarg>
#include <sys/stat.h>

#include <tools.hpp>

using namespace std;

void internal_crash(int line,const char *file,const char *temp,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);
  
  cerr<<"ERROR at line "<<line<<" of file "<<file<<": "<<buffer<<endl;
  exit(1);
}

string combine(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  return buffer;
}

int file_exists(string path)
{
  int status=1;
  
  FILE *f=fopen(path.c_str(),"r");
  if(f!=NULL)
    {
      status=1;
      fclose(f);
    }
  else status=0;
  
  return status;
}

int dir_exists(string path)
{
  struct stat info;
  int rc=stat(path.c_str(),&info);
  int is=(info.st_mode&S_IFDIR);
  return (rc==0)&&is;
}
