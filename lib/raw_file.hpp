#ifndef _RAW_FILE_HPP
#define _RAW_FILE_HPP

#include <file.hpp>

#include <complex>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <map>
#include <macros.hpp>
#include <string>
#include <tools.hpp>
#include <typeinfo>
#include <vector>

using namespace std;

//! open file, basic
class raw_file_t
{
  //! handle
  FILE *file;
  
  //! stored path
  string path;
  
  //! copy constructor - made private, as for any class stream
  raw_file_t(const raw_file_t &oth);
  
public:
  //! open the file with error check
  void open(const string &_path,const string &mode)
  {
    path=_path;
    //cout<<"Opening "<<path<<endl;
    file=fopen(path.c_str(),mode.c_str());
    if(file==nullptr)
      {
	perror("for this reason:");
	CRASH("Unable to open \"%s\" with mode %s",path.c_str(),mode.c_str());
      }
  }
  
  //! return the path
  string get_path()
  {
    return path;
  }
  
  //! close the file
  void close()
  {
    //cout<<"Should close "<<path<<": "<<(file!=nullptr)<<endl;
    if(file!=nullptr)
      {
	//cout<<"Closing "<<path<<endl;
	fclose(file);
      }
    file=nullptr;
  }
  
  //! default constructor
  raw_file_t()
  {
    file=nullptr;
  }
  
  //! move constructor
  raw_file_t(raw_file_t&& oth) : raw_file_t()
  {
    swap(file,oth.file);
    swap(path,oth.path);
  };
  
  //! creator with name
  raw_file_t(const string &_path,string mode)
  {
    open(_path,mode);
  }
  
  //! destructor
  ~raw_file_t()
  {
    close();
  }
  
  //! default
  //template <class T,typename=void> void bin_write(const T &out);
  
  //! binary write, non-vector case
  template <class T>
  auto bin_write(const T &out) const -> enable_if_t<is_pod<T>::value>
  {
    if(fwrite(&out,sizeof(T),1,file)!=1)
      CRASH("Writing to file");
  }
  
  //! specialization for vector
  template <class T>
  auto bin_write(const T &out) const -> enable_if_t<is_vector<T>::value>
  {
    for(auto &it : out)
      bin_write(it);
  }
  
  //! binary read, non-vector case
  template <class T>
  auto bin_read(T &out) const -> enable_if_t<is_pod<T>::value>
  {
    int rc=fread(&out,sizeof(T),1,file);
    if(rc!=1)
      CRASH("Reading from file %s, rc: %d",path.c_str(),rc);
  }
  
  //! binary read, complex case
  template <class T>
  auto bin_read(complex<T> &out) const -> enable_if_t<is_pod<T>::value>
  {
    for(size_t ri=0;ri<2;ri++)
      bin_read(((T*)&out)[ri]);
  }
  
  //! specialization for vector
  template <class T>
  auto bin_read(T &out) const -> enable_if_t<is_vector<T>::value>
  {
    for(auto &it : out)
      bin_read(it);
  }
  
  //! return what read
  template <class T>
  T bin_read()
  {
    T out;
    bin_read(out);
    return out;
  }
  
  //! check that the token is found
  void expect(const char *tok) const
  {
    word_t rea;
    int rc=fscanf(file,"%s",rea);
    if(rc!=1) CRASH("Obtained %d while expecting token %s",rc,tok);
    if(strcasecmp(tok,rea)) CRASH("Obtained %s while expecting %s",rea,tok);
  }
  
  //! expect a list of tokens
  void expect(const initializer_list<string> &toks) const
  {for(auto &tok : toks) expect(tok.c_str());}
  
  //! set the position to the passed value
  long get_pos() {return ftell(file);}
  
  //! set the position to the passed value
  void set_pos(long offset) {fseek(file,offset,SEEK_SET);}
  
  //! reset the position to the head
  void go_to_head() {set_pos(0);}
  
  //! reset the position to the end
  void go_to_end() {fseek(file,0,SEEK_END);}
  
  //! get the size of the file
  long size()
  {
    long ori=get_pos();
    //go to the end and get position
    go_to_end();
    long length=get_pos();
    //return
    set_pos(ori);
    
    return length;
  }
  
  //! named or unnamed read
  template <class T>
  typename return_type<T>::type read(const char *name=nullptr) const
  {
    //check tag
    if(name!=nullptr) expect(name);
    
    //read a type
    typename working_type<T>::type out;
    int rc=fscanf(file,format_str<T>::value(),get_ptr(out));
    
    //check and return
    if(rc!=1) CRASH("Unable to read %s (%s) with name \"%s\" from file",format_str<T>::value(),typeid(T).name(),name);
    return (typename return_type<T>::type)out;
  }
  
  //! read with check
  template <class T,class=enable_if_t<is_pod<T>::value>>
  void read(T &out,const char *name=nullptr) const
  {
    out=read<T>(name);
  }
  
  //! read a string
  void read(string &out,const char *name=nullptr) const
  {
    out=read<string>(name);
  }
  
  //! read a vector
  template <class TV,class TS=typename TV::base_type,class=enable_if_t<is_vector<TV>::value>>
  void read(TV &out,const char *name=nullptr) const
  {
    if(name)
      expect(name);
    for(auto &o : out)
      read(o);
  }
  
  //! read a line and eliminate trailing new line
  char *get_line(line_t line)
  {
    char *ret=fgets(line,max_line_length,file);
    if(ret)
      {
	char *pos=strchr(line,'\n');
	if(pos) *pos='\0';
      }
    
    return ret;
  }
  
  //! skip a line
  void skip_line(size_t nskip=1)
  {
    line_t skip;
    for(size_t iskip=0;iskip<nskip;iskip++) get_line(skip);
  }
  
  //! write the data with a possible tag
  template <class T>
  int write(const T &out,const char *name=nullptr) const
  {
    int rc=0;
    
    //check tag
    if(name!=nullptr) rc+=fprintf(file,format_str<T>::value(),get_ptr(out));
    
    //write a type
    rc+=fprintf(file,format_str<T>::value(),get_ptr(out));
    
    return rc;
  }
  
  //! write a string
  int write(const char *in,const char *name=nullptr) const
  {
    return write<string>(in,name);
  }
  
  //! return whether we are at the end of file
  bool feof()
  {
    return std::feof(file);
  }
  
  //! formatted print
  int printf(const char *format,...)
  {
    va_list ap;
    va_start(ap,format);
    int ret=vfprintf(file,format,ap);
    va_end(ap);
    
    return ret;
  }
};

//! return the size of the passed path
inline long file_size(const string &path)
{
  return raw_file_t(path,"r").size();
}

//! open a file for reading, skip commented lines, put columns one after the other
class input_file_t : public raw_file_t
{
public:
  //! construct with a path
  input_file_t(const string &path) : raw_file_t(path,"r") {};
};

#endif
