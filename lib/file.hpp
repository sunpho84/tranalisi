#ifndef _FILE_HPP
#define _FILE_HPP

#include <cstdio>
#include <cstring>
#include <iostream>
#include <jack.hpp>
#include <map>
#include <macros.hpp>
#include <string>
#include <tools.hpp>
#include <typeinfo>
#include <vector>

using namespace std;

const size_t max_word_length=128;
const size_t max_line_length=1024;
using word_t=char[max_word_length];
using line_t=char[max_line_length];

//! class to format data
template <class T> class format_str {public: static const enable_if<is_void<T>::value,char> *value(){return "";}};
template <> class format_str<int> {public: static const char *value(){return "%d";}};
template <> class format_str<size_t> {public: static const char *value(){return "%zu";}};
template <> class format_str<double> {public: static const char *value(){return "%lg";}};
template <> class format_str<char> {public: static const char *value(){return "%c";}};
template <> class format_str<string> {public: static const char *value(){return "%s";}};

//! class to handle return type from reader
template<class T> struct return_type {typedef T type;};
template<size_t N> struct return_type<char[N]> {typedef string type;};

//! class to handle working type from reader
template<class T> struct working_type {typedef T type;};
template<> struct working_type<string> {typedef word_t type;};

//! class to handle working type from reader
template<class T> T *get_ptr(T &in) {return &in;}
template<class T> T *get_ptr(T *in) {return in;}
template<class T,size_t n> T *get_ptr(T in[n]) {return in;}

//! open file, basic
class raw_file_t
{
  FILE *file;
  
public:
  //! open the file with error check
  void open(string s,string mode)
  {
    file=fopen(s.c_str(),mode.c_str());
    if(file==NULL) CRASH("Unable to open %s with mode %s",s.c_str(),mode.c_str());
  }
  
  //! close the file
  void close()
  {if(file) fclose(file);}
  
  //! default constructor
  raw_file_t() {file=NULL;}
  
  //! copy constructor
  raw_file_t(const raw_file_t& oth)=default;
  
  //! move constructor
  raw_file_t(raw_file_t&& oth)=default;
  
  //! creator with name
  raw_file_t(string s,string mode)
  {open(s,mode);}
  
  //! destructor
  ~raw_file_t()
  {close();}
  
  //! default
  //template <class T,typename=void> void bin_write(const T &out);
  
  //! binary write, non-vector case
  template <class T> auto bin_write(const T &out) const -> enable_if_t<is_pod<T>::value>
  {if(fwrite(&out,sizeof(T),1,file)!=1) CRASH("Writing to file");}
  
  //! specialization for vector
  template <class T> auto bin_write(const T &out) const -> enable_if_t<is_vector<T>::value>
  {for(auto &it : out) bin_write(it);}
  
  //! check that the token is found
  void expect(const char *tok)
  {
    word_t rea;
    int rc=fscanf(file,"%s",rea);
    if(rc!=1) CRASH("Obtained %d while expecting token %s",rc,tok);
    if(strcasecmp(tok,rea)) CRASH("Obtained %s while expecting %s",rea,tok);
    //cout<<"Discarding "<<rea<<endl;
  }
  
  //! named or unnamed read
  template <class T> typename return_type<T>::type read(const char *name=NULL)
  {
    //check tag
    if(name) expect(name);
    
    //read a type
    typename working_type<T>::type out;
    int rc=fscanf(file,format_str<T>::value(),get_ptr(out));
    
    //check and return
    if(rc!=1) CRASH("Unable to read %s (%s) with name \"%s\" from file",format_str<T>::value(),typeid(T).name(),name);
    return (typename return_type<T>::type)out;
  }
  
  //! read with check
  template <class T> void read(T &out,const char *name=NULL)
  {out=read<T>(name);}
  
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
  
  //! return whether we are at the end of file
  bool feof()
  {return std::feof(file);}
};

///////////////////////////////////////////////////// file reading observables /////////////////////////////////////

//! open a file for reading, skip commented lines, put columns one after the other
class obs_file_t : public raw_file_t
{
  //! total number of columns
  size_t ntot_col;
  
  //! number of visible columns
  size_t nvis_col;
  
  //! view of cols
  map<size_t,vector<size_t>> col_contr;
  
public:
  //! init
  obs_file_t(size_t ntot_col=1) : ntot_col(ntot_col) {set_col_view(vector<size_t>{0});}
  
  //! init reading
  obs_file_t(const char *path,size_t ntot_col=1) : obs_file_t(ntot_col)
  {
    set_col_view(vector<size_t>{0});
    open(path);
  }
  
  //! set the view on columns
  void set_col_view(const vector<size_t> &cols)
  {
    col_contr.clear();
    nvis_col=cols.size();
    
    //push back in the list to which the col contribute
    for(size_t it=0;it<cols.size();it++)
      {
	if(cols[it]>=ntot_col) CRASH("Col=%d >= NtotCol=%d",cols[it],ntot_col);
	col_contr[cols[it]].push_back(it);
      }
  }
  
  //! set the view on columns
  void set_col_view(size_t icol)
  {set_col_view(vector<size_t>{icol});}
  
  //! open for reading obs
  void open(const char *path)
  {raw_file_t::open(path,"r");}
  
  //! read
  vector<double> read(size_t nlines=1)
  {
    //returned obj
    vector<double> data(nvis_col*nlines);
    
    size_t iline=0;
    do
      {
	//read a line
	line_t line;
	char *rc=get_line(line);
	
	//feed the line to the tokenizer
	char *tok=NULL;
	char *saveptr;
	if(rc!=NULL) tok=strtok_r(line," \t",&saveptr);
	
	//skip blank line and comments
	if(tok!=NULL && strcasecmp(tok,"#"))
	  {
	    //parse the line
	    size_t nread_col=0;
	    vector<double> temp(ntot_col);
	    do
	      {
		//read a double from the line
		int rc=sscanf(tok,"%lg",&temp[nread_col]);
		if(rc!=1) CRASH("Parsing col %d, rc %d from %s, line %s",nread_col,rc,tok,line);
		//check not exceeding ntot_col
		if(nread_col>=ntot_col) CRASH("nread_col=%d exceeded ntot_col=%d",nread_col,ntot_col);
		
		//search next tok
		tok=strtok_r(NULL," \t",&saveptr);
		nread_col++;
	      }
	    while(tok);
	    
	    //store
	    if(nread_col==nvis_col)
	      {
		for(auto &col_list : col_contr)
		  for(auto &icol : col_list.second)
		    data[icol+nvis_col*iline]=temp[col_list.first];
		iline++;
	      }
	  }
      }
    while(!feof() && iline<nlines);
    
    //invalidate failed reading
    if(iline<nlines) data.clear();
    
    return data;
  }
};

//////////////////////////////////////////////// reading input files /////////////////////////////////////////

//! open a file for reading, skip commented lines, put columns one after the other
class input_file_t : public raw_file_t
{
public:
  //! construct with a path
  input_file_t(const char *path) : raw_file_t(path,"r") {};
};

#endif
