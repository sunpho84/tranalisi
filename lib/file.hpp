#ifndef _FILE_HPP
#define _FILE_HPP

#include <cstdio>
#include <cstring>
#include <jack.hpp>
#include <map>
#include <string>
#include <tools.hpp>
#include <vector>

using namespace std;

const int max_line_length=1024;
using line_t=char[max_line_length];

//! class to format data
template <class T> class format_str {public: static const char *value(){return "";}};
template <> class format_str<int> {public: static const char *value(){return "%d";}};
template <> class format_str<double> {public: static const char *value(){return "%lg";}};
template <> class format_str<char> {public: static const char *value(){return "%c";}};
template <> class format_str<char*> {public: static const char *value(){return "%s";}};

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
  
  //! default creator
  raw_file_t() {file=NULL;}
  
  //! creator with name
  raw_file_t(string s,string mode)
  {open(s,mode);}
  
  //! destructor
  ~raw_file_t()
  {close();}
  
  //! named or unnamed read
  template <class T> T read(const char *name=NULL)
  {
    //! unnamed read
    T out;
    int rc=fscanf(file,(name==NULL)?format_str<T>::value():((string)name+format_str<T>::value()).c_str(),&out);
    if(rc!=1) CRASH("Unbale to read %s with name \"%s\" from file",format_str<T>::value(),name);
    return out;
  }
  
  //! read with check
  template <class T> void read(T &out,const char *name=NULL)
  {out=read<T>(name);}
  
  //! read a line and eliminate trailing new line
  char* get_line(line_t line)
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
	if(rc!=NULL) tok=strtok(line," ");
	
	//skip blank line and comments
	if(tok!=NULL && strcasecmp(tok,"#"))
	  {
	    //parse the line
	    size_t nread_col=0;
	    vector<double> temp(ntot_col);
	    do
	      {
		//read a double from the line
		if(sscanf(tok,"%lg",&temp[nread_col])!=1) CRASH("Parsing col %d",nread_col);
		//check not exceeding ntot_col
		if(nread_col>=ntot_col) CRASH("nread_col=%d exceeded ntot_col=%d",nread_col,ntot_col);
		
		//search next tok
		tok=strtok(NULL," ");
		nread_col++;
	      }
	    while(tok);
	    
	    //store
	    if(nread_col==ntot_col)
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
    if(iline<nlines)
      {
	data.clear();
	//cout<<"Failed reading, invalidating"<<endl;
      }
    
    return data;
  }
};

//! read from a set of confs
djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,int nlines);

#endif
