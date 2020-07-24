#ifndef _OBS_FILE_HPP
#define _OBS_FILE_HPP

#include <raw_file.hpp>

///////////////////////////////////////////////////// file reading observables /////////////////////////////////////

//! open a file for reading, skip commented lines, put columns one after the other
class obs_file_t : public raw_file_t
{
  //! total number of columns
  size_t ntot_cols;
  
  //! view of cols
  vector<size_t> cols;
  
public:
  //! init
  obs_file_t(size_t ntot_cols=1) : ntot_cols(ntot_cols) {set_col_view(vector<size_t>{0});}
  
  //! init reading
  obs_file_t(const string &path,size_t ntot_cols=1,const vector<size_t> &ext_cols={0}) : obs_file_t(ntot_cols)
  {
    set_col_view(ext_cols);
    open(path);
  }
  
  //! set the view on columns
  void set_col_view(const vector<size_t> &ext_cols)
  {
    if(cols.size()>=ntot_cols) CRASH("Col=%zu >= NtotCol=%zu",cols.size(),ntot_cols);
    cols=ext_cols;
  }
  
  //! measure according to ncols
  size_t length(size_t nlines)
  {
    //read the position and go to the head
    long ori=get_pos();
    go_to_head();
    
    //read until the end
    size_t out=0,cur_length;
    do
      {
	cur_length=this->read(nlines,false).size();
	out+=cur_length;
      }
    while(cur_length);
    
    //go back to previous position
    set_pos(ori);
    
    return out;
  }
  
  //! set the view on columns
  void set_col_view(size_t icol)
  {set_col_view(vector<size_t>{icol});}
  
  //! open for reading obs
  void open(const string &path)
  {raw_file_t::open(path,"r");}
  
  //! read
  vector<double> read(size_t nlines=1,bool verbosity=VERBOSE)
  {
    //returned obj
    vector<double> data(cols.size()*nlines);
    
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
	if(tok!=NULL and strcasecmp(tok,"#"))
	  {
	    //parse the line
	    size_t nread_col=0;
	    vector<double> temp(ntot_cols);
	    do
	      {
		//read a double from the line
		int rc=sscanf(tok,"%lg",&temp[nread_col]);
		if(rc!=1) CRASH("Parsing col %d, rc %d from %s, line %s",nread_col,rc,tok,line);
		//check not exceeding ntot_cols
		if(nread_col>=ntot_cols) CRASH("nread_col=%d exceeded ntot_cols=%d",nread_col,ntot_cols);
		
		//search next tok
		tok=strtok_r(NULL," \t",&saveptr);
		nread_col++;
	      }
	    while(tok);
	    
	    //store
	    for(size_t icol=0;icol<cols.size();icol++)
	      data[iline+nlines*icol]=temp[cols[icol]];
	    
	    iline++;
	  }
      }
    while(!feof() and iline<nlines);
    
    //invalidate failed reading
    if(iline<nlines)
      {
	if(verbosity==VERBOSE) cout<<"Read "<<iline<<" lines instead of "<<nlines<<" in file "<<get_path()<<endl;
	data.clear();
      }
    
    return data;
  }
};

#endif
