#ifndef _GRACE_HPP
#define _GRACE_HPP

#include <fstream>
#include <ave_err.hpp>
#include <tools.hpp>

using namespace std;

//! class to write a grace file
class grace_file_t : public ofstream
{
  size_t iset;
public:
  //! default constructor
  grace_file_t(const string &path) : ofstream(path),iset(0) {if(!this->good()) CRASH("Unable to open grace file %s",path.c_str());}
  
  //! shift iset
  void shift_iset(size_t how_many=1)
  {iset+=how_many;}
  
  //! start a new set
  void new_set()
  {
    (*this)<<"&"<<endl;
    shift_iset();
  }
  
  //! set a property
  void set_prop(string what){(*this)<<"@s"<<iset<<" "<<what<<endl;}
  
  //! set line style
  void set_line_style(size_t how){set_prop("line type "+to_string(how));}
  
  //! set no line
  void no_line(){set_line_style(0);}

};

//! write a vector of average and error
grace_file_t &operator<<(grace_file_t &out,const vec_ave_err_t &data);

#endif
