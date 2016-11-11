#ifndef _GRACE_HPP
#define _GRACE_HPP

#include <fstream>
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
};

//! write a vector
template <class T> // enable_if_t<has_method_ave_err<T>::value,grace_file_t>
grace_file_t
&operator<<(grace_file_t &out,const vector<T> &data)
{
  out<<"@type xydy"<<endl;
  for(size_t i=0;i<data.size();i++) out<<i<<" "<<data[i].ave_err()<<endl;
  return out;
}

#endif
