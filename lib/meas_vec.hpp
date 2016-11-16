#ifndef _MEAS_VEC_HPP
#define _MEAS_VEC_HPP

#include <boot.hpp>
#include <file.hpp>
#include <tools.hpp>

//! type defining a vector of measures
template <class meas_t> class vmeas_t : public vector<meas_t>
{
public:
  //! bind the base type of meas_t
  using base_type=typename meas_t::base_type;
  
  //! constructor specifying nel only (avoid copy constructor)
  explicit vmeas_t(size_t nel=0) : vector<meas_t>(nel) {}
  
  //! init from vector of vector
  vmeas_t(const vector<vector<base_type>> &o) : vmeas_t(o.size())
  {
    auto start=take_time();
    cout<<"Initializing vmeas_t from vector of vector"<<endl;
    for(size_t it=0;it<o.size();it++) (*this)[it]=o[it];
    cout<<elapsed_time(start)<<" to init vmeas_t"<<endl;
  }
  
  //! move constructor
  vmeas_t(vmeas_t&& oth) : vector<meas_t>(forward<vector<meas_t>>(oth)) {cout<<"vec move const"<<endl;}
  
  //! copy constructor
  vmeas_t(const vmeas_t &oth) : vector<meas_t>(oth) {cout<<"vec copy const"<<endl;}
  
  //! range constructor
  template <class InputIterator> vmeas_t(InputIterator first,InputIterator last,enable_if_t<!is_integral<InputIterator>::value>* =0) : vector<meas_t>(first,last) {}
  
  //! move assignement
  vmeas_t &operator=(vmeas_t &&oth)=default;
  
  //! copy assignement
  vmeas_t &operator=(const vmeas_t &oth)// =default;
  {vector<meas_t>::operator=(oth);cout<<"vec copy"<<endl;return *this;}
  
  //! assign from a scalar
  vmeas_t& operator=(const base_type &oth)
  {for(auto &it : *this) it=oth;return *this;}
  
  //! compute average and error
  vec_ave_err_t ave_err() const
  {
    vec_ave_err_t out(this->size());
    for(size_t it=0;it<this->size();it++) out[it]=(*this)[it].ave_err();
    return out;
  }
  
  //! write to a stream
  void bin_write(const raw_file_t &out)
  {out.bin_write(*this);}
  
  //! wrapper with name
  void bin_write(const char *path)
  {bin_write(raw_file_t(path,"w"));}
  
  //! wrapper with name
  void bin_write(const string &path)
  {bin_write(path.c_str());}
  
  //! read to a stream
  void bin_read(const raw_file_t &out)
  {out.bin_read(*this);}
  
  //! wrapper with name
  void bin_read(const char *path)
  {bin_read(raw_file_t(path,"r"));}
  
  //! wrapper with name
  void bin_read(const string &path)
  {bin_read(path.c_str());}
  
  //! assign from a scalar
  vmeas_t& operator=(const meas_t &oth) {for(auto &it : *this) it=oth;return *this;}
  
  //! return the averaged
  vmeas_t symmetrized(int par)
  {
    size_t nel=this->size();
    size_t nelh=nel/2;
    
    if(abs(par)!=1) CRASH("Unknown parity %d",par);
    
    if(nel%2) CRASH("Size %zu odd",nel);
    
    //! prepare output copying the first half+1
    vmeas_t out(&((*this)[0]),&((*this)[nelh+1]));
    
    //sum the mirror
    for(size_t iel=1;iel<nelh;iel++)
      out[iel]=(out[iel]+par*((*this)[nel-iel]))/2;
    
    return out;
  }
};

//! typically we use double
using djvec_t=vmeas_t<jack_t<double>>;
using dbvec_t=vmeas_t<boot_t<double>>;

//! read a binary from file
template <class T> T read_vec_meas(raw_file_t &file,size_t nel,size_t ind=0)
{
  T out(nel);
  if(nel==0) CRASH("Nel==0");
  file.set_pos(ind*sizeof(typename T::base_type)*out[0].size()*nel);
  out.bin_read(file);
  return out;
}

//! read a binary from the path
template <class T> T read_vec_meas(string path,size_t nel,size_t ind=0)
{
  raw_file_t file(path,"r");
  return read_vec_meas<T>(file,nel,ind);
}

//! read a djvec from path
inline djvec_t read_djvec(string path,size_t nel,size_t ind=0)
{return read_vec_meas<djvec_t>(path,nel,ind);}

//! read a dbvec from path
inline dbvec_t read_dbvec(string path,size_t nel,size_t ind=0)
{return read_vec_meas<dbvec_t>(path,nel,ind);}

//! read a djvec from path
inline djvec_t read_djvec(raw_file_t &file,size_t nel,size_t ind=0)
{return read_vec_meas<djvec_t>(file,nel,ind);}

//! read a dbvec from path
inline dbvec_t read_dbvec(raw_file_t &file,size_t nel,size_t ind=0)
{return read_vec_meas<dbvec_t>(file,nel,ind);}

//! initialize from a jvec
template <class T> vmeas_t<boot_t<T>> bvec_from_jvec(const boot_init_t &iboot_ind,const vmeas_t<jack_t<T>> &jvec)
{
  vmeas_t<boot_t<T>> out(jvec.size());
  for(size_t it=0;it<jvec.size();it++) out[it].fill_from_jack(iboot_ind,jvec[it]);
  return out;
}

//! read from a set of confs
djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,size_t nlines=1);

#endif
