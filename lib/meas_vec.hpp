#ifndef _MEAS_VEC_HPP
#define _MEAS_VEC_HPP

#include <boot.hpp>

//! type defining a vector of measures
template <class meas_t> class vmeas_t : public vector<meas_t>
{
  //! bind the base type of meas_t
  using base_type=typename meas_t::base_type;
  
public:
  //! constructor specifying nel only (avoid copy constructor)
  explicit vmeas_t(size_t nel=0) : vector<meas_t>(nel) {}
  
  //! init from vector of vector
  vmeas_t(const vector<vector<base_type>> &o) : vmeas_t(o.size())
  {for(size_t it=0;it<o.size();it++) (*this)[it]=o[it];}
  
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
  {
    out.bin_write<size_t>(njacks);
    out.bin_write(this->size());
    out.bin_write(*this);
  }
  
  //! wrapper with name
  void bin_write(const char *path)
  {bin_write(raw_file_t(path,"w"));}
  
  //! wrapper with name
  void bin_write(const string &path)
  {bin_write(path.c_str());}
  
  //! assign from a scalar
  vmeas_t& operator=(const meas_t &oth) {for(auto &it : *this) it=oth;return *this;}
};

//! typically we use double
using djvec_t=vmeas_t<jack_t<double>>;
using dbvec_t=vmeas_t<boot_t<double>>;

//! initialize from a jvec
template <class T> vmeas_t<boot_t<T>> bvec_from_jvec(const boot_init_t &iboot_ind,const vmeas_t<jack_t<T>> &jvec)
{
  vmeas_t<boot_t<T>> out(jvec.size());
  for(size_t it=0;it<jvec.size();it++) out[it].fill_from_jack(iboot_ind,jvec[it]);
  return out;
}

//! read from a set of confs
djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,int nlines);

#endif
