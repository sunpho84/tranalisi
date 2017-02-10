#ifndef _MEAS_VEC_HPP
#define _MEAS_VEC_HPP

#include <boot.hpp>
#include <file.hpp>
#include <grace.hpp>
#include <tools.hpp>

//! type defining a vector of measures
template <class meas_t> class vmeas_t : public valarray<meas_t>
{
public:
  //! bind the base type of meas_t
  using base_type=meas_t;
  
  //! constructor specifying nel only (avoid copy constructor)
  explicit vmeas_t(size_t nel=0) : valarray<meas_t>(nel) {}
  
  //! constructor specifying nel and an element
  explicit vmeas_t(size_t nel,const meas_t &in) : valarray<meas_t>(in,nel) {}
  
  //! construct from jvec_t
  explicit vmeas_t(const boot_init_t &bi,const vmeas_t<jack_t<typename base_type::base_type>> &oth) : vmeas_t(oth.size())
  {for(size_t it=0;it<this->size();it++) (*this)[it]=boot_t<typename base_type::base_type>(bi,oth[it]);}
  
  //! init from vector of vector
  vmeas_t(const vector<vector<base_type>> &o) : vmeas_t(o.size())
  {
    auto start=take_time();
    cout<<"Initializing vmeas_t from vector of vector"<<endl;
    for(size_t it=0;it<o.size();it++) (*this)[it]=o[it];
    cout<<elapsed_time(start)<<" to init vmeas_t"<<endl;
  }
  
  //! move constructor
  vmeas_t(vmeas_t&& oth) : valarray<meas_t>(forward<valarray<meas_t>>(oth)) {// cout<<"vec move const"<<endl;
  }
  
  //! copy constructor
  vmeas_t(const vmeas_t &oth) : valarray<meas_t>(oth) {// cout<<"vec copy const"<<endl;
  }
  
  //! construct from sliced array
  vmeas_t(slice_array<meas_t> &&slice) : valarray<meas_t>(forward<valarray<meas_t>>(slice)) {}
  
  //! construct from sliced array
  vmeas_t(gslice_array<meas_t> &&gslice) : valarray<meas_t>(forward<valarray<meas_t>>(gslice)) {}
  
  //! construct from expr
  template<class _Dom> vmeas_t(const _Expr<_Dom,meas_t> &oth) : valarray<meas_t>(oth) {}
  
  //! move assignement
  vmeas_t &operator=(vmeas_t &&oth) noexcept =default;
  
  //! expression assignement
  //template<class op> vmeas_t &operator=(const _Expr<op,meas_t> &oth) {*this=oth;}
  
  //! copy assignement
  vmeas_t &operator=(const vmeas_t &oth)// =default;
  {valarray<meas_t>::operator=(oth);// cout<<"vec copy"<<endl;
    return *this;}
  
  //! assign from a scalar
  vmeas_t& operator=(const typename base_type::base_type &oth)
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
  
  //! clusterize each element
  void clusterize(size_t clust_size=1)
  {for(auto &a : *this) a.clusterize(clust_size);}
  
  //! assign from a scalar
  vmeas_t& operator=(const meas_t &oth) {for(auto &it : *this) it=oth;return *this;}
  
  //! return a subset including end
  vmeas_t subset(size_t beg,size_t end)
  {
    if(beg>end||beg>this->size()||beg>this->size()) CRASH("Asked to extract from %zu to %zu a vector of length %zu",beg,end,this->size());
    return (*this)[slice(beg,end-beg+1,1)];
  }
  
  //! return the tail-backked
  vmeas_t inverse()
  {
    size_t s=this->size();
    vmeas_t out(s);
    for(size_t it=0;it<s;it++) out[s-1-it]=(*this)[it];
    
    return out;
  }
  
  //! return the symmetric
  vmeas_t symmetric()
  {
    size_t s=this->size();
    vmeas_t out(s);
    for(size_t it=0;it<s;it++) out[(s-it)%s]=(*this)[it];
    
    return out;
  }
  
  //! return the averaged
  vmeas_t symmetrized(int par=1)
  {
    size_t nel=this->size();
    size_t nelh=nel/2;
    
    if(abs(par)!=1 && par!=0) CRASH("Unknown parity %d",par);
    
    if(nel%2) CRASH("Size %zu odd",nel);
    
    return vmeas_t(this->subset(0,nelh)+par*this->symmetric().subset(0,nelh))/(1.0+abs(par));//*(this->subset(0,nelh)+this->subset(nelh,nel).inverse());
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
  
  file.set_pos(ind*sizeof(typename T::base_type::base_type)*out[0].size()*nel);
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
