#ifndef _MEAS_VEC_HPP
#define _MEAS_VEC_HPP

#include <boot.hpp>
#include <file.hpp>
#include <grace.hpp>
#include <obs_file.hpp>
#include <stopwatch.hpp>
#include <tools.hpp>

//! type defining a vector of measures
template <class meas_t> class vmeas_t : public valarray<meas_t>
{
public:
  //! safe resize
  void resize(size_t sz,meas_t C=meas_t())
  {
    if(sz!=this->size())
      {
	auto old=*this;
	this->valarray<meas_t>::resize(sz);
	for(size_t i=0;i<min(this->size(),old.size());i++) (*this)[i]=old[i];
      }
  }
  
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
  
  //! init (as for STL containers)
  meas_t* begin() {return &((*this)[0]);}
  const meas_t* begin() const {return &((*this)[0]);}
  
  //! end (as for STL containers)
  meas_t* end() {return &((*this)[0])+this->size();}
  const meas_t* end() const {return &((*this)[0])+this->size();}
  
  //! compute average and error
  vec_ave_err_t ave_err() const
  {
    vec_ave_err_t out(this->size());
    for(size_t it=0;it<this->size();it++) out[it]=(*this)[it].ave_err();
    return out;
  }
  
  //! compute average
  vector<double> ave() const
  {
    vector<double> out(this->size());
    for(size_t it=0;it<this->size();it++) out[it]=(*this)[it].ave();
    return out;
  }
  
  //! write to a stream
  void bin_write(raw_file_t &out) const
  {out.bin_write(*this);}
  
  //! wrapper with name
  void bin_write(const char *path) const
  {
    raw_file_t fout(path,"w");
    bin_write(fout);
  }
  
  //! wrapper with name
  void bin_write(const string &path) const
  {bin_write(path.c_str());}
  
  //! read from a stream
  void bin_read(const raw_file_t &in)
  {in.bin_read(*this);}
  
  //! wrapper with name
  void bin_read(const char *path)
  {
    raw_file_t fin(path,"r");
    bin_read(fin);
  }
  
  //! wrapper with name
  void bin_read(const string &path)
  {bin_read(path.c_str());}
  
  //! clusterize each element
  vmeas_t &clusterize(size_t clust_size=1)
  {
    for(auto &a : *this) a.clusterize(clust_size);
    return *this;
  }
  
  //! assign from a scalar
  vmeas_t& operator=(const meas_t &oth) {for(auto &it : *this) it=oth;return *this;}
  
  //! return a subset including end
  vmeas_t subset(size_t beg,size_t end) const
  {
    return (*this)[slice(std::max(size_t(0),beg),std::min(end-beg+1,this->size()),1)];
  }
  
  //! return a filtered vector
  template <class filter_fun_t>
  vmeas_t filter(const filter_fun_t &filter_fun) const
  {
    size_t n=0;
    for(size_t i=0;i<this->size();i++) n+=filter_fun(i);
    
    vmeas_t out(n);
    size_t j=0;
    for(size_t i=0;i<this->size();i++)
      {
	bool f=filter_fun(i);
	if(f) out[j++]=(*this)[i];
      }
    
    return out;
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
  vmeas_t symmetric() const
  {
    size_t s=this->size();
    vmeas_t out(s);
    for(size_t it=0;it<s;it++) out[(s-it)%s]=(*this)[it];
    
    return out;
  }
  
  //! return the averaged
  vmeas_t symmetrized(int par=1) const
  {
    size_t nel=this->size();
    size_t nelh=nel/2;
    
    if(abs(par)!=1 and par!=0) CRASH("Unknown parity %d",par);
    
    if(nel%2) CRASH("Size %zu odd",nel);
    
    return vmeas_t(this->subset(0,nelh)+par*this->symmetric().subset(0,nelh))/(1.0+abs(par));//*(this->subset(0,nelh)+this->subset(nelh,nel).inverse());
  }
  
  //! self-assign the symmetrized
  vmeas_t symmetrize(int par=1)
  {
    (*this)=symmetrized(par);
    return *this;
  }
  
  //! get a specific "slice" of events
  vector<double> get_all_events(size_t i) const
  {
    vector<double> out(this->size());
    transform(this->begin(),this->end(),out.begin(),[i](const base_type &b){return b[i];});
    return out;
  }
  
  //! put a specific "slice" of events
  void put_all_events(const vector<double> &in,size_t i)
  {for(size_t iel=0;iel<in.size();iel++) (*this)[iel][i]=in[iel];}
};

//! traits for vmeas_t
template <class TS> class vector_traits<vmeas_t<TS>> : public true_vector_traits<TS> {};

//! typically we use double
using djvec_t=vmeas_t<jack_t<double>>;
using dbvec_t=vmeas_t<boot_t<double>>;

//! read a binary from file
template <class TV> TV read_vec_meas(raw_file_t &file,size_t nel,size_t ind=0)
{
  TV out(nel);
  if(nel==0) CRASH("Nel==0");

  file.set_pos(ind*sizeof(base_type_t<base_type_t<TV>>)*out[0].size()*nel);
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

//! read from a list of confs
djvec_t read_conf_set_t(const string &template_path,vector<size_t> &id_list,size_t ntot_col,const vector<size_t> &cols,size_t nlines,bool verbosity=VERBOSE);

//! read from a range of confs
djvec_t read_conf_set_t(const string &template_path,const range_t &range,size_t ntot_col,const vector<size_t> &cols,size_t nlines=1,bool verbosity=VERBOSE);

//! integrate
template <class TS>
TS integrate_corr_up_to(const valarray<TS> &corr,size_t &upto,size_t ord=1)
{
  if(upto>=corr.size()) CRASH("Upto=%zu >= corr.size()=%zu",upto,corr.size());
  
  valarray<double> weight;
  switch(ord)
    {
    case 0: weight={1};       break;
    case 1: weight={1,1};     break;
    case 2: weight={1,4,1};   break;
    case 3: weight={3,9,9,3}; break;
    default: CRASH("Unknwown order %zu",ord);
    }
  
  size_t len=weight.size()-1;
  double norm=weight.sum()/len;
  
  TS out;
  out=0.0;
  upto=size_t(upto/len)*len;
  for(size_t t=0;t<upto;t+=len)
    for(size_t j=0;j<=len;j++)
      out+=corr[t+j]*weight[j];
  
  out/=norm;
  
  return out;
}

template <class TS>
TS integrate_corr(const valarray<TS> &corr,size_t &from,size_t &upto,size_t ord=1)
{
  if(upto>=corr.size()) CRASH("Upto=%zu >= corr.size()=%zu",upto,corr.size());
  
  valarray<double> weight;
  switch(ord)
    {
    case 0: weight={1};       break;
    case 1: weight={1,1};     break;
    case 2: weight={1,4,1};   break;
    case 3: weight={3,9,9,3}; break;
    default: CRASH("Unknwown order %zu",ord);
    }
  
  size_t len=weight.size()-1;
  double norm=weight.sum()/len;
  
  TS out;
  out=0.0;
  upto=size_t(upto/len)*len;
  for(size_t t=from;t<upto;t+=len)
    for(size_t j=0;j<=len;j++)
      out+=corr[t+j]*weight[j];
  
  out/=norm;
  
  return out;
}

//! determines the vector type, forward declaration
template <class T>
class meas_vec_of;

//! determine the vector type for a djack_t
template <>
class meas_vec_of<djack_t>
{
public:
  typedef djvec_t type;
};

//! determine the vector type for a dboot_t
template <>
class meas_vec_of<dboot_t>
{
public:
  typedef dbvec_t type;
};

//! helper class
template <class T>
using meas_vec_of_t=typename meas_vec_of<T>::type;

#endif
