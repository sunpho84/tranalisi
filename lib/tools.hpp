#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <Eigen/Dense>

#include <array>
#include <iostream>
#include <macros.hpp>
#include <string>
#include <valarray>
#include <vector>

#include <traits.hpp>

using namespace std;

DEFINE_HAS_METHOD(size);

DEFINE_HAS_METHOD(ave_err);

DEFINE_HAS_METHOD(is_printable);

//! crashes emitting the message
void internal_crash(int line,const char *file,const char *fun,const char *temp,...);

//! check if two quantities have the same sign
template <class T>
bool same_sign(const T &a,const T &b)
{return (a<=0 and b<=0) or (a>=0 and b>=0);}

//! check it two quantities have opposite sign
template <class T>
bool opposite_sign(const T &a,const T &b)
{return !same_sign(a,b);}

//! combine arguments in a single string
string combine(const char *format,...);

//! close gently with a message
void close_with_mess(const char *format,...);

//! flush the unused memory
void flush_unused_memory();

//!check if a file exists
int file_exists(string path);

//! check if a directoy exists
int dir_exists(string path);

//! handle signals
void signal_handler(int sig);

//! check that a vector is orderd
template <class T>
void check_ordered(const initializer_list<T> &vec)
{
  auto prev=vec.begin();
  auto cur=prev+1;
  size_t it=1;
  do
    {
      if(*prev>*cur)
	CRASH("Element %zu has value %s, greater than element %zu, %s",it-1,to_string(*prev).c_str(),it,to_string(*cur).c_str());
      prev=cur++;
      it++;
    }
  while(cur!=vec.end());
}

//! check agreement of sizes of two vectors
template <class T1,class T2,class=enable_if_t<is_vector<T1>::value and is_vector<T2>::value>>
void check_match_size(const T1 &first,const T2 &second)
{if(first.size()!=second.size()) CRASH("Vectors do not agree in size, %d vs %d",first.size(),second.size());}

//! return a range of int
class range_t : private array<size_t,3>
{
public:
  
  //! copy construct
  range_t(const range_t& oth)
  {
    static_cast<array<size_t,3>>(*this)=static_cast<array<size_t,3>>(oth);
  }
  
  const range_t& operator=(const range_t &oth)
  {
    static_cast<array<size_t,3>>(*this)=static_cast<array<size_t,3>>(oth);
    
    return oth;
  }
  
  //bind
  size_t &start=(*this)[0];
  size_t &each=(*this)[1];
  size_t &end=(*this)[2];
  
  //! get an element
  size_t operator()(size_t i) const
  {
    size_t id=i*each+start;
    if(id>end) CRASH("Asking an element %zu beyond end %zu",id,end);
    return id;
  }
  
  //! number of elements
  size_t size() const
  {return (end-start)/each+1;}
  
  //! default constructor
  range_t() {}
  
  //! init from triplet
  range_t(initializer_list<size_t> list)
  {
    if(list.size()!=3) CRASH("list size %d while expecting 3",list.size());
    copy(list.begin(),list.end(),this->begin());
  }
};

//! convert crashing if not working
inline int to_int(string s)
{
  int out=0;
  try{out=stoi(s);}
  catch(exception &e){CRASH("Exception %s while converting to int string \"%s\"",e.what(),s.c_str());}
  
  return out;
}

//! return a filled vector ranging from offset(0) to max (excluded) each stride(1)
template <class T>
vector<T> vector_up_to(const T max,const T offset=0,const T stride=1)
{
  vector<T> v;
  for(T x=offset;x<max;x+=stride) v.push_back(x);
  return v;
}

//! return a vector filled witha a grid ranging from xmin to xmax (included), with n points (n-1 intervals)
template <class T>
vector<T> vector_grid(const T xmin,const T xmax,const size_t n)
{
  if(n==0) CRASH("n must be different from 0");
  const T dx=(xmax-xmin)/(n-1);
  vector<T> y(n);
  for(size_t i=0;i<n;i++) y[i]=xmin+dx*i;
  return y;
}

//! set to zero a double
inline void set_to_zero(double &x)
{x=0;}

//! set to zero an eigen matr
template <typename D>
void set_to_zero(const Eigen::MatrixBase<D> &x)
{x.Zero();}

//! return a percentage
template <class T>
double percentage(const T &num,const T& den)
{return num*100.0/den;}

//! sum from n0 to infinity, with a given tolerance
template <class T,class Fun>
T series(int n0,int dn,const Fun &fun,T tol)
{
  T out=0;
  T contr;
  
  //positive contribution
  int n=n0,i=0;
  bool get_out;
  do
    {
      //compute, check nan and increment
      contr=fun(n);
      if(isnan(contr)) CRASH("nan found when evaluating %d in %s",n0,__PRETTY_FUNCTION__);
      out+=contr;
      
      //check if we go out
      if(fabs(out)==0.0) get_out=(i>10);
      else               get_out=(fabs(contr/out)<tol);
      
      //debug
      //cout<<"  "<<n<<" "<<out<<" "<<contr<<" "<<get_out<<endl;
      
      //increment
      n+=dn;
      i++;
    }
  while(not get_out);
  
  return out;
}

#endif
