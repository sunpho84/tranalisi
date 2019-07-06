#ifndef _AVE_ERR_HPP
#define _AVE_ERR_HPP

#include <cmath>
#include <fstream>
#include <math.hpp>
#include <macros.hpp>
#include <tools.hpp>
#include <valarray>
#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////// average and error /////////////////////////////////////////////////

//! average and error
class ave_err_t : public pair<double,double>
{
public:
  //! rebind base constructor
  ave_err_t(const double &a=0,const double &b=0) :
    pair<double,double>(a,b) {};
  
  //! rebind average
  double ave() const {return first;}
  
  //! rebind error
  double err() const {return second;}
  
  //! rebind average
  double &ave() {return first;}
  
  //! rebind error
  double &err() {return second;}
  
  //! return average-error
  double ave_minus_err() const {return ave()-err();}
  
  //! return average+error
  double ave_plus_err() const {return ave()+err();}
  
  //! check that is printable
  bool is_printable() const {return isfinite(ave()) and isfinite(err());}
  
  //! return the significativity
  double significativity() const
  {
    return fabs(err()/ave());
  }
};

//! compute average and stddev of a vector
template <class V>
ave_err_t range_ave_stddev(const V &v,int size=-1)
{
  if(size==-1)
    size=v.size();
  
  ave_err_t ae;
  
  for(int i=0;i<size;i++)
    {
      auto x=v[i];
      ae.ave()+=x;
      ae.err()+=sqr(x);
    }
  ae.ave()/=size;
  ae.err()/=size;
  ae.err()-=sqr(ae.ave());
  ae.err()=sqrt(fabs(ae.err()));
  
  return ae;
}

//! a vector of ave_err_t
class vec_ave_err_t : public vector<ave_err_t>
{
public:
  vec_ave_err_t() {}
  
  explicit vec_ave_err_t(size_t in) : vector<ave_err_t>(in) {}
  
  //! write to a grace file
  void write(const vector<double> &x,const string &path) const;
  void write(const string &path)
  {write(vector_up_to<double>(this->size()),path);}
  
  vector<double> significativity() const
  {
    vector<double> out(this->size());
    transform(this->begin(), this->end(), out.begin(), [](const ave_err_t& x){return x.significativity();});
    
    return out;
  }
};

//! traits for jack_t
template <>
class vector_traits<vec_ave_err_t> : public true_vector_traits<vec_ave_err_t> {};

//! output of ave_err_t
ostream& operator<<(ostream &out,const ave_err_t &ae);

//! print the average and errors considering all rounding
string smart_print(double ave,const vector<double> &errors,int ndigits=2);

inline string smart_print(double ave,double error,int ndigits=2)
{return smart_print(ave,vector<double>({error}),ndigits);}

inline string smart_print(const ave_err_t &ae,int ndigits=2)
{return smart_print(ae.ave(),vector<double>({ae.err()}),ndigits);}

template <class T,class=enable_if_t<has_method_ave_err<T>::value>>
string smart_print(const T &o,int ndigits=2)
{return smart_print(o.ave_err(),ndigits);}

template <class V>
ave_err_t err_with_err(V x,int n=-1)
{
  if(n==-1)
    n=x.size();
  
  const double a=x.ave();
  
  //amplify the fluctuation
  for(int i=0;i<n;i++)
    x[i]=x[i]*n-a*(n-1);
  
  //compute sum and sum of the squares
  double s=0;
  double s2=0;
  for(int j=0;j<n;j++)
    {
      s+=x[j];
      s2+=sqr(x[j]);
    }
  
  double se=0,s2e=0;
  for(int i=0;i<n;i++)
    {
      const double si=(s-x[i])/(n-1);
      const double s2i=(s2-sqr(x[i]))/(n-1);
      
      const double e=sqrt(fabs(s2i-si*si)/(n-1));
      
      se+=e;
      s2e+=e*e;
    }
  
  se/=n;
  s2e/=n;
  
  s2e=sqrt(fabs(s2e-se*se)*(n-1));
  
  return {se,s2e};
}

template <class V>
ave_err_t skewness(V x,int n=-1)
{
  if(n==-1)
    n=x.size();
  
  const double a=x.ave();
  
  //amplify the fluctuation
  for(int i=0;i<n;i++)
    x[i]=x[i]*n-a*(n-1);
  
  double ssk=0,s2sk=0;
  for(int i=0;i<n;i++)
    {
      double s=0;
      for(int j=0;j<n;j++)
	if(i!=j)
	  s+=x[j];
      s/=(n-1);
      
      double s2=0,s3=0;
      for(int j=0;j<n;j++)
	if(i!=j)
	  {
	    s2+=sqr(x[j]-s);
	    s3+=cube(x[i]-s);
	  }
      s2/=(n-2);
      s3/=(n-1);
      
      const double sk=s3/pow(s2,1.5);
      
      ssk+=sk;
      s2sk+=sk*sk;
    }
  
  ssk/=n;
  s2sk/=n;
  
  s2sk=sqrt(fabs(s2sk-ssk*ssk)/(n-1));
  
  return {ssk,s2sk};
}

#endif
