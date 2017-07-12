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
};

//! compute average and stddev of a vector
template <class T> ave_err_t range_ave_stddev(const valarray<T> &v,size_t size)
{
  ave_err_t ae;
  
  for(size_t i=0;i<size;i++)
    {
      double x=v[i];
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
  explicit vec_ave_err_t(size_t in) : vector<ave_err_t>(in) {}
  
  //! write to a grace file
  void write(const string &path) const;
};


//! traits for jack_t
template <> class vector_traits<vec_ave_err_t> : public true_vector_traits<vec_ave_err_t> {};

//! output of ave_err_t
ostream& operator<<(ostream &out,const ave_err_t &ae);

//! print the average and errors considering all rounding
string smart_print(double ave,const vector<double> &errors,int ndigits=2);

inline string smart_print(double ave,double error,int ndigits=2)
{return smart_print(ave,vector<double>({error}),ndigits);}

inline string smart_print(const ave_err_t &ae,int ndigits=2)
{return smart_print(ae.ave(),vector<double>({ae.err()}),ndigits);}

template <class T,class=enable_if_t<has_method_ave_err<T>::value>> string smart_print(const T &o,int ndigits=2)
{return smart_print(o.ave_err(),ndigits);}

#endif
