#ifndef _AVE_ERR_HPP
#define _AVE_ERR_HPP

#include <cmath>
#include <fstream>
#include <math.hpp>
#include <valarray>
#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////// average and error /////////////////////////////////////////////////

//! average and error
class ave_err_t : pair<double,double>
{
public:
  //! rebind base constructor
  ave_err_t(double a=0,double b=0) : pair<double,double>(a,b) {};
  
  //! move constructor
  ave_err_t(ave_err_t&& oth)=default;
  
  //! copy constructor
  ave_err_t(const ave_err_t &oth)=default;
  
  //! move assignement
  ave_err_t &operator=(ave_err_t &&)=default;
  
  //! copy assignement
  ave_err_t &operator=(const ave_err_t &oth)
  {pair<double,double>::operator=(oth);return *this;}
  
  //! rebind average
  double &ave=first;
  
  //! rebind error
  double &err=second;
  
  //! return average-error
  double ave_minus_err(){return ave-err;}
  
  //! return average+error
  double ave_plus_err(){return ave+err;}
  
  //! check that is printable
  bool is_printable() const {return isfinite(ave) and isfinite(err);}
};

//! compute average and stddev of a vector
template <class T> ave_err_t range_ave_stddev(const valarray<T> &v,size_t size)
{
  ave_err_t ae;
  
  for(size_t i=0;i<size;i++)
    {
      double x=v[i];
      ae.ave+=x;
      ae.err+=sqr(x);
    }
    ae.ave/=size;
    ae.err/=size;
    ae.err-=sqr(ae.ave);
    ae.err=sqrt(fabs(ae.err));
    
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

//! output of ave_err_t
ostream& operator<<(ostream &out,const ave_err_t &ae);

#endif
