#ifndef _AVE_ERR_HPP
#define _AVE_ERR_HPP

#include <cmath>
#include <fstream>
#include <math.hpp>
#include <oper.hpp>
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
  
  //! check that is printable
  bool is_printable() const {return isfinite(ave) and isfinite(err);}
};

//! compute average and stddev of a vector
template <class T> ave_err_t range_ave_stddev(const vector<T> &v,size_t size)
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
    ae.err=sqrt(ae.err);
    
    return ae;
}

using vec_ave_err_t=vector<ave_err_t>;

//! output of ave_err_t
ostream& operator<<(ostream &out,const ave_err_t &ae);

#endif
