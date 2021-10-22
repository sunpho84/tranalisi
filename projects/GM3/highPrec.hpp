#include <tranalisi.hpp>

#include <mpfr.h>

//#define FAKE_HP

struct PrecFloat
{
  static void setDefaultPrecision(const int& n)
  {
#ifndef FAKE_HP
    mpfr_set_default_prec(n);
#endif
  }
  
#ifndef FAKE_HP
  mpfr_t data{};
#else
  double data;
#endif
  
  double get() const
  {
#ifdef FAKE_HP
    return data;
#else
    return
      mpfr_get_d(data,MPFR_RNDD);
#endif
  }
  
  PrecFloat& operator=(const double& in)
  {
#ifdef FAKE_HP
    data=in;
#else
    mpfr_set_d(data,in,MPFR_RNDD);
#endif
    
    return
      *this;
  }
  
  PrecFloat& operator=(const PrecFloat& oth)
  {
#ifdef FAKE_HP
    data=oth.data;
#else
    mpfr_set(data,oth.data,MPFR_RNDD);
#endif
    
    return
      *this;
  }
  
  PrecFloat()
  {
#ifdef FAKE_HP
#else
    mpfr_init(data);
#endif
  }
  
  PrecFloat(const PrecFloat& oth)
  {
#ifdef FAKE_HP
    data=oth.data;
#else
    mpfr_init_set(data,oth.data,MPFR_RNDD);
#endif
  }
  
  PrecFloat(const double& in)
  {
#ifdef FAKE_HP
    data=in;
#else
    mpfr_init_set_d(data,in,MPFR_RNDD);
#endif
  }
  
  PrecFloat(const int& in)
  {
#ifdef FAKE_HP
    data=in;
#else
    mpfr_init_set_si(data,in,MPFR_RNDD);
#endif
  }
  
  PrecFloat(const unsigned long int& in)
  {
#ifdef FAKE_HP
    data=in;
#else
    mpfr_init_set_ui(data,in,MPFR_RNDD);
#endif
  }
  
  ~PrecFloat()
  {
#ifdef FAKE_HP
#else
    mpfr_clear(data);
#endif
  }
  
  bool operator<(const PrecFloat& oth) const
  {
    return
#ifdef FAKE_HP
      data<oth.data;
#else
      mpfr_less_p(this->data,oth.data);
#endif
  }
  
  bool operator>(const PrecFloat& oth) const
  {
    return
#ifdef FAKE_HP
      data>oth.data;
#else
      mpfr_greater_p(this->data,oth.data);
#endif
  }
  
  bool operator==(const PrecFloat& oth) const
  {
    return
#ifdef FAKE_HP
      data==oth.data;
#else
      mpfr_cmp(this->data,oth.data)==0;
#endif
  }
  
  bool operator!=(const PrecFloat& oth) const
  {
    return
#ifdef FAKE_HP
      data!=oth.data;
#else
      not ((*this)==oth);
#endif
  }
  
  PrecFloat operator+(const PrecFloat& b) const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out=data+b.data;
#else
    mpfr_add(out.data,this->data,b.data,MPFR_RNDD);
#endif
    
    return
      out;
  }
  
  PrecFloat& operator+=(const PrecFloat& b)
  {
    (*this)=
      (*this)+b;
    
    return
     *this;
  }
  
  PrecFloat operator-(const PrecFloat& b) const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out.data=data-b.data;
#else
    mpfr_sub(out.data,this->data,b.data,MPFR_RNDD);
#endif
    
    return
      out;
  }
  
  PrecFloat& operator-=(const PrecFloat& b)
  {
    (*this)=
      (*this)-b;
    
    return
     *this;
  }
  
  PrecFloat operator-() const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out.data=-data;
#else
    mpfr_neg(out.data,this->data,MPFR_RNDD);
#endif
    
    return
      out;
  }
  
  PrecFloat operator*(const PrecFloat& b) const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out.data=data*b.data;
#else
    mpfr_mul(out.data,this->data,b.data,MPFR_RNDD);
#endif
    
    return
      out;
  }
  
  PrecFloat operator*(const double& b) const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out.data=data*b;
#else
    mpfr_mul_d(out.data,this->data,b,MPFR_RNDD);
#endif
    
    return
      out;
  }
  
  PrecFloat& operator*=(const PrecFloat& b)
  {
    (*this)=
      (*this)*b;
    
    return
      *this;
  }
  PrecFloat operator/(const PrecFloat& b) const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out.data=data/b.data;
#else
    mpfr_div(out.data,this->data,b.data,MPFR_RNDD);
#endif
    
    return
      out;
  }
  
  PrecFloat& operator/=(const PrecFloat& b)
  {
    (*this)=
      (*this)/b;
    
    return
      *this;
  }
};

inline PrecFloat operator*(const double& a,const PrecFloat& b)
{
  return
    b*a;
}

inline PrecFloat exp(const PrecFloat& in)
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=exp(in.data);
#else
  mpfr_exp(out.data,in.data,MPFR_RNDD);
#endif
  
  return
    out;
}

inline PrecFloat abs(const PrecFloat& in)
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=abs(in.data);
#else
  mpfr_abs(out.data,in.data,MPFR_RNDD);
#endif
  
  return
    out;
}

inline PrecFloat sqrt(const PrecFloat& in)
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=sqrt(in.data);
#else
  mpfr_sqrt(out.data,in.data,MPFR_RNDD);
#endif
  
  return
    out;
}

inline PrecFloat erf(const PrecFloat& in)
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=erf(in.data);
#else
  mpfr_erf(out.data,in.data,MPFR_RNDD);
#endif
  
  return
    out;
}

namespace Eigen
{
  template<>
  struct NumTraits<PrecFloat> :
    GenericNumTraits<PrecFloat>
  {
    typedef PrecFloat Real;
    typedef PrecFloat NonInteger;
    typedef PrecFloat Nested;
    
    static inline Real epsilon()
    {
      return 0;
    }
    
    static inline Real dummy_precision()
    {
      return 0;
    }
    
    static inline int digits10()
    {
      return 0;
    }
    
    enum
      {
	IsInteger=0,
	IsSigned=1,
	IsComplex=0,
	RequireInitialization=1,
	ReadCost=6,
	AddCost=150,
	MulCost=100
      };
  };
}
