#include <eigen3/Eigen/Dense>

#include <mpfr.h>

/// Undef this to disable actual usage of the high preceision
//#define FAKE_HP

/// Structure to represent arbitray precision real number
struct PrecFloat
{
  /// Sets the default precision
  static void setDefaultPrecision(const int& n)
  {
#ifndef FAKE_HP
#pragma omp parallel
    mpfr_set_default_prec(n);
#endif
  }
  
  /// Storage
#ifndef FAKE_HP
  mpfr_t data{};
#else
  double data;
#endif
  
  /// Returns the internal data
  double get() const
  {
#ifdef FAKE_HP
    return data;
#else
    return
      mpfr_get_d(data,MPFR_RNDD);
#endif
  }
  
  /// Assignment
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
  
  /// Assign from another number
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
  
  /// Default initialization
  PrecFloat()
  {
#ifdef FAKE_HP
#else
    mpfr_init(data);
#endif
  }
  
  /// Copy constructor
  PrecFloat(const PrecFloat& oth)
  {
#ifdef FAKE_HP
    data=oth.data;
#else
    mpfr_init_set(data,oth.data,MPFR_RNDD);
#endif
  }
  
  /////////////////////////////////////////////////////////////////  
  
#ifdef FAKE_HP
#define PROVIDE_CONVERSION_FROM(TYPE,)			\
  PrecFloat(const TYPE& in)				\
  {							\
    data=in;						\
  }
#else
#define PROVIDE_CONVERSION_FROM(TYPE,MPFR_TAG)		\
  PrecFloat(const TYPE& in)				\
  {							\
  mpfr_init_set_ ## MPFR_TAG(data,in,MPFR_RNDD);	\
  }
#endif
  
  PROVIDE_CONVERSION_FROM(double,d);
  PROVIDE_CONVERSION_FROM(int,si);
  PROVIDE_CONVERSION_FROM(unsigned long int,ui);
  
#undef PROVIDE_CONVERSION_FROM
  
  /// Destructor
  ~PrecFloat()
  {
#ifdef FAKE_HP
#else
    mpfr_clear(data);
#endif
  }
  
  /////////////////////////////////////////////////////////////////
  
#ifdef FAKE_HP
#define BINARY_COMPARISON_OPERATOR_HELPER(NAME,MPFR_NAME)	\
  data NAME in.data
#else
#define BINARY_COMPARISON_OPERATOR_HELPER(NAME,MPFR_NAME)	\
  MPFR_NAME(data,in.data)
#endif
  
#define PROVIDE_BINARY_COMPARISON_OPERATOR(NAME,MPFR_NAME)	\
  								\
  inline bool operator NAME(const PrecFloat& in) const	\
  {								\
    return							\
      BINARY_COMPARISON_OPERATOR_HELPER(NAME,MPFR_NAME);	\
  }

  PROVIDE_BINARY_COMPARISON_OPERATOR(<,mpfr_less_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(<=,mpfr_lessequal_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(>,mpfr_greater_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(>=,mpfr_greaterequal_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(==,mpfr_equal_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(!=,!mpfr_equal_p);
  
#undef BINARY_COMPARISON_OPERATOR_HELPER
#undef BINARY_COMPARISON_OPERATOR
  
  /////////////////////////////////////////////////////////////////
  
  // Providing the binary operator as nonmember function allows to
  // take into account automatically also the cases in which the first
  // operand is not a PrecFloat
  
#ifdef FAKE_HP
#define BINARY_OPERATOR_HELPER(NAME,MPFR_NAME)\
  out.data=in1.data NAME in2.data
#else
#define BINARY_OPERATOR_HELPER(NAME,MPFR_NAME)\
  MPFR_NAME(out.data,in1.data,in2.data,MPFR_RNDD)
#endif
  
#define PROVIDE_SELF_BINARY_OPERATOR(NAME)			\
  								\
  inline PrecFloat& operator NAME ## =(const PrecFloat& in)	\
    {								\
      return							\
      (*this)=(*this)NAME in;					\
    }
  
#define PROVIDE_BINARY_OPERATOR(NAME,MPFR_NAME)			\
  								\
  friend inline PrecFloat operator NAME(const PrecFloat& in1,	\
					const PrecFloat& in2)	\
  {								\
    PrecFloat out;						\
  								\
    BINARY_OPERATOR_HELPER(NAME,MPFR_NAME);			\
								\
    return out;							\
  }								\
  								\
  PROVIDE_SELF_BINARY_OPERATOR(NAME);				\
  
  PROVIDE_BINARY_OPERATOR(+,mpfr_add)
  PROVIDE_BINARY_OPERATOR(-,mpfr_sub)
  PROVIDE_BINARY_OPERATOR(*,mpfr_mul)
  PROVIDE_BINARY_OPERATOR(/,mpfr_div)
  
#undef BINARY_OPERATOR_HELPER
#undef PROVIDE_BINARY_OPERATOR
#undef PROVIDE_SELF_BINARY_OPERATOR
  
  /// Negation
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
};

/////////////////////////////////////////////////////////////////

#ifdef FAKE_HP
#define UNARY_HELPER(NAME,MPFR_NAME)\
  out.data=NAME(in.data)
#else
#define UNARY_HELPER(NAME,MPFR_NAME)\
  MPFR_NAME(out.data,in.data,MPFR_RNDD)
#endif

#define PROVIDE_UNARY_FUNCTION(NAME,MPFR_NAME)	\
						\
  inline PrecFloat NAME(const PrecFloat& in)	\
{						\
  PrecFloat out;				\
						\
  UNARY_HELPER(NAME,MPFR_NAME);			\
						\
  return					\
    out;					\
}

PROVIDE_UNARY_FUNCTION(exp,mpfr_exp)
PROVIDE_UNARY_FUNCTION(abs,mpfr_abs)
PROVIDE_UNARY_FUNCTION(sqrt,mpfr_sqrt)
PROVIDE_UNARY_FUNCTION(asin,mpfr_asin)
PROVIDE_UNARY_FUNCTION(acos,mpfr_acos)
PROVIDE_UNARY_FUNCTION(sin,mpfr_sin)
PROVIDE_UNARY_FUNCTION(cos,mpfr_cos)
PROVIDE_UNARY_FUNCTION(sinh,mpfr_sinh)
PROVIDE_UNARY_FUNCTION(cosh,mpfr_cosh)
PROVIDE_UNARY_FUNCTION(erf,mpfr_erf)
PROVIDE_UNARY_FUNCTION(erfc,mpfr_erfc)

#undef PROVIDE_UNARY_FUNCTION
#undef UNARY_HELPER

/////////////////////////////////////////////////////////////////

/// Precise definition of Pi
inline PrecFloat precPi()
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=M_PI;
#else
  mpfr_const_pi(out.data,MPFR_RNDD);
#endif
  
  return
    out;
}

/////////////////////////////////////////////////////////////////

// Tell Eigen how to deal with PrecFloat numbers

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

using PrecVect=
  Eigen::Matrix<PrecFloat,Eigen::Dynamic,1>;

using PrecMatr=
  Eigen::Matrix<PrecFloat,Eigen::Dynamic,Eigen::Dynamic>;
