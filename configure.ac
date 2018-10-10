AC_INIT([tranalisi],[1.0],[fr.sanfilippo@gmail.com])

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIG_TIME=`date`
AC_DEFINE_UNQUOTED([CONFIG_TIME], "$CONFIG_TIME", [time when configured])
AC_MSG_RESULT([$CONFIG_TIME])

AC_CONFIG_SRCDIR([lib])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(include/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#silent automake
AM_SILENT_RULES([yes])

AC_CHECK_HEADERS([malloc.h])

#gsl
AC_ARG_WITH(gsl,
        AS_HELP_STRING([--with-gsl[=dir]], [Specify where to find gsl]),
        with_gsl="${withval}"
        CPPFLAGS="-I${with_gsl}/include/ $CPPFLAGS"
        LDFLAGS="-L${with_gsl}/lib/ $LDFLAGS")
AC_MSG_RESULT([with gsl... $with_gsl])

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,AC_MSG_ERROR([Unable to find gslcblas]))
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,AC_MSG_ERROR([Unable to find gsl]))

AC_PROG_RANLIB
AC_PROG_CXX
AC_PROG_CXXCPP
AC_LANG([C++])

AC_PROG_YACC
AM_PROG_LEX

AX_CXX_COMPILE_STDCXX(14,noext,mandatory)
AX_CXXFLAGS_WARN_ALL

#check openmp
AX_OPENMP(AC_DEFINE_UNQUOTED([USE_OMP],1,[Use or not OpenMP]))
CPPFLAGS="$CPPFLAGS $OPENMP_CPPFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

#default number of bootstrap
AC_ARG_ENABLE(def_nboots,
        AS_HELP_STRING([--enable-def-nboots=n],[Set default number bootstap sample [Default: 100]]),
        def_nboots="${enableval}",
        def_nboots=100)
if ! printf "%f" "$def_nboots" >/dev/null 2>&1; then AC_MSG_ERROR([String $def_nboots is not a number]); fi
AC_MSG_RESULT([enabling number of default bootstrap samples... $def_nboots])
AC_DEFINE_UNQUOTED([DEF_NBOOTS],$def_nboots,[Default number of bootstrap])

#decide whether to propagate or not the average
AC_DEFINE_UNQUOTED([PROP_MEAN],[0],["Propagate the mean"])
AC_DEFINE_UNQUOTED([DISTR_MEAN],[1],["Compute the mean from distribution"])
AC_ARG_ENABLE(mean,
	AS_HELP_STRING([--enable-mean],[Propagate the original average [prop] or (default) compute average from distribtuion [distr]]),
	enable_mean="${enableval}",enable_mean=distr)
case "$enable_mean" in
     prop) mean_type=PROP_MEAN;;
     distr) mean_type=DISTR_MEAN;;
     *) AC_MSG_ERROR(["Unkwnown mean propagation type ${enableval}"])
esac
AC_DEFINE_UNQUOTED([MEAN_TYPE],[$mean_type],[How to compute mean])
AC_MSG_RESULT([Enable mean type... ${enable_mean}])

#search Minuit
AC_ARG_WITH(Minuit,
        AS_HELP_STRING([--with-Minuit[=dir]], [Specify where to find Minuit]),
        with_Minuit="${withval}"
        CPPFLAGS="-I${with_Minuit}/include/ $CPPFLAGS"
        LDFLAGS="-L${with_Minuit}/lib/ $LDFLAGS",
        with_Minuit=no)
AC_CHECK_HEADERS([TMinuit.h],[],dnl AC_MSG_ERROR(["Unable to find headers of Minuit, good luck!"])
)
AC_MSG_CHECKING([for libMinuit])
SAVED_LIBS=$LIBS
LIBS="$LIBS -lMinuit -lCore"
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([[#include <TMinuit.h>
		]],
		   [[TMinuit minu;]])],
  	       [AC_MSG_RESULT([yes]);Minuit_found=true],
	       [AC_MSG_RESULT([no]);Minuit_found=false])
LIBS=$SAVED_LIBS

#search Minuit2
AC_ARG_WITH(Minuit2,
        AS_HELP_STRING([--with-Minuit2[=dir]], [Specify where to find Minuit2]),
        with_Minuit2="${withval}"
        LDFLAGS="-L${with_Minuit2}/lib/ $LDFLAGS",
        with_Minuit2=no)
AC_CHECK_HEADERS([Minuit2/MinimumBuilder.h],[],dnl AC_MSG_ERROR(["Unable to find headers of Minuit2! Install it using the script config/install_minuit2.sh"])
)
AC_MSG_CHECKING([for libMinuit2])
SAVED_LIBS=$LIBS
LIBS="$LIBS -lMinuit2"
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([[
		#include <Minuit2/FCNBase.h>
		#include <Minuit2/MnPrint.h>
		]],
		   [[class test_t : public ROOT::Minuit2::FCNBase
		   {
		   public:
			double operator()(const std::vector<double>&) const{return 0;}
			double Up() const {return 1;}
		   };
		test_t test;
		ROOT::Minuit2::MnPrint::SetLevel(0);
		]])],
  	       [AC_MSG_RESULT([yes]);Minuit2_found=true],
	       [AC_MSG_RESULT([no]);Minuit2_found=false])
LIBS=$SAVED_LIBS

#find the default minimizer
default_minimizer=none
#take Minuit if available
AC_MSG_RESULT([Minuit_found: ${Minuit_found}])
if test "${Minuit_found}" == true
then
	default_minimizer=Minuit
fi
#overwrite with Minuit2 if available
if test "${Minuit2_found}" == true
then
	default_minimizer=Minuit2
fi
#check that we have a valid default otherwise crash
if test "$default_minimizer" == none
then
	AC_MSG_ERROR(["No minimizer tool found"])
fi
AC_MSG_RESULT([Default minimizer: ${default_minimizer}])

#decide whether to propagate or not the average
AC_DEFINE_UNQUOTED([MINUIT],[0],["Use Minuit"])
AC_DEFINE_UNQUOTED([MINUIT2],[1],["Use Minuit2"])
AC_ARG_ENABLE(minimizer,
	AS_HELP_STRING([--enable-minimizer],[Choose the minimizer]),
	enable_minimizer="${enableval}",
	enable_minimizer="${default_minimizer}")
case "$enable_minimizer" in
     Minuit)
       if test $Minuit_found != true
       then
	   AC_MSG_ERROR([Cannot select not found minimizer Minuit])
       fi
       minimizer_type=MINUIT
       LIBS="$LIBS -lMinuit -lCore"
       ;;
     Minuit2)
       if test $Minuit2_found != true
       then
	   AC_MSG_ERROR([Cannot select not found minimizer Minuit2])
       fi
       minimizer_type=MINUIT2
       LIBS="$LIBS -lMinuit2"
       ;;
     *) AC_MSG_ERROR(["Unkwnown minimizer tool ${enableval}"])
esac
AC_DEFINE_UNQUOTED([MINIMIZER_TYPE],[$minimizer_type],[How to minimize])
AC_MSG_RESULT([Enable minimizer... ${enable_minimizer}])

#activate debug or not
AC_ARG_ENABLE(debug,
        AS_HELP_STRING([--enable-debug=yes[,no]],[Enable debug]),
        debug="${enableval}",
        debug=no)
case "$debug" in
     no) ;;
     yes) AC_DEFINE_UNQUOTED([DEBUG],$debug,[Enable debug mode]);;
     *) AC_MSG_ERROR(["Unkwnown debug method $debug"])
esac
AC_MSG_RESULT([Enabling debug... $debug])

#activate Eigen vectorization
AC_ARG_ENABLE(Eigen_vec_opt,
        AS_HELP_STRING([--enable-Eigen_vec_opt=yes[,no]],[Enable Eigen vectorization (default: no)]),
        Eigen_vec_opt="${enableval}",
        Eigen_vec_opt=no)
case "$Eigen_vec_opt" in
     no) AC_DEFINE_UNQUOTED([EIGEN_DONT_VECTORIZE],,[Disable Eigen_vectorization]);;
     yes) ;;
     *) AC_MSG_ERROR(["Unkwnown Eigen_vec_opt method $Eigen_vec_opt"])
esac
AC_MSG_RESULT([Enabling Eigen vectorization... $Eigen_vec_opt])

#check Eigen
AC_ARG_WITH(Eigen,
        AS_HELP_STRING([--with-Eigen[=dir]], [Specify where to find Eigen folder (include eigen3 in the argument!)]),
        with_Eigen="${withval}"
	CPPFLAGS="-I${with_Eigen} $CPPFLAGS")
AC_MSG_RESULT([with Eigen... $with_Eigen])
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([Eigen/Dense],,AC_MSG_ERROR([Error! Eigen required]))
AC_LANG_POP([C++])

#activate arch optimization
AC_ARG_ENABLE(arch-optimizations,
        AS_HELP_STRING([--enable-arch-optimizations=val[-march=native]],[Enable architecture optimizations]),
        arch_opt="${enableval}",
        arch_opt=-march=native)
case "$arch_opt" in
     no) ;;
     *) CXXFLAGS="${arch_opt} $CXXFLAGS";;
esac
AC_MSG_RESULT([Enabling architecture optimization... $arch_opt])

AC_CONFIG_FILES(Makefile
	include/Makefile
	lib/Makefile
	projects/Makefile
	tools/Makefile)
AC_OUTPUT
