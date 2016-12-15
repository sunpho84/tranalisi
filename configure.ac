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

AC_PROG_RANLIB
AC_PROG_CXX
AC_PROG_CXXCPP
AC_LANG([C++])

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
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

#search Minuit2
AC_CHECK_HEADERS([Minuit2/MinimumBuilder.h],[],AC_MSG_ERROR(["Unable to find headers of Minuit2! Install it using the script config/install_minuit2.sh"]))
LIBS="$LIBS -lMinuit2"
AC_MSG_CHECKING([for libMinuit2])
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([[#include <Minuit2/FCNBase.h>
		]],
		   [[class test : public ROOT::Minuit2::FCNBase {};]])],
  	       [AC_MSG_RESULT([yes])],
	       [AC_MSG_RESULT([no])
	        AC_MSG_ERROR([libMinuit2 is not installed.])])

#activate debug or not bootstrap
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

AC_CONFIG_FILES(Makefile
	include/Makefile
	lib/Makefile
	projects/Makefile
	projects/IB/Makefile
	tools/Makefile)
AC_OUTPUT
