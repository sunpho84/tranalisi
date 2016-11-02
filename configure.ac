AC_INIT([tranalisi],[1.0],[fr.sanfilippo@gmail.com])

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

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)
AX_CXXFLAGS_WARN_ALL

#default number of bootstrap
AC_ARG_ENABLE(def_nboots,
        AS_HELP_STRING([--enable-def-nboots=n],[Set default number bootstap sample [Default: 100]]),
        de_nboots="${enableval}",
        def_nboots=100)
if ! printf "%f" "$def_nboots" >/dev/null 2>&1; then AC_MSG_ERROR([String $def_nboots is not a number]); fi
AC_MSG_RESULT([enabling number of default bootstrap samples... $def_nboots])
AC_DEFINE_UNQUOTED([DEF_NBOOTS],$def_nboots,[Default number of bootstrap])

AC_CONFIG_FILES(Makefile bin/Makefile lib/Makefile)
AC_OUTPUT
