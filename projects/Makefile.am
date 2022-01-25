AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS=
BUILT_SOURCES=
CLEANFILES=

AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'`_ -olex.yy.c
AM_YFLAGS = -Wno-yacc -d -p `basename $* | sed 's,y$$,,'`_

AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include
LDADD=$(top_builddir)/lib/libtran.a

if HAVE_HIGH_PREC
include SPECDENS/Makefile.am
endif
include GM3/Makefile.am
include fitter/Makefile.am
include VIRTUALPH/Makefile.am
include REPH/Makefile.am
include IB/Makefile.am
include HADRON_SCHEME/Makefile.am
include GILBACCO/Makefile.am
include BARYONS/Makefile.am
#include FPIEM/Makefile.am
#include PLNG/Makefile.am
#include MOM2/Makefile.am
include GM2/Makefile.am
#include SMETEST/Makefile.am
#include SPIPOLLO/Makefile.am
#include STEP_SCALING/Makefile.am
