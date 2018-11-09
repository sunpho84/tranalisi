AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS=
BUILT_SOURCES=
CLEANFILES=

AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'`_ -olex.yy.c
AM_YFLAGS = -d -p `basename $* | sed 's,y$$,,'`_

AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include
LDADD=$(top_builddir)/lib/libtran.a

include FPIEM/Makefile.am
include IB/Makefile.am
include MOM2/Makefile.am
include PLNG/Makefile.am
include SMETEST/Makefile.am
include SPIPOLLO/Makefile.am
include STEP_SCALING/Makefile.am
