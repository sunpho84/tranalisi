AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS=

AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include
LDADD=$(top_builddir)/lib/libtran.a

include FPIEM/Makefile.am
include IB/Makefile.am
include MOM
include SMETEST/Makefile.am
include SPIPOLLO/Makefile.am
include STEP_SCALING/Makefile.am
