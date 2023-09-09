AM_CPPFLAGS+=-I$(top_srcdir)/projects/%D%

###

bin_PROGRAMS+=$(top_builddir)/bin/trambacco

__top_builddir__bin_trambacco_SOURCES= \
	%D%/loadCorrCov.cpp \
	%D%/trambacco.cpp

###

bin_PROGRAMS+=$(top_builddir)/bin/tauVus

__top_builddir__bin_tauVus_SOURCES= \
	%D%/loadCorrCov.cpp \
	%D%/tauVus.cpp

# __top_builddir__bin_trambacco_LDADD=$(top_builddir)/lib/libtran.a -lmpfr
