ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include

LDADD=$(top_builddir)/lib/libtran.a

bin_PROGRAMS= \
        $(top_builddir)/bin/B_to_Pi \
        $(top_builddir)/bin/combo_opt \
        $(top_builddir)/bin/exc_cont \
        $(top_builddir)/bin/smerad

__top_builddir__bin_exc_cont_SOURCES= \
	exc_cont.cpp

__top_builddir__bin_B_to_Pi_SOURCES= \
	B_to_Pi.cpp

__top_builddir__bin_smerad_SOURCES= \
	smerad.cpp

__top_builddir__bin_combo_opt_SOURCES= \
	combo_opt.cpp

