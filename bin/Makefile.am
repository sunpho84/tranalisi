ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include

LDADD=$(top_builddir)/lib/libtran.a

bin_PROGRAMS= \
	djvec_load \
	effective_mass \
	main \
	prepare_jackknives

djvec_load_SOURCES=djvec_load.cpp

effective_mass_SOURCES=effective_mass.cpp

main_SOURCES=main.cpp

prepare_jackknives_SOURCES=prepare_jackknives.cpp
