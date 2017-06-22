ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include

LDADD=$(top_builddir)/lib/libtran.a

bin_PROGRAMS= \
	combine_vdoubles \
	djvec_load \
	effective_mass \
	onlinemeas \
	prepare_jackknives \
	simple_effective_mass \
	two_pts_fit

combine_vdoubles_SOURCES=combine_vdoubles.cpp

djvec_load_SOURCES=djvec_load.cpp

effective_mass_SOURCES=effective_mass.cpp

onlinemeas_SOURCES=onlinemeas.cpp

prepare_jackknives_SOURCES=prepare_jackknives.cpp

simple_effective_mass_SOURCES=simple_effective_mass.cpp

two_pts_fit_SOURCES=two_pts_fit.cpp
