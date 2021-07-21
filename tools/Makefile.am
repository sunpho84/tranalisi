ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include

LDADD=$(top_builddir)/lib/libtran.a

bin_PROGRAMS= \
	backstab \
	combine_vdoubles \
	gevp \
	djvec_load \
	djvec_ratio \
	effective_mass \
	linfit \
	onlinemeas \
	planfit \
	prepare_jackknives \
	simple_effective_mass \
	test \
	test_hdf5 \
	two_pts_fit \
	w0

backstab_SOURCES=backstab.cpp

combine_vdoubles_SOURCES=combine_vdoubles.cpp

djvec_load_SOURCES=djvec_load.cpp

djvec_ratio_SOURCES=djvec_ratio.cpp

gevp_SOURCES=gevp.cpp

linfit_SOURCES=linfit.cpp

effective_mass_SOURCES=effective_mass.cpp

onlinemeas_SOURCES=onlinemeas.cpp

planfit_SOURCES=planfit.cpp

prepare_jackknives_SOURCES=prepare_jackknives.cpp

simple_effective_mass_SOURCES=simple_effective_mass.cpp

test_hdf5_SOURCES=test_hdf5.cpp
test_hdf5_CXXFLAGS=@MPI_CXXFLAGS@ @CXXFLAGS@
test_hdf5_LDFLAGS=@MPI_CXXLDFLAGS@ @LDFLAGS@ @MPI_CXXLIBS@ @LIBS@

test_SOURCES=test.cpp

two_pts_fit_SOURCES=two_pts_fit.cpp

w0_SOURCES=w0.cpp
