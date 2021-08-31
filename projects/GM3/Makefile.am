AM_CPPFLAGS+=-I$(top_srcdir)/projects/%D%

bin_PROGRAMS+= \
	$(top_builddir)/bin/gm3

__top_builddir__bin_gm3_SOURCES=%D%/gm3.cpp
__top_builddir__bin_gm3_CXXFLAGS=@MPI_CXXFLAGS@
__top_builddir__bin_gm3_LDFLAGS=@MPI_CXXLDFLAGS@
__top_builddir__bin_gm3_LDADD=@MPI_CXXLIBS@ $(top_builddir)/lib/libtran.a

