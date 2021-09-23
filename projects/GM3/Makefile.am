AM_CPPFLAGS+=-I$(top_srcdir)/projects/%D%

bin_PROGRAMS+= \
	$(top_builddir)/bin/gm3

#	$(top_builddir)/bin/gm3Test

#__top_builddir__bin_gm3Test_SOURCES= 
#	%D%/gm3Test.cpp

__top_builddir__bin_gm3_SOURCES= \
	%D%/amu.cpp \
	%D%/cachedData.cpp \
	%D%/data.cpp \
	%D%/gm3.cpp \
	%D%/mPCAC.cpp \
	%D%/kernels.cpp \
	%D%/params.cpp \
	%D%/rawData.cpp \
	%D%/renoConstants.cpp \
	%D%/VKVKRepresentation.cpp

__top_builddir__bin_gm3_CXXFLAGS=@MPI_CXXFLAGS@
__top_builddir__bin_gm3_LDFLAGS=@MPI_CXXLDFLAGS@
__top_builddir__bin_gm3_LDADD=@MPI_CXXLIBS@ $(top_builddir)/lib/libtran.a

