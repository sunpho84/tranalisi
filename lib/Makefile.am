ACLOCAL_AMFLAGS=-I m4

%D%/git_info.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/get_git_info.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/get_git_info.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:


BUILT_SOURCES=%D%/git_info.hpp
lib_LIBRARIES=libtran.a

libtran_a_SOURCES= \
	ave_err.cpp \
	boot.cpp \
	Dirac.cpp \
	effective.cpp \
	file.cpp \
	fit.cpp \
	functions.cpp \
	gevp.cpp \
	grace.cpp \
	jack.cpp \
	meas_vec.cpp \
	minimizer.cpp \
	oper.cpp \
	systematics.cpp \
	tools.cpp

include_HEADERS= \
	%D%/git_info.hpp
