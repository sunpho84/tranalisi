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

include_HEADERS+= \
	%D%/ave_err.hpp \
	%D%/boot.hpp \
	%D%/compl.hpp \
	%D%/Dirac.hpp \
	%D%/effective.hpp \
	%D%/file.hpp \
	%D%/fit.hpp \
	%D%/functions.hpp \
	%D%/gevp.hpp \
	%D%/grace.hpp \
	%D%/gslWrap.hpp \
	%D%/hashedFunction.hpp \
	%D%/index.hpp \
	%D%/jack.hpp \
	%D%/macros.hpp \
	%D%/minimizer.hpp \
	%D%/math.hpp \
	%D%/meas_vec.hpp \
	%D%/obs_file.hpp \
	%D%/oper.hpp \
	%D%/random.hpp \
	%D%/raw_file.hpp \
	%D%/solve.hpp \
	%D%/stopwatch.hpp \
	%D%/systematics.hpp \
	%D%/tasks.hpp \
	%D%/tools.hpp \
	%D%/traits.hpp \
	%D%/vecFilter.hpp
