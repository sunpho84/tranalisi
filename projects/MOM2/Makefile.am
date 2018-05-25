AM_CPPFLAGS+=-I$(top_srcdir)/projects/%D%

bin_PROGRAMS+= \
	$(top_builddir)/bin/mom

%C%_parser_sources= \
        %D%/parser.ypp \
        %D%/tokenizer.lpp \
        %D%/driver.hpp
%C%_parser_built_sources= \
        %D%/tokenizer.cpp \
        %D%/parser.hpp \
        %D%/parser.cpp \
	%D%/tokenizer.cpp

%C%_tokenizer.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
BUILT_SOURCES+=$(%C%_parser_built_sources)
CLEANFILES+=$(%C%_parser_built_sources)


__top_builddir__bin_mom_SOURCES= \
	$(%C%_parser_sources) \
	%D%/analysis.cpp \
	%D%/corrections.cpp \
	%D%/contractions.cpp \
	%D%/deltam.cpp \
	%D%/geometry.cpp \
	%D%/meson_mass.cpp \
	%D%/mom.cpp \
	%D%/pars.cpp \
	%D%/perens.cpp \
	%D%/pr_bil.cpp \
	%D%/pr_meslep.cpp \
	%D%/prop.cpp \
	%D%/scan.cpp \
	%D%/sigma.cpp \
	%D%/timings.cpp \
	%D%/Zbil.cpp \
	%D%/Zmeslep.cpp \
	%D%/Zq.cpp
