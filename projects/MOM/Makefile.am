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
	%D%/corrections.cpp \
	%D%/contractions.cpp \
	%D%/deltam_cr.cpp \
	%D%/geometry.cpp \
	%D%/ingredients.cpp \
	%D%/meson_masses.cpp \
	%D%/mom.cpp \
	%D%/prop.cpp \
	%D%/read.cpp \
	%D%/timings.cpp \
	%D%/Zbil.cpp \
	%D%/Zmeslep.cpp \
	%D%/Zq.cpp \
	%D%/Zq_sig1.cpp
