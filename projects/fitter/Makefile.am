AM_CPPFLAGS+=-I$(top_srcdir)/projects/%D%

bin_PROGRAMS+= \
	$(top_builddir)/bin/fitter

%C%_fitter_parser_sources= \
         %D%/fitter_parser.ypp \
         %D%/fitter_tokenizer.lpp
%C%_fitter_parser_built_sources= \
         %D%/fitter_tokenizer.cpp \
         %D%/fitter_parser.hpp \
         %D%/fitter_parser.cpp \
	 %D%/fitter_tokenizer.cpp

%C%_fitter_tokenizer.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
BUILT_SOURCES+=$(%C%_fitter_parser_built_sources)
CLEANFILES+=$(%C%_fitter_parser_built_sources)


__top_builddir__bin_fitter_SOURCES= \
	$(%C%_fitter_parser_sources) \
	%D%/fitter.cpp \
	%D%/operations.cpp
