git_hash.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/generate_hash.sh Makefile
	@echo "Generating $@"
	@$(top_srcdir)/config/generate_hash.sh > $@

include_HEADERS= \
        git_hash.hpp

include_HEADERS+= \
	%D%/config.hpp

include_HEADERS+= \
	%D%/tranalisi.hpp
