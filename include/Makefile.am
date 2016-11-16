git_hash.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index Makefile
	@echo "Generating $@"
	@echo '#ifndef GIT_HASH' > $@
	@echo -n ' #define GIT_HASH "' >> $@
	@git rev-parse HEAD|tr -d "\n" >> $@
	@echo '"' >> $@
	@echo -n ' #define GIT_LOG "' >> $@
	@git log -1 --pretty=%B|tr -d "\n" >> $@
	@echo '"' >> $@
	@echo '#endif' >> $@

include_HEADERS= \
        git_hash.hpp
