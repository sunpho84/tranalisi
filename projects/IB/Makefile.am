ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS=-I$(top_srcdir)/lib -I$(top_srcdir)/include

LDADD=$(top_builddir)/lib/libtran.a

bin_PROGRAMS= \
	$(top_builddir)/bin/gm2_IB_sc \
	$(top_builddir)/bin/Kl2_IB \
	$(top_builddir)/bin/manette \
	$(top_builddir)/bin/manette2 \
	$(top_builddir)/bin/meson_mass_IB \
	$(top_builddir)/bin/meson_mass_IB_glb \
	$(top_builddir)/bin/reno \
	$(top_builddir)/bin/reno_old

#	$(top_builddir)/bin/test_antonino

__top_builddir__bin_manette_SOURCES=manette.cpp
__top_builddir__bin_reno_SOURCES=reno.cpp
__top_builddir__bin_reno_old_SOURCES=reno_old.cpp
#__top_builddir__bin_test_antonino_SOURCES=test_antonino.cpp
__top_builddir__bin_Kl2_IB_SOURCES=Kl2_IB.cpp Kl2_IB_FSE.cpp common.cpp
__top_builddir__bin_manette2_SOURCES=manette2.cpp
__top_builddir__bin_meson_mass_IB_SOURCES=meson_mass_IB.cpp
__top_builddir__bin_meson_mass_IB_glb_SOURCES=common.cpp meson_mass_IB_glb.cpp meson_mass_IB_fit.cpp
__top_builddir__bin_gm2_IB_sc_SOURCES=common.cpp gm2_IB_sc.cpp gm2_IB_integrators.cpp
