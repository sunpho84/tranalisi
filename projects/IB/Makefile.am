AM_CPPFLAGS+=-I$(top_srcdir)/projects/%D%

bin_PROGRAMS+= \
	$(top_builddir)/bin/Kl2_IB \
	$(top_builddir)/bin/gm2_IB_sc \
	$(top_builddir)/bin/manette \
	$(top_builddir)/bin/manette2 \
	$(top_builddir)/bin/meson_mass_IB \
	$(top_builddir)/bin/meson_mass_IB_glb \
	$(top_builddir)/bin/reno \
	$(top_builddir)/bin/reno_old \
	$(top_builddir)/bin/zs_fr_zp

#	$(top_builddir)/bin/test_antonino

__top_builddir__bin_manette_SOURCES=%D%/manette.cpp
__top_builddir__bin_reno_SOURCES=%D%/reno.cpp
__top_builddir__bin_reno_old_SOURCES=%D%/reno_old.cpp
#__top_builddir__bin_test_antonino_SOURCES=%D%/test_antonino.cpp
__top_builddir__bin_Kl2_IB_SOURCES=%D%/Kl2_IB.cpp %D%/Kl2_IB_FSE.cpp %D%/common.cpp
__top_builddir__bin_manette2_SOURCES=%D%/manette2.cpp
__top_builddir__bin_zs_fr_zp_SOURCES=%D%/zs_fr_zp.cpp
__top_builddir__bin_meson_mass_IB_SOURCES=%D%/meson_mass_IB.cpp
__top_builddir__bin_meson_mass_IB_glb_SOURCES=%D%/common.cpp %D%/meson_mass_IB_glb.cpp %D%/meson_mass_IB_fit.cpp
__top_builddir__bin_gm2_IB_sc_SOURCES=%D%/common.cpp %D%/gm2_IB_sc.cpp %D%/gm2_IB_integrators.cpp
