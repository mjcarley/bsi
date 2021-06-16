## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libbsi.la

bin_PROGRAMS = \
	bsi-eval \
	bsi-init

bsi_eval_LDADD = \
	$(top_builddir)/src/.libs/libbsi.a 

bsi_init_LDADD = \
	$(top_builddir)/src/.libs/libbsi.a 
