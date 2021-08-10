AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"BSI\"

lib_LTLIBRARIES = libbsi.la

libbsi_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libbsi_la_SOURCES = \
	bsi.c \
	quadrature.c \
	quadrature.h \
	sources.c \
	rbf.c \
	rbf.h \
	interpolate.c \
	integrate.c \
	tetquad.c \
	tetquad.h

include_HEADERS = \
	bsi.h

CLEANFILES = $(BUILT_SOURCES)

