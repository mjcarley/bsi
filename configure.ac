#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

##AC_PREREQ([2.69])
AC_INIT([bsi],[0.1.0],[codes@paraffinalia.co.uk])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/bsi.c])
AC_CONFIG_HEADERS([config.h])

BSI_MAJOR_VERSION=0
BSI_MINOR_VERSION=1
BSI_MICRO_VERSION=0
BSI_INTERFACE_AGE=0
BSI_BINARY_AGE=0
BSI_VERSION=$BSI_MAJOR_VERSION.$BSI_MINOR_VERSION.$BSI_MICRO_VERSION

AC_SUBST(BSI_MAJOR_VERSION)
AC_SUBST(BSI_MINOR_VERSION)
AC_SUBST(BSI_MICRO_VERSION)
AC_SUBST(BSI_VERSION)

# libtool versioning
LT_RELEASE=$BSI_MAJOR_VERSION.$BSI_MINOR_VERSION
LT_CURRENT=`expr $BSI_MICRO_VERSION - $BSI_INTERFACE_AGE`
LT_REVISION=$BSI_INTERFACE_AGE
LT_AGE=`expr $BSI_BINARY_AGE - $BSI_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

dnl Initialize libtool
LT_INIT

# For automake.
VERSION=$BSI_VERSION
PACKAGE=bsi

# Checks for programs.
AC_PROG_CC

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_SEARCH_LIBS([strerror],[cposix])
AC_HEADER_STDC

# Checks for libraries.
AM_PATH_GLIB_2_0(2.4.0, has_glib2=yes, has_glib2=no, gthread gmodule)
if test x$has_glib2 = xyes; then
  glib_cflags=`$PKG_CONFIG glib-2.0 --cflags`
  glib_thread_cflags=`$PKG_CONFIG glib-2.0 --cflags gthread-2.0`
  glib_module_cflags=`$PKG_CONFIG glib-2.0 --cflags gmodule-2.0`
  glib_libs=`$PKG_CONFIG glib-2.0 --libs`
  glib_thread_libs=`$PKG_CONFIG glib-2.0 --libs gthread-2.0`
  glib_module_libs=`$PKG_CONFIG glib-2.0 --libs gmodule-2.0`
  GLIB_LIBS="$glib_libs"
  GLIB_DEPLIBS="$glib_libs"
else
  if test x$with_glib = xyes ; then
    AC_MSG_ERROR([
*** Directory must be specified for --with-glib])
  fi

  if test x$with_glib = x ; then 
    # Look for separately installed glib

    GLIB_REQUIRED_VERSION=1.2.8
    AM_PATH_GLIB($GLIB_REQUIRED_VERSION,,
      AC_MSG_ERROR([
*** GLIB $GLIB_REQUIRED_VERSION or better is required. The latest version of GLIB
*** is always available from ftp://ftp.gtk.org/.]),
      gmodule gthread)

    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`
    GLIB_LIBS="$glib_libs"
    GLIB_DEPLIBS="$glib_libs"
  else
    # Use uninstalled glib (assume they got the version right)

    GLIB_CONFIG=$with_glib/glib-config
    if test -x $GLIB_CONFIG ; then 
      :
    else
      AC_MSG_ERROR([GLIB directory ($with_glib) not present or not configured])
    fi

    # For use in wmpi-config
    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`

    glib_release=`$GLIB_CONFIG --version | sed 's%\\.[[0-9]]*$%%'`

    # canonicalize relative paths
    case $with_glib in 
      /*)
        glib_dir=$with_glib
        ;;
      *)
        glib_dir="\$(top_builddir)/$with_glib"
        ;;
    esac

    GLIB_CFLAGS="-I$glib_dir -I$glib_dir/gmodule"
    GLIB_LIBS="$glib_dir/libglib.la $glib_dir/gmodule.la"
    GLIB_DEPLIBS=

    AC_SUBST(GLIB_CFLAGS)
    AC_SUBST(GLIB_LIBS)
  fi
fi

AC_CHECK_PROG(GTV, gtv-config, yes, no)
if test x$GTV = xyes ; then
   AC_DEFINE(HAVE_GTV, 1, [singular quadrature library available])
   gtv_cflags=`gtv-config --cflags`
   gtv_libs=`gtv-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# General Tetrahedralized Volumes library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/gtv
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(SQT, sqt-config, yes, no)
if test x$SQT = xyes ; then
   AC_DEFINE(HAVE_SQT, 1, [singular quadrature library available])
   sqt_cflags=`sqt-config --cflags`
   sqt_libs=`sqt-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Singular Quadrature on Triangles library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/sqt
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(GQR, gqr-config, yes, no)
if test x$GQR = xyes ; then
   AC_DEFINE(HAVE_GQR, 1, [singular quadrature library available])
   gqr_cflags=`gqr-config --cflags`
   gqr_libs=`gqr-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Gaussian Quadrature Rules library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/gqr
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(MOP, mop-config, yes, no)
if test x$MOP = xyes ; then
   AC_DEFINE(HAVE_MOP, 1,
   		       [multi-variable orthogonal polynomial library available])
   mop_cflags=`mop-config --cflags`
   mop_libs=`mop-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Multi-variable Orthogonal Polynomial library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/mop
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(WBFMM, wbfmm-config, yes, no)
if test x$WBFMM = xyes ; then
   AC_DEFINE(HAVE_WBFMM, 1, [wide band FMM library available])
   wbfmm_cflags=`wbfmm-config --cflags`
   wbfmm_libs=`wbfmm-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Wide Band Fast Multipole Method library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/wbfmm
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(BLASWRAP, blaswrap-config, yes, no)
if test x$BLASWRAP = xyes ; then
   AC_DEFINE(HAVE_BLASWRAP, 1, [BLAS wrapper header available])
   blaswrap_cflags=`blaswrap-config --cflags`
   blaswrap_libs=`blaswrap-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# BLAS wrapper library not installed or not properly installed.
  echo \# You can get it from
  echo \#
  echo \# https://github.com/mjcarley/blaswrap
  echo \#
  echo \#################################################################
  exit
fi

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx -DWBFMM_USE_AVX"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma -DWBFMM_USE_AVX"
fi

AC_CHECK_LIB(m, cos)
AC_CONFIG_MACRO_DIR([m4])

# Checks for library functions.
AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
ACX_BLAS
ACX_LAPACK

blas_libs="$BLAS_LIBS"
lapack_libs="$LAPACK_LIBS"

CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $glib_cflags $wbfmm_cflags $sqt_cflags $gqr_cflags $mop_cflags $gtv_cflags $blaswrap_cflags"
LDFLAGS="$LDFLAGS $glib_libs $wbfmm_libs $sqt_libs $gqr_libs $mop_libs $gtv_libs"

AC_SUBST(blas_libs)
AC_SUBST(lapack_libs)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
#AC_SUBST(DX_RULES)

LIBS="$glib_libs $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 tools/Makefile])
AC_OUTPUT
