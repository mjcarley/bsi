/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define if you have an ATLAS BLAS library. */
/* #undef HAVE_ATLAS_BLAS */

/* Define if avx2 instructions are supported */
/* #undef HAVE_AVX2_INSTRUCTIONS */

/* Define if avx instructions are supported */
#define HAVE_AVX_INSTRUCTIONS 1

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* BLAS wrapper header available */
#define HAVE_BLASWRAP 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define if fma instructions are supported */
/* #undef HAVE_FMA_INSTRUCTIONS */

/* singular quadrature library available */
#define HAVE_GQR 1

/* singular quadrature library available */
#define HAVE_GTV 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* multi-variable orthogonal polynomial library available */
#define HAVE_MOP 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* wide band FMM library available */
#define HAVE_WBFMM 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "bsi"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "codes@paraffinalia.co.uk"

/* Define to the full name of this package. */
#define PACKAGE_NAME "bsi"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "bsi 0.1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "bsi"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.1.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "0.1.0"
