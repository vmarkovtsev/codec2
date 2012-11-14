#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT(codec2, 0.2, david@rowetel.com)
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_BUILD
AC_CANONICAL_SYSTEM
AC_LANG([C])
# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE(1.11 foreign)
AM_SILENT_RULES([yes])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h malloc.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_ARG_ENABLE([libavcodec-fft],
              AS_HELP_STRING([--enable-libavcodec-fft], [Use libavcodec FFT implementation instead of the built-in kiss_fft]),
              [enable_libavcodec_fft=yes CFLAGS="$CFLAGS -DLIBAVCODEC_FFT"],
              [enable_libavcodec_fft=no CFLAGS="$CFLAGS -DKISS_FFT"])
AM_CONDITIONAL(LIBAVCODEC_FFT, test x$enable_libavcodec_fft = xyes) 

AS_IF([test "x$enable_libavcodec_fft" == "xyes"], [
# livavcodec FFT check
PKG_PROG_PKG_CONFIG 
PKG_CHECK_MODULES([LIBAVCODEC], [libavcodec]) 
])  

AC_ARG_ENABLE([neon],
              AS_HELP_STRING([--enable-neon], [Enable ARM NEON SIMD optimizations]),
              [enable_neon=yes CFLAGS="$CFLAGS -mfpu=neon -DNEON"],
              [enable_neon=no])
AS_IF([test "x$enable_neon" == "xyes"], [
	AC_CHECK_HEADERS([arm_neon.h])
])     

AC_CONFIG_FILES([Makefile src/Makefile unittest/Makefile])
AC_OUTPUT
