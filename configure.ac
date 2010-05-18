dnl Process this file with autoconf to produce a configure script.
AC_INIT(dafty, 0, roundyd@physics.oregonstate.edu, dafty)
AC_CONFIG_SRCDIR(src/dafty.cpp)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(src/config.h)
# use c++ as a compiler in the tests
AC_LANG(C++)

dnl Set default CXX flags
if test "x$CXXFLAGS" = "x"; then
  CXXFLAGS='-ansi -g -pipe -W -Wall -O2'
fi

# check whether to enable profiling or coverage testing.
AC_ARG_ENABLE(profiling, [AC_HELP_STRING([--enable-profiling],
                         [use profiling])],
                  use_profiling=$enableval, use_profiling=no)
if test "x$use_profiling" = "xyes"; then
   CXXFLAGS="$CXXFLAGS -pg"
fi
AC_ARG_ENABLE(coverage, [AC_HELP_STRING([--enable-coverage],
                         [use coverage testing])],
                  use_coverage=$enableval, use_coverage=no)
if test "x$use_coverage" = "xyes"; then
   CXXFLAGS="$CXXFLAGS -fprofile-arcs -ftest-coverage"
fi

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# and support large files...
AC_SYS_LARGEFILE

dnl Checks for libraries.

AC_CHECK_LIB(fftw3, fftw_plan_r2r_1d, [I_HAVE_FFTW3=1],
             [AC_MSG_WARN([could not find fftw3 library!])
              I_HAVE_FFTW3=0])
AC_CHECK_LIB(fftw, fftw_create_plan, [I_HAVE_FFTW2=1],
             [AC_CHECK_LIB(fftw -lm, fftw_create_plan, [I_HAVE_FFTW2=1],
                           [AC_MSG_WARN([could not find fftw library!])
                            I_HAVE_FFTW2=0])])

AC_ARG_WITH(fftw2, [AC_HELP_STRING([--without-fftw2],[do not use fftw 2])],
                    with_fftw2=$withval, with_fftw2=yes)
AC_ARG_WITH(fftw3, [AC_HELP_STRING([--with-fftw3],[use fftw 3])],
                    with_fftw3=$withval, with_fftw3=no)
if test "x$with_fftw2" = "xno" -a $I_HAVE_FFTW3 == 1; then
   I_HAVE_FFTW2=0
fi
if test "x$with_fftw3" = "xno" -a $I_HAVE_FFTW2 == 1; then
   I_HAVE_FFTW3=0
fi

AC_MSG_CHECKING([for fft library])
if test $I_HAVE_FFTW2 == 1 -a $I_HAVE_FFTW3 == 1; then
  AC_MSG_RESULT([both fftw and fftw3... using fftw3])
  AC_DEFINE(DFT_USE_FFTW3, 1, [Define to 1 if you want to use fftw3])
  LIBS="$LIBS -lfftw3"
elif test $I_HAVE_FFTW3 == 1; then
  AC_MSG_RESULT([fftw3])
  AC_DEFINE(DFT_USE_FFTW3, 1, [Define to 1 if you want to use fftw3])
  LIBS="$LIBS -lfftw3"
elif test $I_HAVE_FFTW2 == 1; then
  AC_MSG_RESULT([fftw])
  AC_DEFINE(DFT_USE_FFTW2, 1, [Define to 1 if you want to use fftw2])
  LIBS="$LIBS -lfftw"
else
  AC_MSG_RESULT([none!])
  AC_MSG_ERROR([You need some sort of an fft library!])
fi

# Do we want threads?

AC_ARG_WITH(threads,
            [AC_HELP_STRING([--without-threads],
                            [disable pthreads parallelization])],
            with_threads=$withval, with_threads=maybe)


if test "x$with_threads" = "xno"; then
  AC_MSG_WARN([Not using pthreads...])
else
  AC_CHECK_LIB(pthread, pthread_join,
    [with_threads=yes
     LIBS="$LIBS -lpthread"
     AC_DEFINE(DFT_THREAD, 1,
               [Define to 1 if you want to use threads])],
    [if test "x$with_threads" = "xyes"; then
       AC_MSG_ERROR([Cannot find pthreads library!])
     fi]
    )
fi

dnl Checks for header files.
AC_STDC_HEADERS
AC_HAVE_HEADERS(fcntl.h sys/time.h)

AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
