#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([pialign], [0.1.0], [neubig@gmail.com])
AC_CONFIG_SRCDIR(src/lib/pialign.cc)
AC_CONFIG_HEADERS([src/include/pialign/config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    src/include/Makefile
    src/lib/Makefile
    src/bin/Makefile
    src/test/Makefile
])

AM_INIT_AUTOMAKE([-Wall])

# Uncomment to disable shared libraries (for example, for debugging)
AM_DISABLE_SHARED

# Checks for programs.
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_PROG_CC

AC_LANG_CPLUSPLUS
LDFLAGS="$LDFLAGS -lpthread"

# Checks for libraries.

AC_SYS_LARGEFILE

# # Boost libraries...
# orig_CPPFLAGS="$CPPFLAGS"
# orig_LDFLAGS="$LDFLAGS"
# CPPFLAGS="$orig_CPPFLAGS $ICU_CPPFLAGS"
# LDFLAGS="$orig_LDFLAGS $ICU_LDFLAGS"
# 
# BOOST_REQUIRE([1.38])
# BOOST_PROGRAM_OPTIONS
# BOOST_FILESYSTEM
# BOOST_IOSTREAMS
# BOOST_SYSTEM
# 
# CPPFLAGS="$orig_CPPFLAGS"
# LDFLAGS="$orig_LDFLAGS"
# 
# ## save pthred-flag, computed by-product of boost-threads
# PTHREAD_FLAGS=$boost_cv_pthread_flag
# AC_SUBST(PTHREAD_FLAGS)

# Checks for header files.
AC_CHECK_HEADERS([sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday pow sqrt])

AC_OUTPUT
