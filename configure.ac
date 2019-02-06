#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([multishifter], [0.1], [jgg@ucsb.edu])
AC_CONFIG_SRCDIR([src/base/base.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
#These nested directories are a PITA
AM_INIT_AUTOMAKE([subdir-objects 1.9 tar-ustar])

#This is for dependent libraries that get compiled. You might be able to get rid of this though.
LT_INIT

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11

# Checks for libraries. There is no escaping boost. Boost was a mistake.
BOOSTV=1.52
AX_BOOST_BASE([$BOOSTV],[],[AC_MSG_ERROR(You need Boost $BOOSTV libraries or higher.)])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_CHRONO
AX_BOOST_REGEX
# AX_BOOST_PYTHON
AC_CHECK_LIB([dl], [dlopen], [], AC_MSG_ERROR(dlopen from dl library not found!))

# We won't depend on a library. Instead we'll stuff all the dependencies into the executable itself
# AC_CHECK_LIB([casm],[main], [], [AC_MSG_ERROR(The CASM libraries are required!)])

#You will need pybind11 for this operation
#Why won't it find the libraries???
# AC_CHECK_HEADER([pybind11/pybind11.h],[],[AC_MSG_ERROR(Pybind11 is required for this operation)])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
