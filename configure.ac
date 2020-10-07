#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([multishifter], [0.2.1], [jg.goiri@gmail.com])
AC_CONFIG_SRCDIR([src/multishifter.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
#These nested directories are a PITA
AM_INIT_AUTOMAKE([subdir-objects 1.9 tar-ustar])

LT_INIT

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17

# Checks for libraries.
AX_CHECK_ZLIB([],AC_MSG_ERROR(zlib library not found!))
AX_PTHREAD([],AC_MSG_WARN(pthread is required to run tests!))

AC_SEARCH_LIBS(dlopen, dl, [], AC_MSG_ERROR(dl library not found!))

# Why is this such a pain in the butt? I just want to -lstdc++fs, I could be doing this in the
# Makefile, but instead I'm following this example:
# https://github.com/LibreOffice/online/blob/master/configure.ac
AC_MSG_CHECKING(for std::filesystem features in -lstdc++fs)
STD_FILESYSTEM=
    save_CXXFLAGS=$CXXFLAGS
    CXXFLAGS="$CXXFLAGS -Werror"
    save_LIBS=$LIBS
    LIBS="$save_LIBS -lstdc++fs"
    AC_LANG_PUSH([C++])
    AC_LINK_IFELSE([AC_LANG_SOURCE([[
            #include <filesystem>
            int main()
            {
                if (!std::filesystem::temp_directory_path().empty())
                    return 0;
                return 1;
            }
            ]])],[STD_FILESYSTEM=TRUE], AC_MSG_ERROR(stdc++fs not found!))
    AC_LANG_POP([C++])
    CXXFLAGS=$save_CXXFLAGS
    LIBS=$save_LIBS

if test -n "$STD_FILESYSTEM" ; then
    LIBS="$LIBS -lstdc++fs"
    AC_MSG_RESULT(yes)
else
    AC_MSG_ERROR("libstdc++fs library not found!)
fi

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([tests/regress/common.rc])
AC_CONFIG_FILES([tests/regress/slice/run.sh],[chmod +x tests/regress/slice/run.sh])
AC_CONFIG_FILES([tests/regress/stack/run.sh],[chmod +x tests/regress/stack/run.sh])
AC_OUTPUT
