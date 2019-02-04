#Store macros locally
CLOCAL_AMFLAGS = -I build-aux/m4

#===========================================================================#
EXTRA_DIST=
#Values are filled by Makemodue.am files
# lib_LTLIBRARIES=
# pyexec_LTLIBRARIES=
noinst_LIBRARIES=
noinst_LTLIBRARIES=
bin_PROGRAMS=
check_PROGRAMS=
man1_MANS=
dist_bin_SCRIPTS=
TESTS=

AM_CXXFLAGS=

AM_CPPFLAGS= -I$(srcdir)/lib\
			 -I$(srcdir)/cxxopts/include\
			 -I$(srcdir)/pybind11/include\
			 -I$(srcdir)/casm-utilities/include\
			 -I$(srcdir)/casm-utilities/CASMcode/include\
			 -I$(srcdir)/casm-utilities/CASMcode/include/casm/external\
			 -I$(srcdir)/casm-utilities/CASMcode/include/casm/external/qhull/libqhullcpp\
			 -I$(srcdir)/casm-utilities/CASMcode/include/casm/external/qhull/libqhull_r

AM_LDFLAGS=$(BOOST_LDFLAGS)	#I hate you

#===========================================================================#
#Now do the programs
#You may comment out any of these lines to suppress parts of the installation process
include $(srcdir)/lib/Makemodule.am    #this is where the casm, casm-utilities, and multishift libraries are (not installed, just convenience)
include $(srcdir)/src/Makemodule.am    #this is where the actual executables are
include $(srcdir)/tests/Makemodule.am

#===========================================================================#

var-check:
	echo $(LIBS)
	echo $(PYTHON_LIBS)
	echo $(BOOST_SYSTEM_LIB)
	echo $(BOOST_PYTHON_LIB)
	echo $(BOOST_PYTHON_LIBS)
	echo $(BOOST_LDFLAGS)
	echo $(BOOST_CPPFLAGS)
	echo $(PYTHON)
	echo $(PYTHON_CPPFLAGS)
	echo $(AM_V_CXX)
	echo $(CXXCOMPILE)
	echo $(DEPDIR)
	echo "lol"
