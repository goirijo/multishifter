#Store macros locally
CLOCAL_AMFLAGS = -I build-aux/m4

#===========================================================================#
#Values are filled by Makemodue.am files
# lib_LTLIBRARIES=
# pyexec_LTLIBRARIES=
noinst_LIBRARIES=
bin_PROGRAMS=
check_PROGRAMS=
man1_MANS=
dist_bin_SCRIPTS=
TESTS=

AM_CXXFLAGS=

AM_CPPFLAGS= -I$(srcdir)/include\
             -I$(srcdir)/CASMcode/include\
             -I$(srcdir)/CASMcode/include/casm/external/qhull/libqhullcpp\
             -I$(srcdir)/CASMcode/include/casm/external/qhull/libqhull_r

AM_LDFLAGS=$(BOOST_LDFLAGS)

#===========================================================================#
#Now do the programs
#You may comment out any of these lines to suppress parts of the installation process
include $(srcdir)/src/Makemodule.am				#utilities library

#===========================================================================#
#Accept .cxpy as a c++ extension
# $(AM_V_CXX)depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.o$$||'`;\
# $(CXXCOMPILE) -MT $@ -MD -MP -MF -c -o $@ $< 
# SUFFIXES: .cxpy
# .cxpy.o:
# 	$(CXXCOMPILE) -c -o $@ $<

# .cxpy.obj:
# 	$(AM_V_CXX)depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.obj$$||'`;\
# 	$(CXXCOMPILE) -MT $@ -MD -MP -MF -c -o $@ `$(CYGPATH_W) '$<'` 

# .cxpy.lo:
# 	$(AM_V_CXX)depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.lo$$||'`;\
# 	$(LTCXXCOMPILE) -MT $@ -MD -MP -MF -c -o $@ $< 


#===========================================================================#
EXTRA_DIST=

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
