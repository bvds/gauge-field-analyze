#
# Unlike chroma, using a single top-level Makefile
#

#
# All the include files (for distribution)
#
include_HEADERS = \
     lib/qdp2mma.h lib/qdp_utilities.h lib/gauge_io.h

#
# Include path and any other compiler flags 
#
AM_CXXFLAGS  = -I$(srcdir)/lib -I$(top_builddir)/lib
AM_CXXFLAGS += @CXXFLAGS@ @QDPXX_CXXFLAGS@

AM_CFLAGS = -I$(srcdir)/lib -I$(top_builddir)/lib 
AM_CFLAGS += @CFLAGS@

lib_LIBRARIES = lib/libganalyze.a

lib_libganalyze_a_SOURCES = \
    lib/qdp2mma.cc lib/qdp_utilities.cc lib/gauge_io.cc

#
# Include Path
#
INCFLAGS=-I$(srcdir)/programs -I$(srcdir)/lib -I$(top_builddir)/lib
#
AM_LDFLAGS = -L$(top_builddir)/lib  @QDPXX_LDFLAGS@ @LDFLAGS@
LDADD = lib/libganalyze.a @QDPXX_LIBS@ @LIBS@ 


bin_PROGRAMS = programs/ganalyze
programs_ganalyze_SOURCES = programs/gauge-analyze.cc
