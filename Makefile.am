#
# Unlike chroma, using a single top-level Makefile
#

#
# All the include files (for distribution)
#
include_HEADERS = \
     lib/qdp2mma.hpp lib/qdp_utilities.hpp 

#
# Include path and any other compiler flags 
#
AM_CXXFLAGS  = -I$(srcdir)/lib -I$(top_builddir)/lib
AM_CXXFLAGS += @CXXFLAGS@ @QDPXX_CXXFLAGS@

AM_CFLAGS = -I$(srcdir)/lib -I$(top_builddir)/lib 
AM_CFLAGS += @CFLAGS@

lib_LIBRARIES = libganalyze.a

libganalyze_a_SOURCES = \
    lib/qdp2mma.cpp lib/qdp_utilities.cpp 

#
# Include Path
#
INCFLAGS=-I$(srcdir)/programs -I$(srcdir)/lib -I$(top_builddir)/lib
#
AM_LDFLAGS = -L$(BUILDDIR)/lib  @QDPXX_LDFLAGS@ @LDFLAGS@
LDADD = -lganalyze @QDPXX_LIBS@ @LIBS@ 

bin_PROGRAMS = programs/ganalyze
programs_ganalyze_SOURCES = programs/gauge-analyze.cpp
