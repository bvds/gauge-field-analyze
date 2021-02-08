#
#
#
AC_INIT([gauge-field-analyze],[1.0],[bvds@asu.edu])

# Store the auxiliary build tools (e.g., install-sh, config.sub, config.guess)
# in this dir
AC_CONFIG_AUX_DIR(config)

AC_PROG_CXX(g++ cxx CC)
AC_PROG_CC(gcc)
AM_PROG_AR
AC_PROG_RANLIB
# BvdS:  do we need this?
AC_CHECK_TOOL(AR, ar, [ar])

# Initialise automake
AM_INIT_AUTOMAKE([-Wall -Wno-error foreign subdir-objects])

# defines...
AC_CONFIG_HEADERS([include/ganalyze_config_internal.h])

#
#
# --with-xxxx and --enable-xxxx switches 
#

# --with-qdp=DIR
AC_ARG_WITH(qdp,
  AC_HELP_STRING(
    [--with-qdp=DIR],
    [Link to QDP++ libraries installed in DIR [/usr/local]]
  ),
  [QDPXX_HOME="$with_qdp"],
  [QDPXX_HOME="/usr/local"]
)

#
#  Test QDP install (from chroma configure.ac).
#

if test "X${QDPXX_HOME}X" = "XX" ; then
  AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
else
  AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
fi

if test "X${QDPXX_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
fi

AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
AC_SUBST(CXX, "`${QDPXX_CONFIG} --cxx`")
AC_MSG_NOTICE([QDP++ compiler: ${CXX}])
AC_SUBST(QDPXX_CXXFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
AC_MSG_NOTICE([QDP++ compile flags: ${QDPXX_CXXFLAGS}])
AC_SUBST(QDPXX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
AC_MSG_NOTICE([QDP++ linking flags: ${QDPXX_LDFLAGS}])
AC_SUBST(QDPXX_LIBS,     "`${QDPXX_CONFIG} --libs`")
AC_MSG_NOTICE([QDP++ libraries flags: ${QDPXX_LIBS}])
AC_SUBST(QDPXX_ND,       "`${QDPXX_CONFIG} --Nd`")
if test "X${QDPXX_ND}X" = "XX" ; then
  AC_MSG_ERROR([Cannot determine QDP++ spacetime dimension])
else
  AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
fi
AC_SUBST(QDPXX_NC,       "`${QDPXX_CONFIG} --Nc`")
AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])
AC_SUBST(QDPXX_NS,       "`${QDPXX_CONFIG} --Ns`")
if test "X${QDPXX_NS}X" = "XX" ; then
  AC_MSG_ERROR([Cannot determine number of spin components in QDP++])
else
  AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])
fi

AC_SUBST(QDPXX_PARALLEL_ARCH, "`${QDPXX_CONFIG} --parallel-arch`")
AC_MSG_NOTICE([QDP++ parallel arch: ${QDPXX_PARALLEL_ARCH}])

AC_SUBST(QDPXX_PRECISION, "`${QDPXX_CONFIG} --precision`")
AC_MSG_NOTICE([QDP++ precision: ${QDPXX_PRECISION}])

# Try to compile a QDP++ program to check the --with options
# This function is defined in acinclude.m4
AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
PAC_QDPXX_LINK_CXX_FUNC(
  ${QDPXX_CXXFLAGS},
  ${QDPXX_LDFLAGS},
  ${QDPXX_LIBS},
  ,
  ,
  [qdpxx_link_ok=yes],
  [qdpxx_link_ok=no]
)
if test "X${qdpxx_link_ok}X" = "XyesX"; 
then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with QDP++. Use --with-qdp++=<dir> to select a working version.])
fi 

#######################
# Produce output
#######################

#
# Generate makefiles
#
AC_CONFIG_FILES([Makefile])
#
# Finish
#
AC_OUTPUT
