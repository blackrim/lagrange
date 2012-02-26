AC_PREREQ([2.68])
AC_INIT([lagrange], [3.0], [eebsmith@umich.edu])
AC_CONFIG_SRCDIR([BioGeoTree.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC


AC_CHECK_LIB([m], [fabs])

# Checks for libraries.
have_nlopt=no
AC_SEARCH_LIBS([nlopt_create],[nlopt], [have_nlopt=yes])

if test "x${have_nlopt}" = xno; then
   AC_MSG_ERROR([
--------------------------------------------------------
The nlopt library is required for optimization.
--------------------------------------------------------])
fi

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
