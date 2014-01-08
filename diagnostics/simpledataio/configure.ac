#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([simpledataio], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([])
AC_CONFIG_SRCDIR([src/simpledataio.c])

# Check whether to build it with support for parallel netcdf IO
AC_MSG_CHECKING([whether to build using parallel netcdf])
AC_ARG_ENABLE([parallel],
              [AS_HELP_STRING([--enable-parallel],
                              [build with parallel netcdf (requires netcdf4, hdf5 and mpi)])],
							[enable_parallel=yes], [enable_parallel=no])
AC_MSG_RESULT([$enable_parallel])
AM_CONDITIONAL(PARALLEL, [test "x$enable_parallel" = xyes])

#AC_ARG_VAR([CC], [C Compiler])
# Checks for programs. If enable parallel, find the appropriate mpi programs, libraries etc.
AX_PROG_CC_MPI([test "x$enable_parallel" = xyes], [], [
    if test x"$enable_parallel" = xyes; then
      AC_MSG_FAILURE([MPI C compiler requested, but could not use MPI.])
    else
      AC_MSG_WARN([No MPI C compiler found, not using MPI.])
    fi
  ])
AX_PROG_FC_MPI([test "x$enable_parallel" = xyes], [], [
    if test x"$enable_parallel" = xyes; then
      AC_MSG_FAILURE([MPI Fortran compiler requested, but could not use MPI.])
    else
      AC_MSG_WARN([No MPI Fortran compiler found, not using MPI.])
    fi
  ])

AC_FC_MODULE_EXTENSION
AM_CONDITIONAL(HAS_F90_MODULES, [test "x$FC_MODEXT" != x])


LT_INIT
#AC_PROG_FC
#AC_PROG_FPP
#AC_PROG_RANLIB
#AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([hdf5], [H5open])
AC_CHECK_LIB([netcdf], [nc_create])
AC_CHECK_LIB([netcdff], [nc_create])
AC_CHECK_DECL([nc_create], [], AC_MSG_FAILURE([Cannot find function nc_create ]), [
							 #include "netcdf.h"
							 ])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([netcdf.h], [], [echo "netcdf.h not found... suggest defining CFLAGS='-I/location/of/netcdf/include'" && exit])
#AC_CHECK_HEADERS([netcdf.mod])
#AC_CHECK_HEADERS([simpledataio.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC


if test "$enable_parallel" = yes; then
  AC_CHECK_HEADERS([mpi.h], [], [echo "mpi.h not found... suggest defining CFLAGS='-I/location/of/mpi/include'" && exit])
	# Check that we have the function nc_create_par available... i.e. that we 
	# have the parallel netcdf libraries
  AC_CHECK_LIB([netcdf], [nc_create_par], [], 
		AC_MSG_FAILURE([Cannot find symbol nc_create_par: parallel support requires netCDF 4 built on top of parallel HDF5]))
		AC_CHECK_DECL([nc_create_par], [], 
								AC_MSG_FAILURE([Cannot find function nc_create_par in  netcdf_par.h ]), [
								 #include "netcdf.h"
								 #include "netcdf_par.h"
								 ])

fi

AC_PROG_FC_UPPERCASE_MOD
AM_CONDITIONAL(UPPER_CASE_MOD, [test "x$ac_cv_prog_f90_uppercase_mod" = xyes])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
