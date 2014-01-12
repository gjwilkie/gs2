#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([simpledataio], [1.0.0], [edmundhighcock@users.sourceforge.net])
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


# Check whether to build it using mpi compilers (must be used for enable_parallel)
if test x"$enable_parallel" = xyes; then
enable_mpi=yes
else
AC_MSG_CHECKING([whether to build using mpi (automatically enabled by --enable-parallel)])
AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi],
                              [build with mpi version of compilers])],
							[enable_mpi=yes], [enable_mpi=no])
AC_MSG_RESULT([$enable_mpi])
fi
AM_CONDITIONAL(WITHMPI, [test "x$enable_mpi" = xyes])

#AC_ARG_VAR([CC], [C Compiler])
# Checks for programs. If enable parallel, find the appropriate mpi programs, libraries etc.
AX_PROG_CC_MPI([test "x$enable_mpi" = xyes], [], [
    if test x"$enable_mpi" = xyes; then
      AC_MSG_FAILURE([MPI C compiler requested, but could not use MPI.])
    else
      AC_MSG_WARN([No MPI C compiler found, not using MPI.])
    fi
  ])
AX_PROG_FC_MPI([test "x$enable_mpi" = xyes], [], [
    if test x"$enable_mpi" = xyes; then
      AC_MSG_FAILURE([MPI Fortran compiler requested, but could not use MPI.])
    else
      AC_MSG_WARN([No MPI Fortran compiler found, not using MPI.])
    fi
  ])

# If we are not using the GNU compiler, we are probably on an HPC system
# in which case, we force static compilation of the test programs
AX_COMPILER_VENDOR
AM_CONDITIONAL(FORCE_STATIC, [test "x$ax_cv_c_compiler_vendor" != xgnu])
#AC_MSG_WARN([compiler vendor is $ax_cv_c_compiler_vendor])

# Check whether the Fortran compiler produces seperate .mod files
AC_FC_MODULE_EXTENSION
AM_CONDITIONAL(HAS_F90_MODULES, [test "x$FC_MODEXT" != x])

AL_CHECK_ISO_C_BINDING
if test "x$nc_iso_c_binding" = xno; then 
	AC_MSG_WARN("The Fortran compiler does not support the intrinsic iso_c_binding module: the fortran netCDF interface will not work")
fi
AM_CONDITIONAL(HAS_ISO_C_BINDING, test "x$nc_iso_c_binding" = xyes)


LT_INIT
#AC_PROG_FC
#AC_PROG_FPP
#AC_PROG_RANLIB
#AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([hdf5], [H5open])

#AC_CHECK_HEADERS([netcdf.mod])
#AC_CHECK_HEADERS([simpledataio.h])

AC_CHECK_HEADERS([stdlib.h])

# Check whether NetCDF builds OK
AL_CHECK_NETCDF
if test "x$al_netcdf_ok" = xno; then
	# If not,
	# try to automatically find the netcdf prefix
	AC_MSG_RESULT([could not find netcdf... trying to setup netcdf using nc-config])
	if AC_TRY_EVAL("nc-config --prefix > conftest.tmp "); then
		AC_MSG_CHECKING([netcdf prefix])
		al_netcdf_prefix=`nc-config --prefix`
		al_netcdf_found_prefix=yes
		AC_MSG_RESULT([$al_netcdf_prefix])
		al_netcdf_libs=`nc-config --libs`
		al_netcdf_inc=`nc-config --includedir`
		FCFLAGS+=" -I$al_netcdf_inc"
		CFLAGS+=" -I$al_netcdf_inc"
		LDFLAGS+=" $al_netcdf_libs"


	else
		AC_MSG_WARN("Couldn't automatically detect netcdf prefix: nc-config wasn't available")
		al_netcdf_found_prefix=no
	fi
	AL_CHECK_NETCDF
	if test "x$al_netcdf_ok" = xno; then
		AC_MSG_FAILURE("Failed to find all netCDF prerequisites")
	fi
fi

if test "x$al_netcdff_present" = xno; then
	AC_MSG_WARN("Cannot find the fortran netcdf library: the fortran netCDF interface will not work")
fi
AM_CONDITIONAL(HAS_FORTRAN_NETCDF, test "x$al_netcdff_present" = xyes)


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
