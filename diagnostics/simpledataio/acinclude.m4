

# AC_PROG_FC_MOD
# By Michael R Nolta
# ---------------
AC_DEFUN([AC_PROG_FC_UPPERCASE_MOD],
[
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([if Fortran 90 compiler capitalizes .mod filenames])
		    cat <<EOF >conftest.f90
		      module conftest
		      end module conftest
EOF
ac_try='${FC} ${FCFLAGS} conftest.f90 >&AS_MESSAGE_LOG_FD'
AC_TRY_EVAL(ac_try)
if test -f CONFTEST.mod ; then
   ac_cv_prog_f90_uppercase_mod=yes
   rm -f CONFTEST.mod
else
   ac_cv_prog_f90_uppercase_mod=no
fi
AC_MSG_RESULT($ac_cv_prog_f90_uppercase_mod)
#rm -f conftest*
AC_LANG_POP(Fortran)
])


# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_prog_cc_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_CC_MPI([MPI-WANTED-TEST[, ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile C programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/). The macro has to
#   be used instead of the standard macro AC_PROG_CC and will replace the
#   standard variable CC with the found compiler.
#
#   MPI-WANTED-TEST is used to test whether MPI is actually wanted by the
#   user. If MPI-WANTED_TEST is omitted or if it succeeds, the macro will
#   try to find out how to use MPI, if it fails, the macro will call
#   AC_PROG_CC to find a standard C compiler instead.
#
#   When MPI is found, ACTION-IF-FOUND will be executed, if MPI is not found
#   (or MPI-WANTED-TEST fails) ACTION-IF-NOT-FOUND is executed. If
#   ACTION-IF-FOUND is not set, the macro will define HAVE_MPI.
#
#   The following example demonstrates usage of the macro:
#
#     # If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
#     # If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
#     # If --with-mpi=no is used, use a standard C compiler instead.
#     AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
#         [compile with MPI (parallelization) support. If none is found,
#         MPI is not used. Default: auto])
#     ],,[with_mpi=auto])
#     #
#     AX_PROG_CC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
#       use_mpi=no
#       if test x"$with_mpi" = xyes; then
#         AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
#       else
#         AC_MSG_WARN([No MPI compiler found, won't use MPI.])
#       fi
#     ])
#
# LICENSE
#
#   Copyright (c) 2010,2011 Olaf Lenz <olenz@icp.uni-stuttgart.de>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.


AC_DEFUN([AX_PROG_CC_MPI], [
AC_PREREQ(2.50)

# Check for compiler
# Needs to be split off into an extra macro to ensure right expansion
# order.
AC_REQUIRE([_AX_PROG_CC_MPI],[_AX_PROG_CC_MPI([$1])])

AS_IF([test x"$_ax_prog_cc_mpi_mpi_wanted" = xno],
  [ _ax_prog_cc_mpi_mpi_found=no ],
  [
    AC_LANG_PUSH([C])
    # test whether MPI_Init is available
    # We do not use AC_SEARCH_LIBS here, as it caches its outcome and
    # thus disallows corresponding calls in the other AX_PROG_*_MPI
    # macros.
    for lib in NONE mpi mpich; do
      save_LIBS=$LIBS
      if test x"$lib" = xNONE; then
        AC_MSG_CHECKING([for function MPI_Init])
      else
        AC_MSG_CHECKING([for function MPI_Init in -l$lib])
        LIBS="-l$lib $LIBS"
      fi
      AC_LINK_IFELSE([AC_LANG_CALL([],[MPI_Init])],
        [ _ax_prog_cc_mpi_mpi_found=yes ],
        [ _ax_prog_cc_mpi_mpi_found=no ])
      AC_MSG_RESULT($_ax_prog_cc_mpi_mpi_found)
      if test "x$_ax_prog_cc_mpi_mpi_found" = "xyes"; then
        break;
      fi
      LIBS=$save_LIBS
    done

    # Check for header
    AS_IF([test x"$_ax_prog_cc_mpi_mpi_found" = xyes], [
      AC_MSG_CHECKING([for mpi.h])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <mpi.h>])],
        [ AC_MSG_RESULT(yes)],
        [ AC_MSG_RESULT(no)
         _ax_prog_cc_mpi_mpi_found=no
      ])
    ])
    AC_LANG_POP([C])
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test x"$_ax_prog_cc_mpi_mpi_found" = xyes], [
        ifelse([$2],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$2])
        :
],[
        $3
        :
])

])dnl AX_PROG_CC_MPI

dnl _AX_PROG_CC_MPI is an internal macro required by AX_PROG_CC_MPI.
dnl To ensure the right expansion order, the main function AX_PROG_CC_MPI
dnl has to be split into two parts.
dnl
dnl Known MPI C compilers:
dnl  mpicc
dnl  mpixlc_r
dnl  mpixlc
dnl  hcc
dnl  mpxlc_r
dnl  mpxlc
dnl  sxmpicc  NEC SX
dnl  mpifcc   Fujitsu
dnl  mpgcc
dnl  mpcc
dnl  cmpicc
dnl  cc
dnl
AC_DEFUN([_AX_PROG_CC_MPI], [
  AC_ARG_VAR(MPICC,[MPI C compiler command])
  ifelse([$1],,[_ax_prog_cc_mpi_mpi_wanted=yes],[
    AC_MSG_CHECKING([whether to compile using MPI])
    if $1; then
      _ax_prog_cc_mpi_mpi_wanted=yes
    else
      _ax_prog_cc_mpi_mpi_wanted=no
    fi
    AC_MSG_RESULT($_ax_prog_cc_mpi_mpi_wanted)
  ])
  if test x"$_ax_prog_cc_mpi_mpi_wanted" = xyes; then
    if test -z "$CC" && test -n "$MPICC"; then
      CC="$MPICC"
    else
      AC_CHECK_TOOLS([CC], [mpicc mpixlc_r mpixlc hcc mpxlc_r mpxlc sxmpicc mpifcc mpgcc mpcc cmpicc cc gcc])
    fi
  fi
  AC_PROG_CC
])dnl _AX_PROG_CC_MPI


# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_prog_fc_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_FC_MPI([MPI-WANTED-TEST[, ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile Fortran77 programs that use
#   MPI (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/).  The macro has to
#   be used instead of the standard macro AC_PROG_FC and will replace the
#   standard variable FC with the found compiler.
#
#   MPI-WANTED-TEST is used to test whether MPI is actually wanted by the
#   user. If MPI-WANTED_TEST is omitted or if it succeeds, the macro will
#   try to find out how to use MPI, if it fails, the macro will call
#   AC_PROG_CC to find a standard C compiler instead.
#
#   When MPI is found, ACTION-IF-FOUND will be executed, if MPI is not found
#   (or MPI-WANTED-TEST fails) ACTION-IF-NOT-FOUND is executed. If
#   ACTION-IF-FOUND is not set, the macro will define HAVE_MPI.
#
#   The following example demonstrates usage of the macro:
#
#     # If --with-mpi=auto is used, try to find MPI, but use standard FC compiler if it is not found.
#     # If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
#     # If --with-mpi=no is used, use a standard FC compiler instead.
#     AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
#         [compile with MPI (parallelization) support. If none is found,
#         MPI is not used. Default: auto])
#     ],,[with_mpi=auto])
#
#     AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
#       use_mpi=no
#       if test x"$with_mpi" = xyes; then
#         AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
#       else
#         AC_MSG_WARN([No MPI compiler found, won't use MPI.])
#       fi
#     ])
#
# LICENSE
#
#   Copyright (c) 2010,2011 Olaf Lenz <olenz@icp.uni-stuttgart.de>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.


AC_DEFUN([AX_PROG_FC_MPI], [
AC_PREREQ(2.50)

# Check for compiler
# Needs to be split off into an extra macro to ensure right expansion
# order.
AC_REQUIRE([_AX_PROG_FC_MPI],[_AX_PROG_FC_MPI([$1])])

AS_IF([test x"$_ax_prog_fc_mpi_mpi_wanted" = xno],
  [ _ax_prog_fc_mpi_mpi_found=no ],
  [
    AC_LANG_PUSH([Fortran])

    # test whether MPI_INIT is available
    # We do not use AC_SEARCH_LIBS here, as it caches its outcome and
    # thus disallows corresponding calls in the other AX_PROG_*_MPI
    # macros.
    for lib in NONE mpichf90 fmpi fmpich; do
      save_LIBS=$LIBS
      if test x"$lib" = xNONE; then
        AC_MSG_CHECKING([for function MPI_INIT])
      else
        AC_MSG_CHECKING([for function MPI_INIT in -l$lib])
        LIBS="-l$lib $LIBS"
      fi
      AC_LINK_IFELSE([AC_LANG_CALL([],[MPI_INIT])],
        [ _ax_prog_fc_mpi_mpi_found=yes ],
        [ _ax_prog_fc_mpi_mpi_found=no ])
      AC_MSG_RESULT($_ax_prog_fc_mpi_mpi_found)
      if test "x$_ax_prog_fc_mpi_mpi_found" = "xyes"; then
        break;
      fi
      LIBS=$save_LIBS
    done

    # Check for header
    AS_IF([test x"$_ax_prog_fc_mpi_mpi_found" = xyes], [
      AC_MSG_CHECKING([for mpif.h])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[
      include 'mpif.h'
]])],
        [ AC_MSG_RESULT(yes)],
        [ AC_MSG_RESULT(no)
	  _ax_prog_fc_mpi_mpi_found=no
      ])
    ])
    AC_LANG_POP([Fortran])
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test x"$_ax_prog_fc_mpi_mpi_found" = xyes], [
        ifelse([$2],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$2])
        :
],[
        $3
        :
])

])dnl AX_PROG_FC_MPI

dnl _AX_PROG_FC_MPI is an internal macro required by AX_PROG_FC_MPI.
dnl To ensure the right expansion order, the main function AX_PROG_FC_MPI
dnl has to be split into two parts. This part looks for the MPI
dnl compiler, while the other one tests whether an MPI program can be
dnl compiled.
dnl
AC_DEFUN([_AX_PROG_FC_MPI], [
  AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
  ifelse([$1],,[_ax_prog_fc_mpi_mpi_wanted=yes],[
    AC_MSG_CHECKING([whether to compile using MPI])
    if $1; then
      _ax_prog_fc_mpi_mpi_wanted=yes
    else
      _ax_prog_fc_mpi_mpi_wanted=no
    fi
    AC_MSG_RESULT($_ax_prog_fc_mpi_mpi_wanted)
  ])
  if test x"$_ax_prog_fc_mpi_mpi_wanted" = xyes; then
    if test -z "$FC" && test -n "$MPIFC"; then
      FC="$MPIFC"
    else
      AC_CHECK_TOOLS([FC], [mpif95 mpxlf95_r mpxlf95 ftn mpif90 mpxlf90_r mpxlf90 mpf90 cmpif90c sxmpif90 mpif77 hf77 mpxlf_r mpxlf mpifrt mpf77 cmpifc xlf95 pgf95 pathf95 ifort g95 f95 fort ifc efc openf95 sunf95 crayftn gfortran lf95 ftn xlf90 f90 pgf90 pghpf pathf90 epcf90 sxf90 openf90 sunf90 xlf f77 frt pgf77 pathf77 g77 cf77 fort77 fl32 af77])
    fi
  fi
  AC_PROG_FC
])dnl _AX_PROG_FC_MPI
	
# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_compiler_vendor.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_COMPILER_VENDOR
#
# DESCRIPTION
#
#   Determine the vendor of the C/C++ compiler, e.g., gnu, intel, ibm, sun,
#   hp, borland, comeau, dec, cray, kai, lcc, metrowerks, sgi, microsoft,
#   watcom, etc. The vendor is returned in the cache variable
#   $ax_cv_c_compiler_vendor for C and $ax_cv_cxx_compiler_vendor for C++.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.


AC_DEFUN([AX_COMPILER_VENDOR],
[AC_CACHE_CHECK([for _AC_LANG compiler vendor], ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor,
  [# note: don't check for gcc first since some other compilers define __GNUC__
  vendors="intel:     __ICC,__ECC,__INTEL_COMPILER
           ibm:       __xlc__,__xlC__,__IBMC__,__IBMCPP__
           pathscale: __PATHCC__,__PATHSCALE__
           clang:     __clang__
           fujitsu:   __FUJITSU
           gnu:       __GNUC__
           sun:       __SUNPRO_C,__SUNPRO_CC
           hp:        __HP_cc,__HP_aCC
           dec:       __DECC,__DECCXX,__DECC_VER,__DECCXX_VER
           borland:   __BORLANDC__,__TURBOC__
           comeau:    __COMO__
           cray:      _CRAYC
           kai:       __KCC
           lcc:       __LCC__
           sgi:       __sgi,sgi
           microsoft: _MSC_VER
           metrowerks: __MWERKS__
           watcom:    __WATCOMC__
           portland:  __PGI
           unknown:   UNKNOWN"
  for ventest in $vendors; do
    case $ventest in
      *:) vendor=$ventest; continue ;;
      *)  vencpp="defined("`echo $ventest | sed 's/,/) || defined(/g'`")" ;;
    esac
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
      #if !($vencpp)
        thisisanerror;
      #endif
    ])], [break])
  done
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor=`echo $vendor | cut -d: -f1`
 ])
])

