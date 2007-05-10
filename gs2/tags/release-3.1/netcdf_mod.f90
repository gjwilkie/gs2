module netcdf_mod
!
! modern Fortran90 module wrapper around netcdf routines
!
! This provides a somewhat simplier interface to the netcdf routines,
! and also automatically generates error messages without the calling
! program having to check the status variable after every call.
!
! Greg Hammett, February 1999
!
! NetCDF documentation: http://www.unidata.ucar.edu/packages/netcdf/guidef/
! NetCDF Fortran interface summary: 
!         http://www.unidata.ucar.edu/packages/netcdf/guidef/guidef-19.html
!
!######################################################################
!######################################################################
!
! For example, this module defines netcdf_put_var( ncid, varid, var),
! which in turn calls nf_put_var_int, nf_put_var_double, or nf_put_var_text,
! depending on the type of the var argument.
!
! All of the routines here have identical arguments as the netcdf routines
! they call (so see the netcdf documentation for definitions), with the
! one exception that in netcdf_put_att the length of a text string does 
! not need to be passed.
!
! Example of usage:
!
! status = netcdf_def_var( ncid, 'density', nf_double, 5, ion_dim, density_id)
!
! Common arguments are demonstrated for netcdf_def_var:
!
! INTEGER FUNCTION NETCDF_DEF_VAR(INTEGER NCID, CHARACTER*(*) NAME,
!                            INTEGER XTYPE, INTEGER NVDIMS,
!                            INTEGER VDIMS(*), INTEGER varid)
! has the same arguments as
!
! INTEGER FUNCTION NF_DEF_VAR(INTEGER NCID, CHARACTER*(*) NAME,
!                            INTEGER XTYPE, INTEGER NVDIMS,
!                            INTEGER VDIMS(*), INTEGER varid)
!
! Where
!
! Integer NCID         ! netCDF ID of an open netCDF dataset
! Character*(*) NAME   ! name of variable in netcdf file
! Integer XTYPE        ! type of netcdf variable: NF_DOUBLE, NF_INT, ...
! Integer NVDIMS       ! of dimensions for this variabe
! Integer VDIMS(NVDIMS) ! size of this variable in each dimension
! Integer VARID        ! ID # to identify this variable in the netcdf file
!
!######################################################################
! 
! Routines for WRITING to a NetCDF file:
!
! Define a dimension (the dimension id # is returned as dimid):
!
! INTEGER FUNCTION netcdf_DEF_DIM (INTEGER NCID, CHARACTER*(*) NAME,
!                          INTEGER LEN, INTEGER dimid)
!
! Define a variable:
!
! INTEGER FUNCTION netcdf_DEF_VAR(INTEGER NCID, CHARACTER*(*) NAME,
!                             INTEGER XTYPE, INTEGER NVDIMS,
!                             INTEGER VDIMS(*), INTEGER varid)
! 
! Put (write to the file) an attribute (usually text):
!
! INTEGER FUNCTION  netcdf_PUT_ATT (INTEGER NCID, INTEGER VARID,
!                                   CHARACTER*(*) NAME, CHARACTER*(*) TEXT)
! 
! Write a scalar variable or a whole array:
!
! INTEGER FUNCTION netcdf_PUT_VAR(INTEGER NCID, INTEGER VARID,
!                                    [double, int, or Text] VALS(*))
! 
! Write a subsection of an array:
!
! INTEGER FUNCTION netcdf_PUT_VARA(INTEGER NCID, INTEGER VARID,
!                                  INTEGER START(*), INTEGER COUNT(*),
!                                  [double or int] VALS(*)) 
!
!######################################################################
! Routines for READING from a NetCDF file:
!
! Given a dimension's name, find out (inquire) the dimension's id #, dimid
!
! INTEGER FUNCTION netcdf_INQ_DIMID (INTEGER NCID, CHARACTER*(*) NAME,
!                               INTEGER dimid)
!
! Given a variables's name, find out (inquire) the variable's id #, varid
! 
! INTEGER FUNCTION netcdf_INQ_VARID(INTEGER NCID, CHARACTER*(*) NAME, 
!                              INTEGER varid)
!
! Read a netcdf variable from a file (scalar or a whole array):
!
! INTEGER FUNCTION netcdf_GET_VAR(INTEGER NCID, INTEGER VARID,
!                                [double or int or tet] vals(*))
!
!######################################################################
!######################################################################

  use mp, only: proc0, iproc, broadcast
  use constants
  implicit none
  private
  include 'netcdf.inc'

!######################################################################
! Routines for WRITING to a NetCDF file:

  public :: netcdf_init     ! initialize, select serial/parallel
  public :: netcdf_def_dim  ! define a NetCDF dimension
  public :: netcdf_def_var  ! define a NetCDF variable
  public :: netcdf_put_att  ! write a NetCDF attribute to file
  public :: netcdf_put_var  ! write a NetCDF variable or array to file
  public :: netcdf_put_var1  ! write one element of a NetCDF array to file
  public :: netcdf_put_vara ! write a section of a NetCDF array to file

!######################################################################
! Routines for READING from a NetCDF file:

  public :: netcdf_inq_dimid ! find a dimensions's id # from its name
  public :: netcdf_inq_varid ! find a variable's id # from its name
  public :: netcdf_get_var   ! read a NetCDF variable from a file

!######################################################################

  logical, parameter :: debug = .false.  ! true to print diagnostic info
  logical, private, save :: serial_io = .false.
  logical, private, save :: proc_write

  interface netcdf_def_var
     module procedure netcdf_def_varn  ! define n-dimensional variable
     module procedure netcdf_def_var1d ! define 1-D variable
  end interface

! Fortran-2000 may have templates, but until then, we have to define
! a different subroutine for every dimensionality array:

  interface netcdf_put_var
     module procedure netcdf_put_var_int0  ! write a scalar integer
     module procedure netcdf_put_var_int1  ! write a whole 1-D integer array
     module procedure netcdf_put_var_int2  ! etc.
     module procedure netcdf_put_var_int3
     module procedure netcdf_put_var_int4
     module procedure netcdf_put_var_int5
     module procedure netcdf_put_var_int6
     module procedure netcdf_put_var_int7
     module procedure netcdf_put_var_double0
     module procedure netcdf_put_var_double1
     module procedure netcdf_put_var_double2
     module procedure netcdf_put_var_double3
     module procedure netcdf_put_var_double4
     module procedure netcdf_put_var_double5
     module procedure netcdf_put_var_double6
     module procedure netcdf_put_var_double7
     module procedure netcdf_put_var_text
  end interface

  interface netcdf_put_var1
     module procedure netcdf_put_var1_int1  ! write an element of a 1-D integer array
     module procedure netcdf_put_var1_inta  ! write an element of an integer array
     module procedure netcdf_put_var1_double1 ! same for real array...
     module procedure netcdf_put_var1_doublea
  end interface

  interface netcdf_put_vara
!    module procedure netcdf_put_vara_int
     module procedure netcdf_put_vara_double1_11
     module procedure netcdf_put_vara_double1_00
     module procedure netcdf_put_vara_double1_10
     module procedure netcdf_put_vara_double1_01
     module procedure netcdf_put_vara_double2
     module procedure netcdf_put_vara_double3
     module procedure netcdf_put_vara_double4
     module procedure netcdf_put_vara_double5
     module procedure netcdf_put_vara_double6
     module procedure netcdf_put_vara_double7
!    module procedure netcdf_put_vara_text
  end interface

  interface netcdf_put_att
     module procedure netcdf_put_att_text
  end interface

  interface netcdf_get_var
     module procedure netcdf_get_var_int0
     module procedure netcdf_get_var_int1
     module procedure netcdf_get_var_int2
     module procedure netcdf_get_var_int3
     module procedure netcdf_get_var_int4
     module procedure netcdf_get_var_int5
     module procedure netcdf_get_var_int6
     module procedure netcdf_get_var_int7
     module procedure netcdf_get_var_double0
     module procedure netcdf_get_var_double1
     module procedure netcdf_get_var_double2
     module procedure netcdf_get_var_double3
     module procedure netcdf_get_var_double4
     module procedure netcdf_get_var_double5
     module procedure netcdf_get_var_double6
     module procedure netcdf_get_var_double7
     module procedure netcdf_get_var_text
  end interface

contains

!#######################################################################

! serial_io2 = .true. to force serial interaction with netcdf, 
! to avoid bugs in parallel netcdf.
! serial_io2 = .false. if parallel netcdf is every fixed.

  subroutine netcdf_init(serial_io2)
    logical, intent(in) :: serial_io2
    serial_io=serial_io2

    if(serial_io) then
       ! only proc0 will write
       proc_write=proc0
    else
       ! all processors will write
       proc_write=.true.
    endif

  end subroutine netcdf_init

!#######################################################################

  function netcdf_def_dim( ncid, name, length, dimid) result(status)

! has identical arguments to nf_def_dim, but adds some automatic error 
! checking. 

    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(in) :: length
    integer, intent(out) :: dimid
    integer :: status
    status = nf_noerr

!    if(proc0) write(*,*) 'defining NF dimension, length: ', name, length 

    if(serial_io) then
       if(proc0) status = nf_def_dim( ncid, name, length, dimid)
!       call broadcast(dimid)
    else
       status = nf_def_dim( ncid, name, length, dimid)
    endif

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) ' in routine netcdf_def_dim(ncid,name,length,dimid) ='
       write(*,*) ncid,name,length, dimid
    endif

  end function netcdf_def_dim

!#######################################################################

  function netcdf_def_varn( ncid, name, xtype, nvdims, vdims, varid) &
       result(status)

! netcdf_def_varn has identical arguments to nf_def_var, but
! it adds some automatic error checking.

    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(in) :: xtype
    integer, intent(in) :: nvdims
    integer, intent(in), dimension(:) :: vdims
    integer, intent(out) :: varid

    integer :: status
    status = nf_noerr

!    if(proc0) write(*,*) 'defining NF variable, type, dims, vdims: ', &
!              name, xtype, nvdims, vdims

    if(serial_io) then
       if(proc0) status = nf_def_var( ncid, name, xtype, nvdims, vdims, varid)
!       call broadcast(varid)
    else
       status = nf_def_var( ncid, name, xtype, nvdims, vdims, varid)
    endif

    if(debug) then
       write(*,*) &
            'netcdf_def_varn( ncid, name, xtype, nvdims, vdims, varid)', &
            ncid, name, xtype, nvdims, vdims, varid
    endif

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_def_vara(ncid, name, xtype, nvdims, vdims, varid)='
       write(*,*) ncid, name, xtype, nvdims, vdims, varid
    endif

  end function netcdf_def_varn

!#######################################################################
  function netcdf_def_var1d( ncid, name, xtype, nvdims, vdims, varid) &
       result(status)

! handles the case when vdims is a scalar and not an array
! (which can happen if the variable being defined is a scalar or a
! a 1-d vector).

    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(in) :: xtype
    integer, intent(in) :: nvdims
    integer, intent(in) :: vdims
    integer, intent(out) :: varid

    integer :: status
    status = nf_noerr

!    if(proc0) write(*,*) &
!              'defining NF scalar or 1-D vector, type, dims, vdims: ', &
!              name, xtype, nvdims, vdims

    if(serial_io) then
       if(proc0) status = nf_def_var( ncid, name, xtype, nvdims, vdims, varid)
!       call broadcast(varid)
    else
       status = nf_def_var( ncid, name, xtype, nvdims, vdims, varid)
    endif

    if(debug) then
       write(*,*) &
            'netcdf_def_var1d( ncid, name, xtype, nvdims, vdims, varid)',&
            ncid, name, xtype, nvdims, vdims, varid
    endif

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_def_var1d(ncid, name, xtype, nvdims, vdims, varid)='
       write(*,*) ncid, name, xtype, nvdims, vdims, varid
    endif

  end function netcdf_def_var1d

!#######################################################################

  function netcdf_put_att_text( ncid, varid, name, text) result(status)

! has identical arguments to nf_put_att_text, except the "len" argument
! has been removed, since in almost all cases it easier to determine it
! from the length of the character string "text" then have the calling
! program figure it out.

    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: name
    character*(*), intent(in) :: text

    integer :: len_text
    integer :: status
    status = nf_noerr

    len_text=len(text)

!    if(proc0) write(*,*) 'defining NF attribute, varid, name, len, text', &
!         varid, name, len, text

    if(proc0) status = nf_put_att_text(  ncid, varid, name, len_text, text)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_att_text( ncid, varid, name, len, text)='
       write(*,*) ncid, varid, name, len_text, text
    endif

  end function netcdf_put_att_text

!#######################################################################
!#######################################################################
!#######################################################################
!
! Specific routines for generic netcdf_put_var (integer)
!

!#######################################################################

  function netcdf_put_var_int0 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int0

!#######################################################################

  function netcdf_put_var_int1 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int1

!#######################################################################

  function netcdf_put_var_int2 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int2

!#######################################################################

  function netcdf_put_var_int3 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int3

!#######################################################################

  function netcdf_put_var_int4 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int4

!#######################################################################

  function netcdf_put_var_int5 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int5

!#######################################################################

  function netcdf_put_var_int6 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int6

!#######################################################################

  function netcdf_put_var_int7 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:,:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_int7

!#######################################################################
!#######################################################################
!#######################################################################
!
! Specific routines for generic netcdf_put_var (double)
!

!#######################################################################

  function netcdf_put_var_double0 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in) :: var
    integer :: status
    status = nf_noerr

    if(debug) then
       write(*,*) 'netcdf_put_var_double0(ncid, varid, var) =', &
            ncid, varid, var
    endif

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double0

!#######################################################################

  function netcdf_put_var_double1 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double1

!#######################################################################

  function netcdf_put_var_double2 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double2

!#######################################################################

  function netcdf_put_var_double3 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double3

!#######################################################################

  function netcdf_put_var_double4 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double4

!#######################################################################

  function netcdf_put_var_double5 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double5

!#######################################################################

  function netcdf_put_var_double6 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double6

!#######################################################################

  function netcdf_put_var_double7 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(in), dimension(:,:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_double( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_double7

!#######################################################################

  function netcdf_put_var_text ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    character*(*), intent(in) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_text(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var_text( ncid, var_id, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_put_var_text

!#######################################################################

  function netcdf_inq_dimid ( ncid, name, dimid) result(status)

    implicit none
    integer, intent(in)  :: ncid
    character*(*), intent(in) :: name
    integer, intent(out) :: dimid
    integer :: status
    status = nf_noerr

    status = nf_inq_dimid(ncid, name, dimid)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_inq_dimid( ncid, name, dimid) ='
       write(*,*) ncid, name, dimid
       write(*,*) 'dimension ',name,' does not exist in the netcdf file'
       write(*,*) 'or problems reading the file'
       dimid=-1
    endif

  end function netcdf_inq_dimid

!#######################################################################

  function netcdf_inq_varid ( ncid, name, varid) result(status)

    implicit none
    integer, intent(in)  :: ncid
    character*(*), intent(in) :: name
    integer, intent(out) :: varid
    integer :: status
    status = nf_noerr

    status = nf_inq_varid(ncid, name, varid)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_inq_varid( ncid, name, varid) ='
       write(*,*) ncid, name, varid
       write(*,*) 'variable ',name,' does not exist in the netcdf file'
       write(*,*) 'or problems reading the file'
       varid=-1
    endif

  end function netcdf_inq_varid

!#######################################################################
!#######################################################################
!#######################################################################
!
! Specific routines for generic netcdf_get_var (integer)
!

  function netcdf_get_var_int0 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int0

!#######################################################################

  function netcdf_get_var_int1 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int1

!#######################################################################

  function netcdf_get_var_int2 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int2

!#######################################################################

  function netcdf_get_var_int3 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int3

!#######################################################################

  function netcdf_get_var_int4 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int4

!#######################################################################

  function netcdf_get_var_int5 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int5

!#######################################################################

  function netcdf_get_var_int6 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int6

!#######################################################################

  function netcdf_get_var_int7 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(out), dimension(:,:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_int(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_int( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_int7

!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!
! Specific routines for generic netcdf_get_var  (double)
!

  function netcdf_get_var_double0 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double0

!#######################################################################

  function netcdf_get_var_double1 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double1

!#######################################################################

  function netcdf_get_var_double2 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double2

!#######################################################################

  function netcdf_get_var_double3 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double3

!#######################################################################

  function netcdf_get_var_double4 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double4

!#######################################################################

  function netcdf_get_var_double5 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double5

!#######################################################################

  function netcdf_get_var_double6 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double6

!#######################################################################

  function netcdf_get_var_double7 ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    real(dp), intent(out), dimension(:,:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    status = nf_get_var_double(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_double( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_double7

!#######################################################################

  function netcdf_get_var_text ( ncid, varid, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    character*(*), intent(out) :: var
    integer :: status

    integer ndims, dimid, len_var

    status = nf_noerr

    status = nf_inq_varndims( ncid, varid, ndims)
    if(status /= nf_noerr) then
       write(*,*) 'ERROR in netcdf_get_var_text', nf_strerror(status)
    endif
    if(ndims .gt. 1) then
       write(*,*) 'ERROR in netcdf_get_var_text'
       write(*,*) 'ndims for text object should be 1'
    endif
    status = nf_inq_vardimid( ncid, varid, dimid)
    

    status = nf_inq_dimlen( ncid, dimid, len_var)
    if(status /= nf_noerr) then
       write(*,*) 'ERROR in netcdf_get_var_text', nf_strerror(status)
    endif

    if(len_var .gt. len(var)) then
       write(*,*) '***ERROR in netcdf_get_var_text, the variable passed'
       write(*,*) ' to this subroutine isn''t long enough to read the'
       write(*,*) ' netcdf string.'
       stop
    endif

    status = nf_get_var_text(ncid,varid,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_get_var_text( ncid, varid, var) ='
       write(*,*) ncid, varid, var
    endif

  end function netcdf_get_var_text

!#######################################################################
!#######################################################################
!#######################################################################
!
! Specific routines for generic netcdf_put_vara
!

!#######################################################################

  function netcdf_put_vara_double1_11 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count
    real(dp), intent(in), dimension(:) :: var

    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double1_11

!#######################################################################

  function netcdf_put_vara_double1_00 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in) :: start, count
    real(dp), intent(in), dimension(:) :: var

    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double1_00

!#######################################################################

  function netcdf_put_vara_double1_10 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:) :: start
    integer, intent(in) :: count
    real(dp), intent(in), dimension(:) :: var

    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double1_10

!#######################################################################

  function netcdf_put_vara_double1_01 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in) :: start
    integer, intent(in), dimension(:) :: count
    real(dp), intent(in), dimension(:) :: var

    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double1_01

!#######################################################################

  function netcdf_put_vara_double2 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count

    real(dp), intent(in), dimension(:,:) :: var
    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double2

!#######################################################################

  function netcdf_put_vara_double3 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count

    real(dp), intent(in), dimension(:,:,:) :: var
    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double3

!#######################################################################

  function netcdf_put_vara_double4 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count

    real(dp), intent(in), dimension(:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double4

!#######################################################################

  function netcdf_put_vara_double5 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count

    real(dp), intent(in), dimension(:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double5

!#######################################################################

  function netcdf_put_vara_double6 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count

    real(dp), intent(in), dimension(:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double6

!#######################################################################

  function netcdf_put_vara_double7 ( ncid, varid, start, count, var) &
       result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, dimension(:), intent(in) :: start, count

    real(dp), intent(in), dimension(:,:,:,:,:,:,:) :: var
    integer :: status
    status = nf_noerr

    if(proc_write) status = nf_put_vara_double( ncid, varid, start, count, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_vara_double( ncid, varid, start, count, var) ='
       write(*,*) ' (excluding var)'
       write(*,*) ncid, varid, start, count
    endif

  end function netcdf_put_vara_double7

!#######################################################################

!
! Specific routines for generic netcdf_put_var1 (integer)
!

!#######################################################################

  function netcdf_put_var1_int1 ( ncid, varid, index, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid, index, var
    integer :: status
    status = nf_noerr

    status = nf_put_var1_int(ncid, varid, index, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var1_int1( ncid, varid, index, var) ='
       write(*,*) ncid, varid, index, var
    endif

  end function netcdf_put_var1_int1

!#######################################################################

  function netcdf_put_var1_inta ( ncid, varid, index, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:) :: index
    integer, intent(in) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var1_int(ncid, varid, index, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_va1r_inta( ncid, varid, index, var) ='
       write(*,*) ncid, varid, index, var
    endif

  end function netcdf_put_var1_inta

!#######################################################################
!#######################################################################
!#######################################################################
!
! Specific routines for generic netcdf_put_var1 (double)
!

!#######################################################################

  function netcdf_put_var1_double1 ( ncid, varid, index, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid, index
    real(dp), intent(in) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var1_double(ncid, varid, index, var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var1_double1( ncid, var_id, index, var) ='
       write(*,*) ncid, varid, index, var
    endif

  end function netcdf_put_var1_double1

!#######################################################################

  function netcdf_put_var1_doublea ( ncid, varid, index, var) result(status)

    implicit none
    integer, intent(in)  :: ncid, varid
    integer, intent(in), dimension(:)  :: index
    real(dp), intent(in) :: var
    integer :: status
    status = nf_noerr

    status = nf_put_var_double(ncid,varid,index,var)

    if(status /= nf_noerr) then
       write(*,*) '***ERROR in NetCDF: ', nf_strerror(status)
       write(*,*) &
       ' in routine netcdf_put_var1_doublea( ncid, var_id, index, var) ='
       write(*,*) ncid, varid, index, var
    endif

  end function netcdf_put_var1_doublea


end module netcdf_mod
