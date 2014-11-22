!> This module provides a dimension class (derived type)
!! The idea is to make it easier to support different
!! output backends at some point by putting backend
!! specific setup in here.
!! 
!! This module is also the home to the specific instances
!! of dimensions used in the diagnostics.
module diagnostics_dimensions
  implicit none

  private
  
!  public :: diagnostics_dimension_type
  public :: kx_dim, ky_dim, theta_dim, e_dim
  public :: l_dim, spec_dim, ri_dim, time_dim
  public :: vpar_dim, theta_ext_dim
  public :: generic_2_dim, generic_3_dim
  public :: generic_4_dim, generic_5_dim
  public :: xx_dim, yy_dim

  public :: dim_string

  integer, parameter :: string_len=256

  type diagnostics_dimension_type
     private
     logical :: is_unlimited 
     logical :: initialised=.false.
     integer :: length !Size of dimension
     character(len=20) :: dim_name !String representing the dimension 
     character(len=string_len) :: description !String to use as a description
     character(len=string_len) :: units !String description of units
   contains
     private
     procedure, public :: init => dimension_type_init
     procedure, public :: add_to_file => dimension_add_to_file
     procedure, public :: get_dim_name => dimension_get_dim_name
  end type diagnostics_dimension_type

  interface dim_string
     module procedure :: make_dim_string
     module procedure :: make_dim_string_arr
  end interface dim_string

  !Specific dimensions used throughout the diagnostics
  type(diagnostics_dimension_type) :: kx_dim, ky_dim, theta_dim, e_dim
  type(diagnostics_dimension_type) :: l_dim, spec_dim, ri_dim, time_dim
  type(diagnostics_dimension_type) :: vpar_dim, theta_ext_dim
  type(diagnostics_dimension_type) :: generic_2_dim, generic_3_dim
  type(diagnostics_dimension_type) :: generic_4_dim, generic_5_dim
  type(diagnostics_dimension_type) :: xx_dim, yy_dim
contains
  !/////////////////////////////
  !// TYPE BOUND PROCEDURES
  !/////////////////////////////

  !> Populate the dimension type
  subroutine dimension_type_init(self,dim_name_in,length_in,description_in,units_in,is_unlimited_in)
    use simpledataio, only: SDATIO_UNLIMITED
    implicit none
    class(diagnostics_dimension_type), intent(inout) :: self
    character(len=*), intent(in) :: dim_name_in
    character(len=*), intent(in) :: description_in, units_in
    logical, intent(in), optional :: is_unlimited_in
    integer, intent(in) :: length_in

    if(self%initialised) return

    self%dim_name=trim(dim_name_in)
    self%description=trim(description_in)
    self%units=trim(units_in)
    self%length=length_in
    if(present(is_unlimited_in)) then
       if(is_unlimited_in)then
          self%length=SDATIO_UNLIMITED
       endif
    endif
    self%is_unlimited=(self%length.eq.SDATIO_UNLIMITED)
    self%initialised=.true.
  end subroutine dimension_type_init
  
  !>Attach the dimension to file
  subroutine dimension_add_to_file(self,gnostics)
    use diagnostics_config, only: diagnostics_type
    use simpledataio, only: add_dimension
    implicit none
    class(diagnostics_dimension_type), intent(in) :: self
    type(diagnostics_type), intent(in) :: gnostics

    if(.not.self%initialised) return

    call add_dimension(gnostics%sfile,trim(self%dim_name),&
         self%length,trim(self%description),trim(self%units))
  end subroutine dimension_add_to_file

  !>Return the dimensions name
  function dimension_get_dim_name(self)
    implicit none
    class(diagnostics_dimension_type), intent(in) :: self
    character(len=string_len) :: dimension_get_dim_name
    if(.not.self%initialised)then
       dimension_get_dim_name=''
    else
       dimension_get_dim_name=self%dim_name
    endif
  end function dimension_get_dim_name

  !/////////////////////////////
  !// STANDARD PROCEDURES
  !/////////////////////////////

  !>Returns a string representing the dimension, to pass to
  !!the file/io routines when creating variables etc.
  function make_dim_string(dim)
    type(diagnostics_dimension_type), intent(in) :: dim
    character(len=string_len) :: make_dim_string
    make_dim_string=trim(dim%get_dim_name())
  end function make_dim_string

  !>Returns a string representing the dimensions, to pass to
  !!the file/io routines when creating variables etc.
  function make_dim_string_arr(dims)
    type(diagnostics_dimension_type), dimension(:), intent(in) :: dims
    character(len=string_len) :: make_dim_string_arr
    character(len=*), parameter :: join_char=''
    integer :: ndim, i
    ndim=size(dims)
    make_dim_string_arr=''
    if(ndim.lt.1) return
    make_dim_string_arr=trim(dims(1)%get_dim_name())
    if(ndim.eq.1) return
    do i=2,ndim
       make_dim_string_arr=trim(make_dim_string_arr)//join_char//trim(dims(i)%get_dim_name())
    enddo
  end function make_dim_string_arr

end module diagnostics_dimensions
