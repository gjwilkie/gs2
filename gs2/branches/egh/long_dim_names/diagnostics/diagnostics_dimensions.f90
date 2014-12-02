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
  
  public :: dim_string, diagnostics_dimension_list_type

  integer, parameter :: string_len=256, dim_string_len=20

  !> This type is used to describe a single dimension
  type diagnostics_dimension_type
     private
     logical :: is_unlimited 
     logical :: initialised=.false.
     integer :: length !Size of dimension
     character(len=dim_string_len) :: dim_name !String representing the dimension 
     character(len=string_len) :: description !String to use as a description
     character(len=string_len) :: units !String description of units
   contains
     private
     procedure, public :: init => dimension_type_init
     procedure, public :: add_to_file => dimension_add_to_file
     procedure, public :: get_dim_name => dimension_get_dim_name
  end type diagnostics_dimension_type

  !> This type is a container for all the specific dimension instances
  type diagnostics_dimension_list_type
     type(diagnostics_dimension_type) :: kx
     type(diagnostics_dimension_type) :: ky
     type(diagnostics_dimension_type) :: theta
     type(diagnostics_dimension_type) :: theta_ext
     type(diagnostics_dimension_type) :: xx
     type(diagnostics_dimension_type) :: yy
     type(diagnostics_dimension_type) :: energy
     type(diagnostics_dimension_type) :: lambda
     type(diagnostics_dimension_type) :: species
     type(diagnostics_dimension_type) :: vpar
     type(diagnostics_dimension_type) :: time
     type(diagnostics_dimension_type) :: ri
     type(diagnostics_dimension_type) :: generic_2
     type(diagnostics_dimension_type) :: generic_3
     type(diagnostics_dimension_type) :: generic_4
     type(diagnostics_dimension_type) :: generic_5
     type(diagnostics_dimension_type) :: input_file_dim
  end type diagnostics_dimension_list_type

  interface dim_string
     module procedure :: make_dim_string
     module procedure :: make_dim_string_arr
  end interface dim_string

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
  subroutine dimension_add_to_file(self,sfile)
    use simpledataio, only: add_dimension, sdatio_file
    implicit none
    class(diagnostics_dimension_type), intent(in) :: self
    type(sdatio_file), intent(in) :: sfile

    if(.not.self%initialised) return

    call add_dimension(sfile,trim(self%dim_name),&
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
!<DD>Comment following as some old compilers don't support allocatable strings
!at some point in the future we can reinstate this.
!    character(len=:), allocatable :: make_dim_string
!This is the replacement
    character(len=dim_string_len) :: make_dim_string
!<DD>Can't do following with older compilers
!    allocate(character(len=len(trim(dim%get_dim_name())))::make_dim_string)
    make_dim_string=trim(dim%get_dim_name())
  end function make_dim_string

  !>Returns a string representing the dimensions, to pass to
  !!the file/io routines when creating variables etc.
  function make_dim_string_arr(dims)
    type(diagnostics_dimension_type), dimension(:), intent(in) :: dims
!<DD>Comment following as some old compilers don't support allocatable strings
!at some point in the future we can reinstate this.
!    character(len=:), allocatable :: make_dim_string_arr, tmp
!This is the replacement
    character(len=dim_string_len*10) :: make_dim_string_arr, tmp
    character(len=*), parameter :: join_char=','
    integer :: ndim, i

    !Count dimensions
    ndim=size(dims)

!    !Make temporary storage
!<DD>Can't do following with older compilers
!    allocate(character(len=dim_string_len*ndim+len(join_char)*(ndim-1))::tmp)
    tmp=''
    if(.not.(ndim.lt.1))then
       tmp=trim(dims(1)%get_dim_name())
       if(.not.(ndim.eq.1)) then
          do i=2,ndim
             tmp=trim(tmp)//join_char//trim(dims(i)%get_dim_name())
          enddo
       endif
    endif
!<DD>Can't do following with older compilers
!    allocate(character(len=len(trim(tmp)))::make_dim_string_arr)
    make_dim_string_arr=trim(tmp)
    !write (*,*) make_dim_string_arr
!<DD>Can't do following with older compilers
!    deallocate(tmp)
  end function make_dim_string_arr
end module diagnostics_dimensions
