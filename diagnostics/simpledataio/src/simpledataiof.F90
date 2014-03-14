
!> Fortran interface to simpledataio, which requires the
!! iso_c_binding intrinsic module introduced in Fortran 2003.
!! All functions of simpledataio are provided in this module
!! with the exception of write_variable which is provided in
!! simpledataio_write
module simpledataio

use netcdf 

#ifdef ISO_C_BINDING
use iso_c_binding
#endif

integer, parameter :: SDATIO_INT= 0
integer, parameter :: SDATIO_FLOAT= 1
integer, parameter :: SDATIO_DOUBLE= 2
integer, parameter :: SDATIO_COMPLEX_DOUBLE= 3

integer, parameter :: SDATIO_UNLIMITED = NF90_UNLIMITED


#ifdef ISO_C_BINDING
type,bind(c) :: sdatio_dimension 
  type(c_ptr) :: name
	integer(c_int) :: size
	integer(c_int) :: nc_id
	integer(c_int) :: start
end type


type :: sdatio_variable 
	character, dimension(:), allocatable :: name
	integer :: nc_id
	integer :: type
  type(c_ptr) :: dimension_list
  type(c_ptr) :: dimension_ids
	integer :: type_size
  type(c_ptr) :: manual_counts
  type(c_ptr) :: manual_starts
  type(c_ptr) :: manual_offsets
end type


type, bind(c) :: sdatio_file 
	integer(c_int) :: nc_file_id
	integer(c_int):: is_parallel
	integer(c_int) :: n_dimensions
  type(c_ptr) ::  dimensions
	integer(c_int)  :: n_variables
  type(c_ptr) :: variables
	integer(c_int) :: data_written
end type

#else

type :: sdatio_file
  integer :: nc_file_id
  integer :: is_parallel
  integer :: n_dimensions  
  integer :: n_variables
  integer :: data_written
end type

#endif

!interface 
!!/* Open a new datafile for writing. fname is the name of the file 
 !!* The stuct sfile is used to store the state information
 !!* of the file.*/
!!/* Create a new dimension in the file sfile. Dimension names must
 !!* be a single letter. */
 !subroutine sdatio_add_dimension(sfile, dimension_name, dimsize, description, units)
   !import sdatio_file
   !type(sdatio_file), intent(in) :: sfile
   !character(*), intent(in) :: dimension_name
   !integer, intent(in) :: dimsize
   !character(*), intent(in) :: description, units
 !end subroutine sdatio_add_dimension


!end interface



!int sdatio_debug

!/* Open a new datafile for writing. fname is the name of the file 
 !* The stuct sfile is used to store the state information
 !* of the file.*/
contains 
  subroutine createfile(sfile, fname)
     type(sdatio_file), intent(out) :: sfile
     character(*), intent(in) :: fname
#ifdef ISO_C_BINDING
#ifdef FORTRAN_NETCDF
     interface
       subroutine sdatio_createfile(sfile, fname) bind(c, name='sdatio_createfile')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         character(c_char) :: fname(*)
       end subroutine sdatio_createfile
     end interface
     call sdatio_createfile(sfile, fname//c_null_char)
#else
     write (*,*) "module simpledataio was built without &
       & the netcdf fortran library and is non-functional. You can use &
       & the function simpledataio_functional to test &
       & this. "
     stop
#endif
#else
     write (*,*) "module simpledataio was built without &
       & ISO_C_BINDING and is non-functional. You can use &
       & the function simpledataio_functional to test &
       & this. "
     stop
#endif
   end subroutine createfile

!#ifdef PARALLEL 
  subroutine createfile_parallel(sfile, fname, comm)
     type(sdatio_file), intent(out) :: sfile
     character(*), intent(in) :: fname
     integer, intent(in) :: comm
#ifdef ISO_C_BINDING
#ifdef FORTRAN_NETCDF
     interface
       subroutine sdatio_createfile_parallel(sfile, fname, comm) bind(c, name='sdatio_createfile_parallel')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         character(c_char) :: fname(*)
         integer(c_int),value :: comm
       end subroutine sdatio_createfile_parallel
     end interface
     call sdatio_createfile_parallel(sfile, fname//c_null_char, comm)
#else
     write (*,*) "module simpledataio was built without &
       & the netcdf fortran library and is non-functional. You can use &
       & the function simpledataio_functional to test &
       & this. "
     stop
#endif
#else
     write (*,*) "module simpledataio was built without &
       & ISO_C_BINDING and is non-functional. You can use &
       & the function simpledataio_functional to test &
       & this. "
     stop
#endif
   end subroutine createfile_parallel
!#endif

!/* Create a new dimension in the file sfile. Dimension names must
 !* be a single letter. */
 subroutine add_dimension(sfile, dimension_name, dimsize, description, units)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: dimension_name
   integer, intent(in) :: dimsize
   character(*), intent(in) :: description, units
#ifdef ISO_C_BINDING
   interface
       subroutine sdatio_add_dimension(sfile, dimension_name, dimsize, description, units) bind(c, name='sdatio_add_dimension')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         character(c_char) :: dimension_name(*)
         integer(c_int), value :: dimsize
         character(c_char) :: units(*)
         character(c_char) :: description(*)
       end subroutine sdatio_add_dimension
   end interface 
   call sdatio_add_dimension(sfile, dimension_name//c_null_char, dimsize, description//c_null_char, units//c_null_char)
#endif
 end subroutine add_dimension
													 !char * dimension_name, 
													 !int size,
													 !char * description,
													 !char * units)


!/* Print out a nice list of all the dimensions defined so far*/
!void sdatio_print_dimensions(struct sdatio_file * sfile)
  subroutine print_dimensions(sfile)
   type(sdatio_file), intent(in) :: sfile
#ifdef ISO_C_BINDING
   interface
     subroutine sdatio_print_dimensions(sfile) bind(c, name='sdatio_print_dimensions')
       use iso_c_binding
       import sdatio_file
       type(sdatio_file) :: sfile
     end subroutine sdatio_print_dimensions
   end interface
   call sdatio_print_dimensions(sfile)
#endif
  end subroutine print_dimensions



!/* Close the file and free all memory associated with sfile*/
!void sdatio_close(struct sdatio_file * sfile)
  subroutine closefile(sfile)
   type(sdatio_file), intent(in) :: sfile
#ifdef ISO_C_BINDING
   interface
     subroutine sdatio_close(sfile) bind(c, name='sdatio_close')
       use iso_c_binding
       import sdatio_file
       type(sdatio_file) :: sfile
     end subroutine sdatio_close
   end interface
   call sdatio_close(sfile)
#endif
  end subroutine closefile

!/* Ensure all variables are written to disk in case of crashes*/
!void sdatio_sync(struct sdatio_file * sfile)
  subroutine syncfile(sfile)
   type(sdatio_file), intent(in) :: sfile
#ifdef ISO_C_BINDING
   interface
     subroutine sdatio_sync(sfile) bind(c, name='sdatio_sync')
       use iso_c_binding
       import sdatio_file
       type(sdatio_file) :: sfile
     end subroutine sdatio_sync
   end interface
   call sdatio_sync(sfile)
#endif
  end subroutine syncfile

!/* Define a variable in the given file. Dimension list 
 !* is a character string listing (in order) the dimension names
 !* (which are all single characters) e.g. "xyx".*/
!void sdatio_create_variable(struct sdatio_file * sfile,
														!int variable_type,
														!char * variable_name,
														!char * dimension_list,
														!char * description,
														!char * units)
 subroutine create_variable(sfile, variable_type, variable_name, dimension_list, description, units)
   implicit none
   type(sdatio_file), intent(in) :: sfile
   integer, intent(in) :: variable_type
   character(*), intent(in) :: variable_name
   character(*), intent(in) :: dimension_list
   character(*), intent(in) :: description, units
   !character, dimension(:), allocatable :: dimension_list_reversed
   character(len(dimension_list)) :: dimension_list_reversed
   integer :: i
#ifdef ISO_C_BINDING
   interface
       subroutine sdatio_create_variable(sfile, variable_type, variable_name, dimension_list, description, units) &
            bind(c, name='sdatio_create_variable')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         integer(c_int), value :: variable_type
         character(c_char) :: variable_name(*)
         character(c_char) :: dimension_list
         character(c_char) :: units(*)
         character(c_char) :: description(*)
       end subroutine sdatio_create_variable
   end interface 
   !allocate(dimension_list_reversed(len(dimension_list)))
   !if (len(dimension_list) .gt. 0) then
     do i = 1,len(dimension_list)
       dimension_list_reversed(i:i) = dimension_list(len(dimension_list)-i+1:len(dimension_list)-i+1)
     end do
   !write (*,*) 'dimension_list ', dimension_list, ' dimension_list_reversed ', dimension_list_reversed
   call sdatio_create_variable(sfile, variable_type,&
     variable_name//c_null_char, dimension_list_reversed//c_null_char, description//c_null_char, units//c_null_char)
#endif
   !else 
   !call sdatio_create_variable(sfile, variable_type,&
     !variable_name//c_null_char, dimension_list//c_null_char, description//c_null_char, units//c_null_char)
   !end if
 end subroutine create_variable



!/* Return a pointer the struct containing all the metadata of the given variable */
!struct sdatio_variable * sdatio_find_variable(struct sdatio_file * sfile, char * variable_name)
!/* Return a pointer the struct containing all the metadata of the given dimension */
!struct sdatio_dimension * sdatio_find_dimension(struct sdatio_file * sfile, char * dimension_name)
 !function find_variable(sfile,variable_name)
   !type(sdatio_file), intent(in) :: sfile
   !character(*), intent(in) :: variable_name
   !type(sdatio_variable), pointer :: find_variable
   !interface
     !type(c_ptr) function sdatio_find_variable(sfile, variable_name) &
            !bind(c, name='sdatio_find_variable' )
         !use iso_c_binding
         !import sdatio_file
         !import sdatio_variable
         !type(sdatio_file) :: sfile
         !character(c_char) :: variable_name(*)
         !!type(sdatio_variable) :: sdatio_find_variable
       !end function sdatio_find_variable
   !end interface 
   !write (*,*) 'calling'
   !call c_f_pointer(sdatio_find_variable(sfile, variable_name//c_null_char), find_variable)
 !end function find_variable

 subroutine set_offset(sfile, variable_name, dimension_name, offset)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character(*), intent(in) :: dimension_name
   integer, intent(in) :: offset
#ifdef ISO_C_BINDING
   interface
       subroutine sdatio_set_offset(sfile, variable_name, dimension_name, offset) &
            bind(c, name='sdatio_set_offset')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         integer(c_int) :: offset
         character(c_char) :: variable_name(*)
         character(c_char) :: dimension_name(*)
       end subroutine sdatio_set_offset
   end interface 
   call sdatio_set_offset(sfile, variable_name//c_null_char, &
                            dimension_name//c_null_char, offset-1)
#endif
 end subroutine set_offset

 subroutine set_count(sfile, variable_name, dimension_name, count)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character(*), intent(in) :: dimension_name
   integer, intent(in) :: count
#ifdef ISO_C_BINDING
   interface
       subroutine sdatio_set_count(sfile, variable_name, dimension_name, count) &
            bind(c, name='sdatio_set_count')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         integer(c_int) :: count
         character(c_char) :: variable_name(*)
         character(c_char) :: dimension_name(*)
       end subroutine sdatio_set_count
   end interface 
   call sdatio_set_count(sfile, variable_name//c_null_char, &
                            dimension_name//c_null_char, count)
#endif
 end subroutine set_count

 subroutine set_start(sfile, variable_name, dimension_name, start)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character(*), intent(in) :: dimension_name
   integer, intent(in) :: start
#ifdef ISO_C_BINDING
   interface
       subroutine sdatio_set_start(sfile, variable_name, dimension_name, start) &
            bind(c, name='sdatio_set_start')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         integer(c_int) :: start
         character(c_char) :: variable_name(*)
         character(c_char) :: dimension_name(*)
       end subroutine sdatio_set_start
   end interface 
   call sdatio_set_start(sfile, variable_name//c_null_char, &
                            dimension_name//c_null_char, start-1) 
#endif
       ! convert from 1-based to 0-based
                     
 end subroutine set_start

 function number_of_dimensions(sfile, variable_name)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer :: number_of_dimensions
#ifdef ISO_C_BINDING
   !integer :: c_number_of_dimensions
   interface
       function sdatio_number_of_dimensions(sfile, variable_name) &
            bind(c, name='sdatio_number_of_dimensions')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         character(c_char) :: variable_name(*)
         integer(c_int) :: sdatio_number_of_dimensions
       end function sdatio_number_of_dimensions
   end interface 
   number_of_dimensions = sdatio_number_of_dimensions(sfile, variable_name//c_null_char)
#endif
 end function number_of_dimensions

 subroutine number_of_unlimited_dimensions(sfile, variable_name, n)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(out) :: n
#ifdef ISO_C_BINDING
   interface
       subroutine sdatio_number_of_unlimited_dimensions(sfile, variable_name, n) &
            bind(c, name='sdatio_number_of_unlimited_dimensions')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         integer(c_int) :: n
         character(c_char) :: variable_name(*)
       end subroutine sdatio_number_of_unlimited_dimensions
   end interface 
   call sdatio_number_of_unlimited_dimensions(sfile, variable_name//c_null_char, n)
#endif
 end subroutine number_of_unlimited_dimensions

 subroutine netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, &
   offsets)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(out) :: fileid, varid
   integer, intent(out), dimension(:) :: starts, counts, offsets
   integer :: i,n
#ifdef ISO_C_BINDING
   type(c_ptr) :: starts_ptr, counts_ptr
   integer(c_size_t), dimension(:), allocatable, target :: starts_c, counts_c, offsets_c
   interface
       subroutine sdatio_netcdf_inputs(sfile, variable_name, fileid, varid, &
            starts_ptr, counts_ptr, offsets_ptr) &
            bind(c, name='sdatio_netcdf_inputs')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         integer(c_int) :: varid, fileid
         character(c_char) :: variable_name(*)
         type(c_ptr), value :: starts_ptr, counts_ptr, offsets_ptr
       end subroutine sdatio_netcdf_inputs
   end interface 
   allocate(starts_c(size(starts)))
   allocate(counts_c(size(counts)))
   allocate(offsets_c(size(counts)))

   call sdatio_netcdf_inputs(sfile, variable_name//c_null_char, fileid, varid, &
     c_loc(starts_c), c_loc(counts_c), c_loc(offsets_c))

   n = size(starts)
   do i = 1,n
     counts(i) = counts_c(n-i+1)
     !counts(i) = counts_c(i)
     starts(i) = starts_c(n-i+1)+1
     offsets(i) = offsets_c(n-i+1)+1
     !starts(i) = starts_c(i)+1
   end do
   
   !write (*,*) variable_name, '  sc', starts, 'c', counts, ' n', n

   deallocate(counts_c, starts_c, offsets_c)
#endif

 end subroutine netcdf_inputs



!/* Print out a nice list of all the variables defined so far*/
!void sdatio_print_variables(struct sdatio_file * sfile)
  subroutine print_variables(sfile)
   type(sdatio_file), intent(in) :: sfile
#ifdef ISO_C_BINDING
   interface
     subroutine sdatio_print_variables(sfile) bind(c, name='sdatio_print_variables')
       use iso_c_binding
       import sdatio_file
       type(sdatio_file) :: sfile
     end subroutine sdatio_print_variables
   end interface
   call sdatio_print_variables(sfile)
#endif
  end subroutine print_variables

!/* Increment the start of the specified infinite dimension */
!void sdatio_increment_start(struct sdatio_file * sfile, char * dimension_name)
  subroutine increment_start(sfile, dimension_name)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: dimension_name
#ifdef ISO_C_BINDING
   interface
     subroutine sdatio_increment_start(sfile, dimension_name) bind(c, name='sdatio_increment_start')
       use iso_c_binding
       import sdatio_file
       type(sdatio_file) :: sfile
       character(c_char) :: dimension_name
     end subroutine sdatio_increment_start
   end interface
   call sdatio_increment_start(sfile, dimension_name//c_null_char)
#endif
  end subroutine increment_start
  

 !>/* Returns .true. if the given variable has already been created, .false. otherwise */
 function variable_exists(sfile, variable_name)
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   logical :: variable_exists
#ifdef ISO_C_BINDING
   integer(c_int) :: c_variable_exists
   interface
       function sdatio_variable_exists(sfile, variable_name) &
            bind(c, name='sdatio_variable_exists')
         use iso_c_binding
         import sdatio_file
         type(sdatio_file) :: sfile
         character(c_char) :: variable_name(*)
         integer(c_int) :: sdatio_variable_exists
       end function sdatio_variable_exists
   end interface 
   c_variable_exists = sdatio_variable_exists(sfile, variable_name//c_null_char)
   if (c_variable_exists==1) then
     variable_exists = .true.
   else 
     variable_exists = .false.
   end if
#endif
 end function variable_exists 

 !> Returns true if simpledataio is functional. simpledataio is designed to
 !! present a uniform interface on any system. If it cannot be built it presents
 !! a non-functioning interface and allows the user to test for this at run time
 !! rather then causing the user's code to fail at compile time
 function simpledataio_functional()
   logical :: simpledataio_functional
   simpledataio_functional = .false.
#ifdef ISO_C_BINDING
#ifdef FORTRAN_NETCDF
   simpledataio_functional = .true.
#endif
#endif
  end function simpledataio_functional

end module simpledataio
