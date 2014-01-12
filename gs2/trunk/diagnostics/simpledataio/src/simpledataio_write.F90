
module simpledataio_write

interface write_variable
  module procedure write_variable_real_0
  module procedure write_variable_real_1
  module procedure write_variable_real_2
  module procedure write_variable_real_3
  module procedure write_variable_real_4
  module procedure write_variable_real_5
  module procedure write_variable_real_6
  module procedure write_variable_double_precision_0
  module procedure write_variable_double_precision_1
  module procedure write_variable_double_precision_2
  module procedure write_variable_double_precision_3
  module procedure write_variable_double_precision_4
  module procedure write_variable_double_precision_5
  module procedure write_variable_double_precision_6
  module procedure write_variable_integer_0
  module procedure write_variable_integer_1
  module procedure write_variable_integer_2
  module procedure write_variable_integer_3
  module procedure write_variable_integer_4
  module procedure write_variable_integer_5
  module procedure write_variable_integer_6
  module procedure write_variable_character_0
  module procedure write_variable_character_1
  module procedure write_variable_character_2
  module procedure write_variable_character_3
  module procedure write_variable_character_4
  module procedure write_variable_character_5
  module procedure write_variable_character_6
  module procedure write_variable_complex_0
  module procedure write_variable_complex_1
  module procedure write_variable_complex_2
  module procedure write_variable_complex_3
  module procedure write_variable_complex_4
  module procedure write_variable_complex_5
  module procedure write_variable_complex_6
  module procedure write_variable_complex_16_0
  module procedure write_variable_complex_16_1
  module procedure write_variable_complex_16_2
  module procedure write_variable_complex_16_3
  module procedure write_variable_complex_16_4
  module procedure write_variable_complex_16_5
  module procedure write_variable_complex_16_6
end interface write_variable

contains

 subroutine write_variable_real_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in) :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+0))
   allocate(counts(n+0))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   (/val/), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_0

 subroutine write_variable_real_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in), dimension(:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+1))
   allocate(counts(n+1))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_1

 subroutine write_variable_real_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in), dimension(:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+2))
   allocate(counts(n+2))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_2

 subroutine write_variable_real_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in), dimension(:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+3))
   allocate(counts(n+3))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_3

 subroutine write_variable_real_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in), dimension(:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+4))
   allocate(counts(n+4))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_4

 subroutine write_variable_real_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in), dimension(:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+5))
   allocate(counts(n+5))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_5

 subroutine write_variable_real_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   real, intent(in), dimension(:,:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+6))
   allocate(counts(n+6))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_real_6

 subroutine write_variable_double_precision_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in) :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+0))
   allocate(counts(n+0))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   (/val/), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_0

 subroutine write_variable_double_precision_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in), dimension(:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+1))
   allocate(counts(n+1))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_1

 subroutine write_variable_double_precision_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in), dimension(:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+2))
   allocate(counts(n+2))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_2

 subroutine write_variable_double_precision_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in), dimension(:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+3))
   allocate(counts(n+3))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_3

 subroutine write_variable_double_precision_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in), dimension(:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+4))
   allocate(counts(n+4))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_4

 subroutine write_variable_double_precision_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in), dimension(:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+5))
   allocate(counts(n+5))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_5

 subroutine write_variable_double_precision_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   double precision, intent(in), dimension(:,:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+6))
   allocate(counts(n+6))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_double_precision_6

 subroutine write_variable_integer_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in) :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+0))
   allocate(counts(n+0))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   (/val/), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_0

 subroutine write_variable_integer_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in), dimension(:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+1))
   allocate(counts(n+1))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_1

 subroutine write_variable_integer_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in), dimension(:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+2))
   allocate(counts(n+2))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_2

 subroutine write_variable_integer_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in), dimension(:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+3))
   allocate(counts(n+3))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_3

 subroutine write_variable_integer_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in), dimension(:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+4))
   allocate(counts(n+4))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_4

 subroutine write_variable_integer_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in), dimension(:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+5))
   allocate(counts(n+5))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_5

 subroutine write_variable_integer_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   integer, intent(in), dimension(:,:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+6))
   allocate(counts(n+6))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_integer_6

 subroutine write_variable_character_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in) :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+0))
   allocate(counts(n+0))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   (/val/), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_0

 subroutine write_variable_character_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in), dimension(:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+1))
   allocate(counts(n+1))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_1

 subroutine write_variable_character_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in), dimension(:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+2))
   allocate(counts(n+2))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_2

 subroutine write_variable_character_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in), dimension(:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+3))
   allocate(counts(n+3))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_3

 subroutine write_variable_character_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in), dimension(:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+4))
   allocate(counts(n+4))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_4

 subroutine write_variable_character_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in), dimension(:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+5))
   allocate(counts(n+5))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_5

 subroutine write_variable_character_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   character, intent(in), dimension(:,:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+6))
   allocate(counts(n+6))
	   
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   val(starts(1):,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_character_6

 subroutine write_variable_complex_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in) :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+0+1))
   allocate(counts(n+0+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1) = real(val)
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval, start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1) = aimag(val)
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval, start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_0

 subroutine write_variable_complex_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in), dimension(:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1, size(val, 1)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+1+1))
   allocate(counts(n+1+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:) = real(val(:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:) = aimag(val(:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_1

 subroutine write_variable_complex_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in), dimension(:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1, size(val, 1), size(val, 2)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+2+1))
   allocate(counts(n+2+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:) = real(val(:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:) = aimag(val(:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_2

 subroutine write_variable_complex_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in), dimension(:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1, size(val, 1), size(val, 2), size(val, 3)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+3+1))
   allocate(counts(n+3+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:) = real(val(:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:) = aimag(val(:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_3

 subroutine write_variable_complex_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in), dimension(:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1, size(val, 1), size(val, 2), size(val, 3), size(val, 4)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+4+1))
   allocate(counts(n+4+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:,:) = real(val(:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:,:) = aimag(val(:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_4

 subroutine write_variable_complex_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in), dimension(:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+5+1))
   allocate(counts(n+5+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:,:,:) = real(val(:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:,:,:) = aimag(val(:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_5

 subroutine write_variable_complex_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex, intent(in), dimension(:,:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 real, dimension(1, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5), size(val, 6)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+6+1))
   allocate(counts(n+6+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:,:,:,:) = real(val(:,:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):,starts(7):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:,:,:,:) = aimag(val(:,:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):,starts(7):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_6

 subroutine write_variable_complex_16_0(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in) :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+0+1))
   allocate(counts(n+0+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1) = real(val)
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval, start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1) = aimag(val)
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval, start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_0

 subroutine write_variable_complex_16_1(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in), dimension(:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1, size(val, 1)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+1+1))
   allocate(counts(n+1+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:) = real(val(:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:) = aimag(val(:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_1

 subroutine write_variable_complex_16_2(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in), dimension(:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1, size(val, 1), size(val, 2)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+2+1))
   allocate(counts(n+2+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:) = real(val(:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:) = aimag(val(:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_2

 subroutine write_variable_complex_16_3(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in), dimension(:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1, size(val, 1), size(val, 2), size(val, 3)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+3+1))
   allocate(counts(n+3+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:) = real(val(:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:) = aimag(val(:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_3

 subroutine write_variable_complex_16_4(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in), dimension(:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1, size(val, 1), size(val, 2), size(val, 3), size(val, 4)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+4+1))
   allocate(counts(n+4+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:,:) = real(val(:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:,:) = aimag(val(:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_4

 subroutine write_variable_complex_16_5(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in), dimension(:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+5+1))
   allocate(counts(n+5+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:,:,:) = real(val(:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:,:,:) = aimag(val(:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_5

 subroutine write_variable_complex_16_6(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf
#endif
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   complex*16, intent(in), dimension(:,:,:,:,:,:)  :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 double precision, dimension(1, size(val, 1), size(val, 2), size(val, 3), size(val, 4), size(val, 5), size(val, 6)) :: realval

#ifdef FORTRAN_NETCDF
   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+6+1))
   allocate(counts(n+6+1))
	 
	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 1)
	 realval(1,:,:,:,:,:,:) = real(val(:,:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):,starts(7):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)


	 call set_count(sfile, variable_name, "r", 1)
	 call set_start(sfile, variable_name, "r", 2)
	 realval(1,:,:,:,:,:,:) = aimag(val(:,:,:,:,:,:))
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   realval(1:,starts(2):,starts(3):,starts(4):,starts(5):,starts(6):,starts(7):), start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)

#endif

 end subroutine write_variable_complex_16_6


end module simpledataio_write

