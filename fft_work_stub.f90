module fft_work

  implicit none

  public :: fft_type, delete_fft
  public :: init_ccfftw, init_crfftw, init_rcfftw
  
  type :: fft_type
     integer :: n, plan, is, type
     real :: scale
  end type fft_type

!  external fftw_f77_create_plan,fftw_f77_one
!  external fftw_f77_destroy_plan

!  external rfftw_f77_create_plan,rfftw_f77_one
!  external rfftw_f77_destroy_plan
  
  integer fftw_in_place
  parameter (fftw_in_place=8)

  integer fftw_estimate,fftw_measure
  parameter (fftw_estimate=0,fftw_measure=1)
  
  integer, parameter, public :: ccfft_table0 = 0
  integer, parameter, public :: ccfft_table1 = 0
  integer, parameter, public :: ccfft_work0 = 0
  integer, parameter, public :: ccfft_work1 = 0

  integer, parameter, public :: csfft_table0 = 0
  integer, parameter, public :: csfft_table1 = 0
  integer, parameter, public :: csfft_work0 = 0
  integer, parameter, public :: csfft_work1 = 0

  private

contains

  subroutine init_ccfftw (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
    
    j = fftw_in_place + fftw_measure
!    call fftw_f77_create_plan(fft%plan,n,is,j)

  end subroutine init_ccfftw

  subroutine init_rcfftw (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0

!    call rfftw_f77_create_plan(fft%plan,N,is,FFTW_MEASURE)

  end subroutine init_rcfftw
  
  subroutine init_crfftw (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0

!    call rfftw_f77_create_plan(fft%plan,N,is,FFTW_MEASURE)


  end subroutine init_crfftw

  subroutine delete_fft(fft)
    
    type (fft_type), intent (in out) :: fft

    if (fft%type == 1) then
!       call fftw_f77_destroy_plan(fft%plan)
    else
!       call rfftw_f77_destroy_plan(fft%plan)
    end if

  end subroutine delete_fft

end module fft_work
