# include "define.inc"

module fft_work

  use constants, only: kind_id

  implicit none

  public :: fft_type, delete_fft
  public :: init_ccfftw, init_crfftw, init_rcfftw, init_z
  public :: FFTW_FORWARD, FFTW_BACKWARD

  private

  interface init_crfftw
     module procedure init_crfftw_1d
     module procedure init_crfftw_2d
  end interface

  interface init_rcfftw
     module procedure init_rcfftw_1d
     module procedure init_rcfftw_2d
  end interface

  type :: fft_type
! TT>
!     integer :: n, plan, is, type
     integer :: n, is, type
     integer (kind_id) :: plan
! <TT
     real :: scale
  end type fft_type

  ! parameters defined in fftw_f77.i
  integer, parameter :: fftw_estimate   =  0
  integer, parameter :: fftw_measure    =  1
  integer, parameter :: fftw_in_place   =  8
  integer, parameter :: fftw_use_wisdom = 16
  integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  integer, parameter :: FFTW_REAL_TO_COMPLEX=FFTW_FORWARD
  integer, parameter :: FFTW_COMPLEX_TO_REAL=FFTW_BACKWARD

contains

  subroutine init_z (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
    
    j = fftw_measure + fftw_use_wisdom
# if FFT == _FFTW_
    call fftw_f77_create_plan(fft%plan,n,is,j)
# endif

  end subroutine init_z

  subroutine init_ccfftw (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
    
    j = fftw_in_place + fftw_measure + fftw_use_wisdom
# if FFT == _FFTW_
    call fftw_f77_create_plan(fft%plan,n,is,j)
# endif

  end subroutine init_ccfftw

  subroutine init_rcfftw_1d (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0

    j = fftw_measure + fftw_use_wisdom
# if FFT == _FFTW_
    call rfftwnd_f77_create_plan(fft%plan,1,N,is,j)
# endif

  end subroutine init_rcfftw_1d
  
  subroutine init_crfftw_1d (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0

    j = fftw_measure + fftw_use_wisdom
# if FFT == _FFTW_
    call rfftwnd_f77_create_plan(fft%plan,1,N,is,j)
# endif

  end subroutine init_crfftw_1d

  subroutine init_rcfftw_2d (fft, is, m, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

    j = fftw_measure + fftw_use_wisdom
# if FFT == _FFTW_
    call rfftw2d_f77_create_plan(fft%plan,m,n,is,j)
# endif

  end subroutine init_rcfftw_2d
  
  subroutine init_crfftw_2d (fft, is, m, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

    j = fftw_measure + fftw_use_wisdom
# if FFT == _FFTW_
    call rfftw2d_f77_create_plan(fft%plan,m,n,is,j)
# endif

  end subroutine init_crfftw_2d

  subroutine delete_fft(fft)
    
    type (fft_type), intent (in out) :: fft

# if FFT == _FFTW_
    if (fft%type == 1) then
       call fftw_f77_destroy_plan(fft%plan)
    else
       call rfftw_f77_destroy_plan(fft%plan)
    end if
# endif

  end subroutine delete_fft

end module fft_work
