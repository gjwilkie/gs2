module fft_work

  implicit none

  public :: fft_type, delete_fft
  public :: init_ccfftw, init_crfftw, init_rcfftw, init_z
  
  interface init_crfftw
     module procedure init_crfftw_1d
     module procedure init_crfftw_2d
  end interface

  interface init_rcfftw
     module procedure init_rcfftw_1d
     module procedure init_rcfftw_2d
  end interface

  type :: fft_type
     integer*8 :: plan
     integer :: n, is, type
     real :: scale
  end type fft_type

  integer, parameter :: fftw_estimate   =  0
  integer, parameter :: fftw_measure    =  1
  integer, parameter :: fftw_in_place   =  8
  integer, parameter :: fftw_use_wisdom = 16

  private

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
    call fftw_f77_create_plan(fft%plan,n,is,j)

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
    call fftw_f77_create_plan(fft%plan,n,is,j)

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
    call rfftwnd_f77_create_plan(fft%plan,1,N,is,j)

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
    call rfftwnd_f77_create_plan(fft%plan,1,N,is,j)

  end subroutine init_crfftw_1d

  subroutine init_rcfftw_2d (fft, is, m, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

    j = fftw_measure + fftw_use_wisdom
    call rfftw2d_f77_create_plan(fft%plan,m,n,is,j)

  end subroutine init_rcfftw_2d
  
  subroutine init_crfftw_2d (fft, is, m, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

    j = fftw_measure + fftw_use_wisdom
    call rfftw2d_f77_create_plan(fft%plan,m,n,is,j)

  end subroutine init_crfftw_2d

  subroutine delete_fft(fft)
    
    type (fft_type), intent (in out) :: fft

    if (fft%type == 1) then
       call fftw_f77_destroy_plan(fft%plan)
    else
       call rfftw_f77_destroy_plan(fft%plan)
    end if

  end subroutine delete_fft

end module fft_work
