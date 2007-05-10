module fft_work

  implicit none

  public :: fft_type, delete_fft
  public :: init_dcft, init_dcrft, init_drcft
  
  type :: fft_type
     integer :: m, n, naux1, naux2, is
     real :: scale
     real, dimension(:), pointer :: aux1
  end type fft_type

  private

contains

  subroutine init_dcft (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    
    real, dimension(:), allocatable :: aux2
    complex :: x, y
    integer :: naux1, naux2

    call dcft_size(1, n, naux1, naux2)

    fft%m = 1
    fft%n = n
    fft%is = is
    fft%naux1 = naux1
    fft%naux2 = naux2
    fft%scale = 1./real(n)
    if (is < 0) fft%scale = 1.

    allocate (fft%aux1(naux1))
    allocate (aux2(naux2))
    
    call dcft(1, x, 1, 0, y, 1, 0, n, 1, is, &
         fft%scale, fft%aux1, naux1, aux2, naux2)

    deallocate (aux2)

  end subroutine init_dcft

  subroutine init_dcrft (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    real, dimension(:), allocatable :: aux2
    complex :: x
    real :: y

    integer :: naux1, naux2

    call dcrft_size(1, n, naux1, naux2)

    fft%m = 1
    fft%n = n
    fft%is = is
    fft%naux1 = naux1
    fft%naux2 = naux2
    fft%scale = 1./real(n)
    if (is < 0) fft%scale = 1.

    allocate (fft%aux1(naux1))
    allocate (aux2(naux2))

    call dcrft(1, x, 0, y, 0, n, 1, is, &
         fft%scale, fft%aux1, naux1, aux2, naux2)

    deallocate (aux2)

  end subroutine init_dcrft

  subroutine init_drcft (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    real, dimension(:), allocatable :: aux2
    real :: x
    complex :: y

    integer :: naux1, naux2

    call drcft_size(1, n, naux1, naux2)

    fft%m = 1
    fft%n = n
    fft%is = is
    fft%naux1 = naux1
    fft%naux2 = naux2
    fft%scale = 1./real(n)
    if (is < 0) fft%scale = 1.

    allocate (fft%aux1(naux1))
    allocate (aux2(naux2))

    call drcft(1, x, 0, y, 0, n, 1, is, &
         fft%scale, fft%aux1, naux1, aux2, naux2)

    deallocate (aux2)

  end subroutine init_drcft

  subroutine delete_fft(fft)
    
    type (fft_type), intent (in out) :: fft

    if(associated(fft%aux1)) deallocate (fft%aux1)
           
  end subroutine delete_fft

  subroutine dcft_size(m, n, nx1, nx2)

    integer, intent (in) :: m, n
    integer, intent (out) :: nx1, nx2

    complex :: x, y
    real, dimension(8) :: aux1, aux2
    external enotrm    

! set small default values
    nx1 = 8
    nx2 = 8

! initialize error table
    call einfo(0)

! make 2015 recoverable, and print no error messages.
    call errset (2015,0,-1,0,enotrm,2015)
    
! call dcft to get correct value of nx1.  Ignore nx2.
    call dcft(1, x, 1, 0, y, 1, 0, n, m, 1, 1., aux1, nx1, &
         aux2, nx2, *100)
    
100 continue

  end subroutine dcft_size

  subroutine dcrft_size(m, n, nx1, nx2)

    integer, intent (in) :: m, n
    integer, intent (out) :: nx1, nx2
    
    complex :: x
    real :: y
    real, dimension(14) :: aux1, aux2
    external enotrm
    
! I have not really done the case of m > 1 yet:

    if (m /= 1) then
       write(*,*) 'm > 1 not allowed.  You have m = ',m
       stop
    end if

! set small default values
    nx1 = 14
    nx2 = 14
    
! initialize error table
    call einfo(0)

! make 2015 recoverable, and print no error messages.
    call errset (2015,0,-1,0,enotrm,2015)
    
! call dcft to get correct value of nx1.  Ignore nx2.
    call dcrft(1, x, 0, y, 0, n, m, 1, 1., aux1, nx1, &
         aux2, nx2, *100)
    
100 continue
    
  end subroutine dcrft_size

  subroutine drcft_size(m, n, nx1, nx2)

    integer, intent (in) :: m, n
    integer, intent (out) :: nx1, nx2
    
    real :: x
    complex :: y
    real, dimension(15) :: aux1, aux2
    external enotrm
    
! I have not really done the case of m > 1 yet:

    if (m /= 1) then
       write(*,*) 'm > 1 not allowed.  You have m = ',m
       stop
    end if

! set small default values
    nx1 = 15
    nx2 = 15
    
! initialize error table
    call einfo(0)

! make 2015 recoverable, and print no error messages.
    call errset (2015,0,-1,0,enotrm,2015)
    
! call dcft to get correct value of nx1.  Ignore nx2.
    call drcft(1, x, 0, y, 0, n, m, 1, 1., aux1, nx1, &
         aux2, nx2, *100)
    
100 continue
    
  end subroutine drcft_size


end module fft_work
