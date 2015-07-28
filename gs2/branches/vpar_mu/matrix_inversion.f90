module matrix_inversion

  implicit none

  public :: invert_matrix

  interface invert_matrix
     module procedure invert_dmatrix
     module procedure invert_zmatrix
  end interface invert_matrix
  
contains

  subroutine invert_dmatrix (A)

    implicit none

    double precision, dimension (:,:), intent (in out) :: A

    integer :: n, lwork, ierr
    double precision :: norm, rcond, dlange

    integer, dimension (:), allocatable :: ipiv
    integer, dimension (:), allocatable :: iwork
    double precision, dimension (:), allocatable :: work

    n = size(A,1)

    lwork = 4*n

    allocate (ipiv(n))
    allocate (work(lwork))
    allocate (iwork(n))

    ! get the LU factorization of A
    call dgetrf (n, n, A, n, ipiv, ierr)
    if (ierr /=0) then
       write (*,*) 'Error in dgetrf: ', ierr
    end if
    
    norm = dlange ('1', n, n, A, n, work)
!    write (*,*) 'norm: ', norm
    
    ! get reciprocal of condition number of matrix A
    ! from LU factorization obtained by dgetrf
    call dgecon ('1', n, A, n, norm, rcond, work, iwork, ierr)
    if (ierr /=0 ) then
       write (*,*) 'Error in dgecon: ', ierr
    end if
    
!    if (abs(rcond) > epsilon(0.0)) then
!       write (*,*) 'condition number: ', 1./rcond
!       write (*,*)
!    end if
    
    ! get inverse of A using LU factorization
    ! obtained by dgetrf
    call dgetri (n, A, n, ipiv, work, lwork, ierr)
    if (ierr /=0 ) then
       write (*,*) 'Error in dgetri: ', ierr
    end if

    deallocate (ipiv, work, iwork)
    
  end subroutine invert_dmatrix

  subroutine invert_zmatrix (A)

    implicit none

    complex (kind=8), dimension (:,:), intent (in out) :: A

    integer :: n, lwork, ierr
    double precision :: norm, rcond, zlange

    integer, dimension (:), allocatable :: ipiv
    double precision, dimension (:), allocatable :: rwork
    double precision, dimension (:), allocatable :: work
    complex (kind=8), dimension (:), allocatable :: cwork

    n = size(A,1)

    lwork = 2*n

    allocate (ipiv(n))
    allocate (work(1)) ! dummy array not actually accessed
    allocate (cwork(lwork))
    allocate (rwork(lwork))

    ! get the LU factorization of A
    call zgetrf (n, n, A, n, ipiv, ierr)
    if (ierr /=0) then
       write (*,*) 'Error in zgetrf: ', ierr
    end if
    
    norm = zlange ('1', n, n, A, n, work)
!    write (*,*) 'norm: ', norm
    
    ! get reciprocal of condition number of matrix A
    ! from LU factorization obtained by dgetrf
    call zgecon ('1', n, A, n, norm, rcond, cwork, rwork, ierr)
    if (ierr /=0 ) then
       write (*,*) 'Error in zgecon: ', ierr
    end if
    
!    if (abs(rcond) > epsilon(0.0)) then
!       write (*,*) 'condition number: ', 1./rcond
!       write (*,*)
!    end if
    
    ! get inverse of A using LU factorization
    ! obtained by dgetrf
    call zgetri (n, A, n, ipiv, cwork, lwork, ierr)
    if (ierr /=0 ) then
       write (*,*) 'Error in zgetri: ', ierr
    end if

    deallocate (ipiv, work, cwork, rwork)
    
  end subroutine invert_zmatrix

end module matrix_inversion
