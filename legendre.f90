module legendre

! Taken from Numerical Recipes, incorporated here by Tomo Tatsuno
! Aug, 2005

  implicit none

  public :: nrgauleg

  private

contains

  subroutine nrgauleg (x1, x2, x, w, eps)

    real, intent(in) :: x1, x2
    real, dimension(:), intent(out) :: x, w
    real, intent(in) :: eps

    integer :: its, j, m, n
    integer, parameter :: maxit=100
    double precision :: xl, xm, pi
    double precision, dimension((size(x)+1)/2) :: p1, p2, p3, pp, z, z1
    logical, dimension((size(x)+1)/2) :: unfinished

    n = size(x)
    pi = asin(real(1.0,kind(pi)))*2.0
    m = (n+1)/2

    xm = real(0.5,kind(xm)) * (x1+x2)   ! middle of the section
    xl = real(0.5,kind(xl)) * (x2-x1)   ! signed half length of the section
    z = (/ (cos(pi*(j-0.25)/(n+0.5)), j=1,m) /)
    unfinished = .true.

    do its=1, maxit
       where(unfinished)
          p1 = real(1.0,kind(p1(1)))
          p2 = real(0.0,kind(p2(1)))
       end where
       do j=1, n
          where (unfinished)
             p3 = p2
             p2 = p1
             p1 = ((2*j-1) * z * p2 - (j-1) * p3) / j
          end where
       end do
! p1 now contains the desired legendre polynomials.
       where (unfinished)
          pp = n * (z * p1 - p2) / (z**2 - 1.0)
          z1 = z
          z = z1 - p1 / pp
          unfinished = (abs(z-z1) > eps)
       end where
       if (.not. any(unfinished)) exit
    end do

    if (its == maxit+1) then
       print*, 'too many iterations in nrgauleg'
       stop
    end if
    x(1:m) = xm - xl * z
    x(n:n-m+1:-1) = xm + xl * z
    w(1:m) = 2.0 * abs(xl) / ((1.0 - z**2) * pp**2)
    w(n:n-m+1:-1) = w(1:m)

  end subroutine nrgauleg
  
end module legendre
