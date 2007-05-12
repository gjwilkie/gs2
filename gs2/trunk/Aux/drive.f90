program drive

  use regression
  implicit none

  integer, parameter :: nr = 1, np = 1000
  real, dimension(nr) :: x, y, a, b, eps
  type (reg_type), dimension(nr) :: r
  integer :: i, j
  real :: avg

  do 
     write (*,*) 'x, y'
     read (*,*,err=100) x(1), y(1)
     call regress (r, x, y)
     write (*,*) r%a, r%b
  end do

100 continue

  write (*,*) 'x, y'
  read (*,*,err=100) x(1), y(1)
  call regress (r, x, y, a, b)

  write (*,*) 'Final values:'
  write (*,*) 'a(1) = ',a(1),'     b(1) = ',b(1)
!  write (*,*) 'a(2) = ',a(2),'     b(2) = ',b(2)

  x(1) = 10.
!  x(2) = 10.
  
  write (*,*) 
  write (*,*) 'x(1) = ',x(1)!,' x(2) = ',x(2)
  write (*,*) 'y= ',yp(r, x)

  y(1) = 36.
!  y(2) = -3.

  write (*,*) 
  write (*,*) 'y(1) = ',y(1)!,' y(2) = ',y(2)
  write (*,*) 'x= ',xp(r, y)


end program drive
