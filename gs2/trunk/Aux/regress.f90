program regressp

  use regression
  type (reg_type) :: var, varvec

  integer, parameter :: npoints = 100

  real, dimension(npoints) :: epsilon
  real :: x, y, a, b
  real, dimension (npoints) :: xv, yv


  call random_number (epsilon)
  
  do i=1,npoints
     x = real(i)
     y = 3.*real(i)+epsilon(i)
     call regress (var, x, y)
  end do

  y = yp (var, 10.)

  write (*,*) y


  write (*,*) 
  write (*,*) 'And now with vectors...'

! now do it again with vectors

  call random_number (epsilon)
  
  do i=1,npoints
     xv(i) = real(i)
  end do
     
  yv = 3.*xv
  call regress (varvec, xv, yv)

  y = yp (var, 10.)

  write (*,*) y

  

end program regressp
