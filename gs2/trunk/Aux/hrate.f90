program hrate

  use command_line, only: cl_getarg, cl_iargc

  real, dimension (:), allocatable :: time
  real, dimension(:,:), allocatable :: h_rate, hrateavg
  character (500) :: runname, c
  integer, dimension (1) :: max_loc, jlo, jhi
  integer :: i, imax, ierr

  if (cl_iargc() /= 0) then
     call cl_getarg (1, runname, 500, ierr)
  else
     write (*,*) 'runname = ?'
     read (*,*) runname
  end if

  open (unit=12, file=trim(runname)//'.heat')
  imax=0
  do 
     read (12, *, end=100)     
     imax = imax + 1
  end do

100 rewind(12)
  
  allocate (time(imax), h_rate(imax, 2), hrateavg(0:imax, 2))
  
  hrateavg(0,:) = 0. 
  do i=1,imax
     read (12, *) c, time(i), c, h_rate(i,1), h_rate(i,2)
     hrateavg(i,:) = (hrateavg(i-1,:)*(i-1) + h_rate(i,:))/real(i)
  end do
  close (12) 


  i=imax
  do i=1,imax
     write (*,*) time(i), hrateavg(i,1), hrateavg(i,2), hrateavg(i,1)/hrateavg(i,2)
  end do

end program hrate
