program om

  use command_line, only: cl_getarg, cl_iargc

  real, dimension(:), allocatable :: w, A
  character (500) :: runname, c
  integer, dimension (1) :: max_loc, jlo, jhi
  integer :: i, imax, ierr, len
  complex :: omega
  real :: wr, wi

  if (cl_iargc() /= 0) then
     len = 500
     call cl_getarg (1, runname, len, ierr)
  else
     write (*,*) 'runname = ?'
     read (*,*) runname
  end if

  open (unit=12, file=trim(runname)//'.w')
  imax=0
  do 
     imax = imax + 1
     read (12, *, end=100)     
  end do

100 rewind(12)
  
  allocate (w(imax), A(imax))
  
  do i=1,imax-1
     read (12, *) c, w(i), c, A(i)
  end do
  close (12) 

  max_loc = maxloc (A)

  jlo = minloc ((A(1:max_loc(1))-0.5*A(max_loc(1)))**2)
  jhi = minloc ((A(max_loc(1):imax-1)-0.5*A(max_loc(1)))**2)+max_loc(1)

  wr = w(max_loc(1))
  wi = (w(jhi(1))-w(jlo(1)))/4.

! spot testing (direct evaluation of (wr, wi)) shows need 
! for these corrections.  Should investigate this further.

  wr = wr * 0.98
  wi = wi * 1.11

  omega = cmplx(wr, -wi)
  write (*,*) 'omega/k_par v_A ~ ',omega
  write (*,*) 'Log10 (omega/k_par v_A) ~ ',alog10(real(omega))
  write (*,*) 'Log10 (gamma/omega) = ',alog10(wi/wr)

end program om
