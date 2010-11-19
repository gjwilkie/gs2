program kp

  use command_line, only: cl_getarg, cl_iargc

  real, dimension(:,:,:), allocatable :: vx2, vy2
  real, dimension (:), allocatable :: kpar, kx, ky, vp2
  character (500) :: runname, c
  integer, dimension (1) :: max_loc, jlo, jhi
  integer :: i, imax, ierr
  complex :: omega
  real :: wr, wi
  integer :: ntheta0, naky, nkpar

  if (cl_iargc() /= 0) then
     call cl_getarg (1, runname, 500, ierr)
  else
     write (*,*) 'runname = ?'
     read (*,*) runname
  end if

  open (unit=12, file=trim(runname)//'.gs')
!  imax=0
!  do 
!     imax = imax + 1
!     read (12, *, end=100)     
!  end do

!100 rewind(12)
  
  write (*,*) 'how many kx values?'
  read (*,*) ntheta0
  
  write (*,*) 'how many ky values?'
  read (*,*) naky
  
  write (*,*) 'how many kpar values?'
  read (*,*) nkpar
  

  allocate (kpar(nkpar), vp2(nkpar))
  allocate (vx2(nkpar, ntheta0, naky))
  allocate (vy2(nkpar, ntheta0, naky))
  allocate (kx(ntheta0))
  allocate (ky(naky))
  
  do ik = 1, naky
     do it=1, ntheta0
        do ig=1,nkpar
           read (12, *) kpar(ig), ky(ik), kx(it), vx2(ig, it, ik), vy2(ig, it, ik)
        end do
        read (12, *)
     end do
  end do
  close (12) 

  do ig=1,nkpar
     vp2(ig) = 0.
     do ik=1, naky
        fac = 0.5
        if (ky(ik) < 2.*epsilon(0.0)) fac = 1.0
        do it = 1, ntheta0
           vp2(ig) = vp2(ig) + (vx2(ig,it,ik) + vy2(ig,it,ik))*fac
        end do
     end do
  end do

  do ig=nkpar/2+2, nkpar-1
     vp2(ig) = vp2(ig)+vp2(nkpar-ig+1)
  end do

  do ig=nkpar/2+1, nkpar-1
     write (*,*) kpar(ig), vp2(ig)
  end do

end program kp
