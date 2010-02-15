program ffttester
!
! CMR, 15/2/2010:
! ffttester has simple objective to test GS2 Fourier Transform Routines
!        initialise gs2 with user input file, but with NSTEP=0 to avoid advance
!        initialise g with cos(jx x)*cos(jy y) 
!
  use gs2_main, only: run_gs2, finish_gs2
  use gs2_diagnostics, only: finish_gs2_diagnostics
  use mp, only: finish_mp
  use kt_grids, only: naky, nx, ny, ntheta0, box, theta0

  implicit none
  integer:: ix, iy

  call run_gs2(nofinish=.true.)
  if ( .not. box ) then
     write(6,*) 'Quitting as this is not a BOX run'
     stop
  endif
  do ix=0,ntheta0/2-1
     do iy=0,naky-1
        call ffttest(ix,iy)
     enddo
  enddo
  call finish_gs2_diagnostics(0)
  call finish_gs2
  call finish_mp
end program ffttester




subroutine ffttest(jx,jy,debug)
! CMR, 15/2/2010:  fft testing routine
!                  set initial g=cos(jx x)*cos(jy y) 
!                  and use to test FFTs
!
  use dist_fn_arrays, only: g
  use constants,only: twopi
  use gs2_layouts, only: g_lo, yxf_lo, accelx_lo, accel_lo, idx_local, idx
  use gs2_layouts, only: layout
  use gs2_transforms, only: init_transforms
  use gs2_transforms, only: transform2, inverse2
  use species, only: nspec
  use theta_grid, only: ntgrid
  use kt_grids, only: akx, aky, naky, ikx, nx, ny, ntheta0, box, theta0
  use le_grids, only: negrid, nlambda
  use mp, only: nproc
  implicit none
  integer, intent(in):: jx,jy
  logical, optional, intent(in):: debug
!  ik,is,ie,il,it are indices for ky, species, energy, lambda and kx respectively
  integer:: ik,is,ie,il,it,ia
  integer:: isgn, ig, iglo, index
  logical:: accelerated
  logical,save:: alloc=.true.
  logical:: printlots, fail
  logical, save :: first=.true. 

  real, save, dimension (:,:), allocatable :: gr  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc
  real, save, dimension (:,:,:), allocatable :: gra  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc
  real:: exact, err

  printlots=.false. ; fail=.false.
  if (present(debug)) printlots=debug
  g=0
  if (printlots) write(6,fmt='("nx=",i6,": ny=",i6,": ntheta0=",i6,": naky=",i6)') ,nx,ny,ntheta0,naky
  if (jy .ge. naky) then 
     write(6,*) "ffttest: quit as jy>=naky --- jy, naky=",jy,naky
     return
  endif
  if (jx .gt. (ntheta0-1)/2) then 
     write(6,fmt='("ffttest: quit as jx>ntheta0 --- jx, ntheta0=",2i6)') jx,ntheta0
     return
  endif
  
!  set g independent of v space = cos(jx x) * cos(jy y)
   do is=1,nspec
      do ie=1,negrid
         do il=1,nlambda
            ik=jy+1
            it=jx+1
            iglo = idx(g_lo, ik, it, il, ie, is)
!
! CMR: factor 0.25 here to return FFT = cos(jx x) * cos(jy y)
! 
            if (idx_local (g_lo, iglo)) then
               g(:,:,iglo)=0.25*cmplx(1.0,0.0)
! CMR: if jx=0 and jy/=0 components of g must be doubled
               if (jx .eq.0) g(:,:,iglo)=2.0*g(:,:,iglo)
            endif
            if (jx .gt. 0) then
               it=ntheta0 +1 - jx
               iglo = idx(g_lo, ik, it, il, ie, is)
               if (idx_local (g_lo, iglo)) then
                  g(:,:,iglo)=0.25*cmplx(1.0,0.0)
               endif
            endif
         end do
      end do
   end do

! Now initialise the FFT routines

   call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)

   if (alloc) then
      if (accelerated) then
         if (printlots) then
            write(6,fmt='(2A20)') "accelx_lo","accel_lo"
            write(6,fmt='(2I20,A20)') accelx_lo%iproc,accel_lo%iproc,"iproc"
            write(6,fmt='(2I20,A20)') accelx_lo%ntgrid,accel_lo%ntgrid,"ntgrid"
            write(6,fmt='(2I20,A20)') accelx_lo%nsign,accel_lo%nsign,"nsign"
            write(6,fmt='(2I20,A20)') accelx_lo%naky,accel_lo%naky,"naky"
            write(6,fmt='(2I20,A20)') -999,accel_lo%ndky,"ndky"
            write(6,fmt='(2i20,A20)') accelx_lo%ny,accel_lo%ny,"ny"
            write(6,fmt='(2i20,A20)') accelx_lo%ntheta0,accel_lo%ntheta0,"ntheta0"
            write(6,fmt='(2i20,A20)') accelx_lo%nx,accel_lo%nx,"nx"
            write(6,fmt='(2i20,A20)') accelx_lo%nxny,accel_lo%nxny,"nxny"
            write(6,fmt='(2i20,A20)') -999,accel_lo%nxnky,"nxnky"
            write(6,fmt='(2i20,A20)') accelx_lo%nlambda,accel_lo%nlambda,"nlambda"
            write(6,fmt='(2i20,A20)') accelx_lo%negrid,accel_lo%negrid,"negrid"
            write(6,fmt='(2i20,A20)') accelx_lo%nspec,accel_lo%nspec,"nspec"
            write(6,fmt='(2i20,A20)') -999,accel_lo%nia,"nia"
            write(6,fmt='(2i20,A20)') accelx_lo%llim_world,accel_lo%llim_world,"llim_world"
            write(6,fmt='(2i20,A20)') accelx_lo%ulim_world,accel_lo%ulim_world,"ulim_world"
            write(6,fmt='(2i20,A20)') accelx_lo%blocksize,accel_lo%blocksize,"blocksize"
            write(6,fmt='(2i20,A20)') accelx_lo%llim_proc,accel_lo%llim_proc,"llim_proc"
            write(6,fmt='(2i20,A20)') accelx_lo%ulim_proc,accel_lo%ulim_proc,"ulim_proc"
            write(6,fmt='(2i20,A20)') accelx_lo%ulim_alloc,accel_lo%ulim_alloc,"ulim_alloc"
         endif
         allocate (gra(-ntgrid:ntgrid, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))
         gra=0.
      else
         if (printlots) then
            write(6,fmt='(2A20)') "yxf_lo"
            write(6,fmt='(A20,I20)') "yxf_lo%iproc",yxf_lo%iproc
            write(6,fmt='(A20,I20)') "yxf_lo%ntgrid",yxf_lo%ntgrid
            write(6,fmt='(A20,I20)') "yxf_lo%nsign",yxf_lo%nsign
            write(6,fmt='(A20,I20)') "yxf_lo%naky",yxf_lo%naky
            write(6,fmt='(A20,I20)') "yxf_lo%ny",yxf_lo%ny
            write(6,fmt='(A20,I20)') "yxf_lo%ntheta0",yxf_lo%ntheta0
            write(6,fmt='(A20,I20)') "yxf_lo%nx",yxf_lo%nx
            write(6,fmt='(A20,I20)') "yxf_lo%nlambda",yxf_lo%nlambda
            write(6,fmt='(A20,I20)') "yxf_lo%negrid",yxf_lo%negrid
            write(6,fmt='(A20,I20)') "yxf_lo%nspec",yxf_lo%nspec
            write(6,fmt='(A20,I20)') "yxf_lo%llim_world",yxf_lo%llim_world
            write(6,fmt='(A20,I20)') "yxf_lo%ulim_world",yxf_lo%ulim_world
            write(6,fmt='(A20,I20)') "yxf_lo%llim_proc",yxf_lo%llim_proc
            write(6,fmt='(A20,I20)') "yxf_lo%ulim_proc",yxf_lo%ulim_proc
            write(6,fmt='(A20,I20)') "yxf_lo%ulim_alloc",yxf_lo%ulim_alloc
         endif
         allocate (gr(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc))
         gr=0.
      end if
      alloc = .false.
   end if

! CMR: Perform FFT from k space to real space
!      NB multiply fftw result by 2 .... MUST CHECK
!         presume due to normalisations in fftw complex->real
! 
   if (accelerated) then
      call transform2 (g, gra, ia)
      gra=2.0*gra
   else
      call transform2 (g, gr)
      gr=2.0*gr
   end if

! NB Here we CHECK result at ONLY one point in (velocity,species,theta) space
! at theta=0 along the field line

   is=1 ; ie=1 ; il=1 ; isgn=1 ; ig=0

   if (accelerated) then
      do ik=1,accelx_lo%ny
         do it=1,accelx_lo%nx
            index = idx(accelx_lo, ik, it, il, ie, is)
            if ( idx_local(accelx_lo, index) ) then
               exact=cos(jx*real(it-1)/real(accelx_lo%nx)*twopi) * cos(jy *real(ik-1)/real(accelx_lo%ny)*twopi)
               err=gra(ig,isgn,index)-exact
               if (err .gt. accelx_lo%nx*accelx_lo%ny*epsilon(err)) then
                  write(6,fmt='("GS2FFT, exact, it, ik=",2e12.4,2I8)') gra(ig,isgn,index),exact,it,ik
                  fail=.true. 
               endif
            endif
         enddo
      enddo
   else
      do ik=1,yxf_lo%ny
         do it=1,yxf_lo%nx
            index = idx(yxf_lo, ig, isgn, it, il, ie, is)
            if ( idx_local(yxf_lo, index) ) then
               exact=cos(jx*real(it-1)/real(yxf_lo%nx)*twopi) * cos(jy *real(ik-1)/real(yxf_lo%ny)*twopi)
               err=abs(gr(ik,index)-exact)
               if (err .gt. yxf_lo%nx*yxf_lo%ny*abs(epsilon(err))) then
                  write(6,fmt='("GS2FFT, exact, it, ik=",2e12.4,2I8)') gr(ik,index),exact,it,ik
                  fail=.true. 
               endif
            endif
         enddo
      enddo
   end if

   if (first) write(6,fmt='(7a8,a12)') "jx","jy","nx","ny","nproc","layout","accel","Pass/Fail"
   if (fail) then
      write(6,fmt='(5i8,a8,l8,a12)') jx,jy, nx, ny, nproc, layout, accelerated, "!!FAIL!!"
   else 
      write(6,fmt='(5i8,a8,l8,a12)') jx,jy, nx, ny, nproc, layout, accelerated, "PASS"
   endif
   first=.false.
end

