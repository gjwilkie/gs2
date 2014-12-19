
!> A program that tests the nonlinear terms module. It incorporates ffttest,
!! an earlier program  which initialises the distribution
!! function and fields with specific functions and sees if the poisson bracket
!! is correctly calculated.
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
!!               Colin Roach (colin.m.roach@ccfe.ac.uk)

program test_nonlinear_terms
  use unit_tests
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils, run_name
  use species, only: init_species, nspec, spec
  use dist_fn, only: init_dist_fn, finish_dist_fn
  use constants, only: pi
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  use dist_fn_arrays, only: g
  use nonlinear_terms, only: nonlinear_terms_unit_test_time_add_nl
  use kt_grids, only: ntheta0, naky
#ifdef MPI
  use mpi
#endif
  implicit none
  real :: eps
    character (500), target :: cbuff
  real :: tstart
  logical :: dummy=.false.
  integer, dimension(:), allocatable :: sizes
  real, dimension(:,:,:), allocatable :: energy_results
  real :: energy_min
  real :: vcut_local
  integer :: i

  complex, dimension (:,:,:), allocatable :: integrate_species_results
  complex, dimension (:,:,:), allocatable :: g1
  complex, dimension (:,:,:), allocatable :: phi, apar, bpar


  ! General config
  eps = 1.0e-7

  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy, name="gs")
       if (proc0) then
          cbuff = trim(run_name)
       end if
       
       call broadcast (cbuff)
       if (.not. proc0) run_name => cbuff

  


  call announce_module_test('nonlinear_terms')


  call init_dist_fn

  allocate(g1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  !allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  allocate(phi(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar(-ntgrid:ntgrid,ntheta0,naky))

  call announce_test('run nonlinear terms')
  call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
  call process_test(.true., 'run nonlinear terms')

  call announce_test('fast fourier transforms in nonlinear term')
  call process_test(test_ffts(),'fast fourier transforms in nonlinear term')


  call finish_dist_fn

  call close_module_test('nonlinear_terms')

  call finish_mp


contains

function test_ffts()
  use unit_tests
  use kt_grids, only: ntheta0, naky
  implicit none
  logical :: test_ffts
  logical ::dummy
  integer :: ix, iy
  character(len=22) :: message = ''
  dummy = .true.
  

  !ix =0
  !iy=0
    !write(message, fmt="(A12, I2, A6, I2)") 'ffttest ix =', ix, ' iy = ', iy
    !write(message, fmt="(A12, I2, A6, I2)") 'ffttest ix =', ix, ' iy = ', iy
  do ix=0,ntheta0/2-1
    do iy=0,naky-1
       !write(message, '(22A1)') ' '

       !write(*, fmt="(A12, I2, A6, I2)") 'ffttest ix =', ix, ' iy = ', iy
       !write(*, *) message
       write(message, fmt="(A12, I2, A6, I2)") 'ffttest ix =', ix, ' iy = ', iy
       !write(message, *) 'ffttest ix =', ix, ' iy = ', iy
       !write(*, *) message
      !write(message, '("ffttest ix = ")') 
      !message = 'Hello'
      call announce_check(message)
      call process_check(test_ffts, ffttest(ix,iy, .false.), message)
    enddo
  enddo

  test_ffts = dummy

end function test_ffts

function ffttest (jx,jy,debug)
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
  use le_grids, only: integrate_moment
  use mp, only: iproc, nproc, proc0, barrier
#ifdef MPI
  use mpi
#endif
  implicit none
  integer, intent(in):: jx,jy
  logical, optional, intent(in):: debug
  logical :: ffttest
!  ik,is,ie,il,it are indices for ky, species, energy, lambda and kx respectively
  integer:: ik,is,ie,il,it,ia
  integer:: isgn, ig, iglo, index, mpierr
  logical:: printlots, fail, anyfail
  logical, save :: accelerated
  logical,save:: first=.true. 

  real, save, dimension (:,:), allocatable :: gr  ! yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc
  real, save, dimension (:,:,:), allocatable :: gra  ! 2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc
  real:: exact, err
  complex:: fexact

  complex, dimension (:,:,:,:), allocatable :: density,cdens
  complex, dimension (:,:,:), allocatable :: d3d,c3d
  real, dimension (:,:,:), allocatable :: rdens  
  integer:: nnx,nny


  ffttest = .true.


  printlots=.false. 
  if (present(debug)) printlots=debug

!CMR, 5-D FFTs
  g=0
  if (printlots) write(6,fmt='("nx=",i6,": ny=",i6,": ntheta0=",i6,": naky=",i6)') nx,ny,ntheta0,naky
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
            if (idx_local(g_lo, iglo)) then
!if (printlots) write(6,fmt='("Processor,iglo,lower/upper bounds=",4i8)') iproc,iglo,lbound(g,3),ubound(g,3)
!if (printlots) write(6,fmt='("Processor, iglo, ik, it, il, ie, is= ",7i5)') iproc,iglo,ik,it,il,ie,is
!if (printlots) flush(6)
               g(:,:,iglo)=0.25*cmplx(1.0,0.0)
!if (printlots) write(6,fmt='("Processor ",i4," successful write")') iproc
!if (printlots) flush(6)

! CMR: if jx=0 and jy/=0 components of g must be doubled
               if (jx .eq.0) g(:,:,iglo)=2.0*g(:,:,iglo)
            endif
            if (jx .gt. 0) then
               it=ntheta0 +1 - jx
               iglo = idx(g_lo, ik, it, il, ie, is)
               if (idx_local(g_lo, iglo)) then
                  g(:,:,iglo)=0.25*cmplx(1.0,0.0)
               endif
            endif
if (printlots) call barrier
         end do
      end do
   end do
   if (first) then
! initialise FFT routines
      call init_transforms (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
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
     !write (*,*) accelerated, 'accelerated', iproc
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
   end if

   if (first .and. proc0 .and. printlots) write(6,fmt='(9a8,a12)') "Dim","Dir","jx","jy","nx","ny","nproc","layout","accel","Pass/Fail"

! CMR, compute 3-D FFT
   fail = .false.
   allocate (d3d(-ntgrid:ntgrid,ntheta0,naky),c3d(-ntgrid:ntgrid,ntheta0,naky))
   d3d=0.0
   d3d(:,jx+1,jy+1)=0.25*cmplx(1.0,0.0)
   if (jx .gt. 0) d3d(:,ntheta0+1-jx,jy+1)=0.25*cmplx(1.0,0.0)
   if (jx .eq. 0) d3d=2.0*d3d

   nnx=2*nx; nny=2*ny
   allocate (rdens(nnx,nny,-ntgrid:ntgrid))
   call transform2(d3d,rdens,nny,nnx)
   call inverse2(rdens,c3d,nny,nnx)
   rdens=2.0*rdens

! CHECK results at theta=0 
   ig=0
   do ik=1,nny
      do it=1,nnx
         exact=cos(jx*real(it-1)/real(nnx)*twopi) * cos(jy *real(ik-1)/real(nny)*twopi)
         err=abs(rdens(it,ik,ig)-exact)
         if (err .gt. nnx*nny*abs(epsilon(err))) then
            if (printlots) write(6,fmt='("ffttest: transform2_3d, exact, err, it, ik=",3e16.8,2I8)') rdens(it,ik,ig),exact,err,it,ik
            fail=.true. 
         endif
      enddo
   enddo
#ifdef MPI
   call mpi_reduce(fail,anyfail,1,MPI_LOGICAL,MPI_LOR,0,mpi_comm_world,mpierr)
   if (mpierr .ne. MPI_SUCCESS) write(6,*) "ffttest: mpi_reduce gives mpierr=",mpierr
#else
   anyfail = fail
#endif
   if (proc0 .and. printlots) then
      if (anyfail) then
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "3-D","c2r",jx,jy, nx, ny, nproc, layout, .false., "!!FAIL!!"
      else 
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "3-D","c2r",jx,jy, nx, ny, nproc, layout, .false., "PASS"
      endif
   endif
   if (anyfail) ffttest = .false.
! check inverse
   fail=.false.
   do ik=1,naky
      do it=1,ntheta0
         err=abs(d3d(ig,it,ik)-c3d(ig,it,ik))
         if (err .gt. nnx*nny*abs(epsilon(err))) then
            if (printlots) write(6,fmt='("ffttest: inverse2_3d, exact, err, it, ik=",5e16.8,2I8)') c3d(ig,it,ik),d3d(ig,it,ik),err,it,ik
            fail=.true. 
         endif
      enddo
   enddo
#ifdef MPI
   call mpi_reduce(fail,anyfail,1,MPI_LOGICAL,MPI_LOR,0,mpi_comm_world,mpierr)
   if (mpierr .ne. MPI_SUCCESS) write(6,*) "ffttest: mpi_reduce gives mpierr=",mpierr
#else
   anyfail = fail
#endif
   if (proc0 .and. printlots) then
      if (anyfail) then
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "3-D","r2c",jx,jy, nx, ny, nproc, layout, .false., "!!FAIL!!"
      else 
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "3-D","r2c",jx,jy, nx, ny, nproc, layout, .false., "PASS"
      endif
   endif
   if (anyfail) ffttest = .false.

   deallocate(d3d,rdens,c3d)


! CMR, compute 4-D FFT
   fail = .false.
   allocate (density(-ntgrid:ntgrid,ntheta0,naky,nspec),cdens(-ntgrid:ntgrid,ntheta0,naky,nspec))
   density=0.0
   density(:,jx+1,jy+1,:)=0.25*cmplx(1.0,0.0)
   if (jx .gt. 0) density(:,ntheta0+1-jx,jy+1,:)=0.25*cmplx(1.0,0.0)
   if (jx .eq. 0) density=2.0*density

   nnx=2*nx; nny=2*ny
   allocate (rdens(nnx,nny,-ntgrid:ntgrid))
   call transform2(density,rdens,nny,nnx)
   rdens=2.0*rdens

! CHECK results at theta=0 for 1st species
   is=1 ; ig=0
   do ik=1,nny
      do it=1,nnx
         exact=cos(jx*real(it-1)/real(nnx)*twopi) * cos(jy *real(ik-1)/real(nny)*twopi)
         err=abs(rdens(it,ik,ig)-exact)
         if (err .gt. nnx*nny*abs(epsilon(err))) then
            if (printlots) write(6,fmt='("ffttest: transform2_4d, exact, err, it, ik=",3e16.8,2I8)') rdens(it,ik,ig),exact,err,it,ik
            fail=.true. 
         endif
      enddo
   enddo
#ifdef MPI
   call mpi_reduce(fail,anyfail,1,MPI_LOGICAL,MPI_LOR,0,mpi_comm_world,mpierr)
   if (mpierr .ne. MPI_SUCCESS) write(6,*) "ffttest: mpi_reduce gives mpierr=",mpierr
#else
   anyfail = fail
#endif
   if (proc0 .and. printlots) then
      if (anyfail) then
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "4-D","c2r",jx,jy, nx, ny, nproc, layout, .false., "!!FAIL!!"
      else 
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "4-D","c2r",jx,jy, nx, ny, nproc, layout, .false., "PASS"
      endif
   endif
   if (anyfail) ffttest = .false.
   deallocate(density,rdens,cdens)

! CMR, Now for 5-D FFTs
! CMR: Perform FFT from k space to real space
!      NB multiply fftw result by 2 .... MUST CHECK
!         presume due to normalisations in fftw complex->real
! 

   fail=.false.
   if (accelerated) then
     !write (*,*) accelerated, 'accelerated', iproc
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
                  if (printlots) write(6,fmt='("ffttest: GS2FFT, exact, it, ik=",2e12.4,2I8)') gra(ig,isgn,index),exact,it,ik
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
                  if (printlots) write(6,fmt='("ffttest: GS2FFT, exact, it, ik=",2e12.4,2I8)') gr(ik,index),exact,it,ik
                  fail=.true. 
               endif
            endif
         enddo
      enddo
   end if
#ifdef MPI
   call mpi_reduce(fail,anyfail,1,MPI_LOGICAL,MPI_LOR,0,mpi_comm_world,mpierr)
   if (mpierr .ne. MPI_SUCCESS) write(6,*) "ffttest: mpi_reduce gives mpierr=",mpierr
#else
   anyfail = fail
#endif

   if (proc0 .and. printlots) then
      if (anyfail) then
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "5-D","c2r",jx,jy, nx, ny, nproc, layout, accelerated, "!!FAIL!!"
      else 
         write(6,fmt='(2a8,5i8,a8,l8,a12)') "5-D","c2r",jx,jy, nx, ny, nproc, layout, accelerated, "PASS"
      endif
   endif
   if (anyfail) ffttest = .false.

! NOW check the INVERSE transform

   if (accelerated) then
      gra=0.5*gra
      call inverse2 (gra, g, ia)
   else
      gr=0.5*gr
      call inverse2 (gr, g)
   end if

! CHECK result at ONLY one point in (velocity,species,theta) space
! at theta=0 along the field line

   is=1 ; ie=1 ; il=1 ; isgn=1 ; ig=0 ; fail=.false.
   do ik=1,g_lo%naky
      do it=1,g_lo%ntheta0
         fexact=cmplx(0.0,0.0)
         if (it .eq. jx+1 .and. ik.eq.jy+1) fexact=0.25*cmplx(1.0,0.0)
         if (it .eq. jx+1 .and. ik.eq.jy+1 .and. jx.eq.0) fexact=0.5*cmplx(1.0,0.0)
         if (jx .gt. 0  .and. ik.eq.jy+1 .and. it .eq. g_lo%ntheta0+1-jx) fexact=0.25*cmplx(1.0,0.0)
         iglo = idx(g_lo, ik, it, il, ie, is)
         if (idx_local(g_lo, iglo)) then
            err=abs(g(ig,isgn,iglo)-fexact)
            if (err .gt. epsilon(err)) then
               if (printlots) write(6,fmt='("ffttest: it, ik, g, fexact=",2I8,4e12.4)') it,ik,g(ig,isgn,iglo),fexact
               fail=.true. 
            endif
         endif
      enddo
   enddo
#ifdef MPI
   call mpi_reduce(fail,anyfail,1,MPI_LOGICAL,MPI_LOR,0,mpi_comm_world,mpierr)
   if (mpierr .ne. MPI_SUCCESS) write(6,*) "ffttest: mpi_reduce gives mpierr=",mpierr
#else
   anyfail = fail
#endif

    if (anyfail) then
      if(proc0 .and. printlots) write(6,fmt='(2a8,5i8,a8,l8,a12)') "5-D","r2c",jx,jy, nx, ny, nproc, layout, accelerated, "!!FAIL!!"
    else 
      if(proc0 .and. printlots) write(6,fmt='(2a8,5i8,a8,l8,a12)') "5-D","r2c",jx,jy, nx, ny, nproc, layout, accelerated, "PASS"
    endif

   if (anyfail) ffttest = .false.

   ! anyfail is only set on proc0... hence we broadcast ffttest
   call broadcast(ffttest)

   first=.false.
end function ffttest

end program test_nonlinear_terms
