
!> A program that repeatedly calls add_nl for benchmarking the ffts and
!! transposes
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program time_ffts
  use mpi
  use FIPC_module
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_abort, proc0, iproc, broadcast, shm_info
  use file_utils, only: init_file_utils, run_name
  use species, only: init_species, nspec, spec
  use constants, only: pi
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  use dist_fn_arrays, only: g
  use nonlinear_terms, only: nonlinear_terms_unit_test_time_add_nl
  use kt_grids, only: ntheta0, naky
  use job_manage, only: time_message
  implicit none
  real :: eps
    character (500), target :: cbuff
  real :: tstart
  logical :: dummy=.false.
  integer, dimension(:), allocatable :: sizes
  real, dimension(:,:,:), allocatable :: energy_results
  real :: energy_min
  real :: vcut_local
  real :: time_taken(2) = 0.0
  real (kind(0.d0)) ts, te
  integer :: i, is, isn, ie, ien, nruns, ierr

  complex, dimension (:,:,:), allocatable :: integrate_species_results
  complex, dimension (:,:,:), pointer :: g1 => null()
  complex, dimension (:,:,:), allocatable :: phi, apar, bpar

  real ng, ng1 ! g g1 norms**2

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
  
  call announce_module_test('time_ffts')
  
  call init_nonlinear_terms

  isn = shm_info%g_lo_se(1, 0)
  ien = shm_info%g_lo_se(2, shm_info%size - 1)
  is = shm_info%g_lo_se(1, shm_info%id) 
  ie = shm_info%g_lo_se(2, shm_info%id)
  allocate(shm_info%g_lo_ptr(2))
  Call FIPC_seg_create( FIPC_ctxt_world, (/ (2*ntgrid+1), 2, (ien - isn + 1)  /),shm_info%g_lo_ptr(1)%p , ierr)
  if (ierr /= FIPC_success ) then 
     write(0,*)' FIPC error 1', ierr
     call mp_abort("time_ffts: error in shared segment allocation")
  endif

  shm_info%g_lo_ptr(1)%p => remap_bounds(-ntgrid, 1, isn, shm_info%g_lo_ptr(1)%p)
 
  g1 => shm_info%g_lo_ptr(1)%p(:,:, is:ie)
  g1 => remap_bounds(-ntgrid, 1, g_lo%llim_proc, g1)
  !write(0,*) 'remap g1', iproc,lbound(g1),ubound(g1)

  Call FIPC_seg_create( FIPC_ctxt_world, (/ (2*ntgrid+1), 2, ien - isn + 1 /),shm_info%g_lo_ptr(2)%p, ierr)
  if (ierr /= FIPC_success ) then
     write(0,*)' FIPC error 2', ierr
     call mp_abort("time_ffts: error in shared segment allocation")
  endif
  shm_info%g_lo_ptr(2)%p => remap_bounds(-ntgrid, 1, isn, shm_info%g_lo_ptr(2)%p)
  g =>  shm_info%g_lo_ptr(2)%p(:,:, is:ie)
  g => remap_bounds(-ntgrid, 1, g_lo%llim_proc, g)
  !write(0,*) 'remap g', iproc, lbound(g),ubound(g)


  !allocate(g1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  !allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  allocate(phi(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar(-ntgrid:ntgrid,ntheta0,naky))

  !if (proc0) call time_message(.false., time_taken, "FFT time")
  ts = mpi_wtime()

  nruns = 5

  do i = 1, nruns
    call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
    !if (proc0) write (*,*) 'Finished nonlinear_terms_unit_test_time_add_nl ', i
  end do

  te = mpi_wtime()
  call time_statistics(ts,te, nruns)
  !if (proc0) then
  !  call time_message(.false., time_taken, "FFT time")
  !  write(*, '(" Time for nonlinear_terms: ",F6.1," s")') time_taken(1)
  !  write(*,*)
  !end if

  ! compute the norms of g and g1
  call compute_global_norm(g,ng)
  call compute_global_norm(g1,ng1)
  
  if (proc0) then
     write(*,*) "Norms g", ng, "g1", ng1
  endif

  call finish_nonlinear_terms

  call FIPC_seg_free(shm_info%g_lo_ptr(1)%p, ierr)
  call FIPC_seg_free(shm_info%g_lo_ptr(2)%p, ierr)
  call close_module_test('time_ffts')

  call finish_mp

  contains

    FUNCTION remap_bounds(lb1,lb2,lb3,array) RESULT(ptr)
      INTEGER, INTENT(IN)                            :: lb1,lb2,lb3
      complex, DIMENSION(lb1:,lb2:,lb3:), INTENT(IN), TARGET :: array
      complex, DIMENSION(:,:,:), POINTER                  :: ptr
      ptr => array
    END FUNCTION remap_bounds

    subroutine compute_global_norm(g, norm)
      use mp, only : sum_reduce
      implicit none
      complex, intent(in) :: g(:,:,:)
      real, intent(out) :: norm
      
      real ln
      
      ln = sum(abs(g)*abs(g))
      
      call sum_reduce(ln,0)
      
      norm = ln
      
    end subroutine compute_global_norm


    subroutine time_statistics(ts, te, nruns)
      use mp, only : min_reduce, max_reduce, sum_reduce, proc0, nproc
      implicit none

      real(kind(0.d0)), intent(in) :: ts, te
      integer, intent(in) :: nruns
      integer ierr
      real t, tsum, tmin, tmax

      t = te - ts
      tmin = t
      call min_reduce(tmin, 0)
      tmax = t
      call max_reduce(tmax, 0)
      tsum = t
      call sum_reduce(tsum, 0)

      if (proc0) then 
         write(*,'(a,/,i6,1x,3(E10.3,1x))')"# FFT times per run: nproc, average, min, max ",&
              nproc, tsum/(nruns*real(nproc)), tmin/nruns, tmax/nruns
      endif

    end subroutine time_statistics


end program time_ffts
