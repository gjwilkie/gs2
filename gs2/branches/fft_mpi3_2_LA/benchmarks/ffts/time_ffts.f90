
!> A program that repeatedly calls add_nl for benchmarking the ffts and
!! transposes
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program time_ffts
  use mpi
  use shm_mpi3
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_abort, proc0, iproc, broadcast
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

  call shm_alloc(g1, (/-ntgrid, ntgrid, 1, 2, g_lo%llim_proc, g_lo%ulim_proc/))
  call shm_alloc(g, (/-ntgrid, ntgrid, 1, 2, g_lo%llim_proc, g_lo%ulim_proc/))
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

  call shm_free(g1)
  call shm_free(g)
  call close_module_test('time_ffts')

  call finish_mp

  contains

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
