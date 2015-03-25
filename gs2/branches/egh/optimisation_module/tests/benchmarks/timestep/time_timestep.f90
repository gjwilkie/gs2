
!> A program that repeatedly calls add_nl for benchmarking the ffts and
!! transposes
!!
!! This is free software released under the MIT licence
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program time_timestep
  use unit_tests
  use mp, only: init_mp, finish_mp, proc0, broadcast, nproc
  use file_utils, only: init_file_utils, run_name
  use file_utils, only: append_output_file, close_output_file
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
  !use runtime_tests, only: get_svn_rev, get_compiler_name
  use benchmarks, only: benchmark_identifier
  use fields, only: init_fields, advance
  use fields_implicit, only: skip_initialisation
  use fft_work, only: save_wisdom, load_wisdom
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
  real :: time_init(2) = 0.0
  integer :: i
  integer :: timing_unit

  complex, dimension (:,:,:), allocatable :: integrate_species_results
  complex, dimension (:,:,:), allocatable :: g1
  complex, dimension (:,:,:), allocatable :: phi, apar, bpar


  skip_initialisation = .true.


  !call load_wisdom("wisdom_file")

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



  call announce_module_test('time_timestep')  

  if (proc0) call time_message(.false., time_init, "dummy initialisation time")
  call init_fields
  if (proc0) then
    call time_message(.false., time_init, "dummy initialisation time")
    write(*, '(" Time for dummy initialisation on ",I6," procs: ",F5.1," s")') nproc, time_init(1)
  end if

  if (proc0) call time_message(.false., time_taken, "advance time")
  do i = 1,10
    call advance(i)
  end do
  if (proc0) then
    call time_message(.false., time_taken, "advance time")
    write(*, '(" Time for 10 advance steps on ",I6," procs: ",F5.1," s")') nproc, time_taken(1)
    call append_output_file(timing_unit, &
      benchmark_identifier())
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
    call close_output_file(timing_unit)
  end if

  !call init_nonlinear_terms

  !allocate(g1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  !allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  !allocate(phi(-ntgrid:ntgrid,ntheta0,naky))
  !allocate(apar(-ntgrid:ntgrid,ntheta0,naky))
  !allocate(bpar(-ntgrid:ntgrid,ntheta0,naky))

  !if (proc0) call time_message(.false., time_taken, "FFT time")

  !do i = 1,50
    !call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
    !!if (proc0) write (*,*) 'Finished nonlinear_terms_unit_test_time_add_nl ', i
  !end do

  !if (proc0) then
    !call time_message(.false., time_taken, "FFT time")
    !write(*, '(" Time for nonlinear_terms on ",I6," procs: ",F3.1," s")') nproc, time_taken(1)
    !write(*,*)
    !call append_output_file(timing_unit, &
      !benchmark_identifier())
    !write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
    !call close_output_file(timing_unit)
  !end if

  !call finish_nonlinear_terms

  call close_module_test('time_timestep')

  !call save_wisdom("wisdom_file")

  call finish_mp

end program time_timestep
