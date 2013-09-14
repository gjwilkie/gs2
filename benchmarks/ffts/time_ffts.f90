
!> A program that repeatedly calls add_nl for benchmarking the ffts and
!! transposes
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program time_ffts
  use unit_tests
  use mp, only: init_mp, finish_mp, proc0, broadcast
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



  call announce_module_test('time_ffts')

  call init_nonlinear_terms

  allocate(g1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  allocate(phi(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar(-ntgrid:ntgrid,ntheta0,naky))

  do i = 1,100
    call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
    if (proc0) write (*,*) 'Finished nonlinear_terms_unit_test_time_add_nl ', i
  end do

  call finish_nonlinear_terms

  call close_module_test('time_ffts')

  call finish_mp

end program time_ffts
