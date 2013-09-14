
!> A program that repeatedly calls add_nl
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program time_ffts
  use unit_tests
  use le_grids, only: init_le_grids, finish_le_grids
  !use egrid
  !use general_f0, only: init_general_f0
  use mp, only: init_mp, finish_mp, proc0, broadcast
  use file_utils, only: init_file_utils, run_name
  use species, only: init_species, nspec, spec
  use constants, only: pi
  !use fields, only: init_fields
  !use fields_arrays, only: phi
  !use dist_fn, only: init_dist_fn
  use dist_fn_arrays, only: g
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  !use gs2_transforms, only: init_transforms
  use nonlinear_terms, only: init_nonlinear_terms, finish_nonlinear_terms
  use dist_fn_arrays, only: g
  use nonlinear_terms, only: nonlinear_terms_unit_test_time_add_nl
  use kt_grids, only: ntheta0, naky
  !use gs2_main, only: finish_gs2
  !use parameter_scan, only: init_parameter_scan, allocate_target_arrays
  !use gs2_diagnostics, only: init_gs2_diagnostics
  !use gs2_time, only: init_tstart
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

  !complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: app, bpar
  !complex, dimension (:,:,:), allocatable :: g
  
  !real, dimension(:), allocatable :: species_weights
  !real :: ene
  !integer :: iglo
  !integer :: ie



  ! General config
  eps = 1.0e-7

  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy, name="gs")
  !call init_fields
  !call init_species
  !call init_general_f0
       if (proc0) then
          !call time_message(.false., time_init,' Initialization')
          cbuff = trim(run_name)
       end if
       
       call broadcast (cbuff)
       if (.not. proc0) run_name => cbuff



  call announce_module_test('time_ffts')

  !call init_dist_fn
  !call init_theta_grid
  if (proc0) write (*,*) 'init_theta_grid'
  !call init_kt_grids
  if (proc0) write (*,*) 'kt_grids'
  !call init_gs2_layouts
  if (proc0) write (*,*) 'init_gs2_layouts'
  call init_nonlinear_terms
  if (proc0) write (*,*) 'init_nonlinear_terms'
  !call init_le_grids(dummy, dummy)
       !call init_parameter_scan
       !call init_fields
       !call init_gs2_diagnostics (dummy, 10)
       !call allocate_target_arrays ! must be after init_gs2_diagnostics
       !call init_tstart (tstart)   ! tstart is in user units 
  !write (*,*) 

  allocate(g1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  allocate(phi(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar(-ntgrid:ntgrid,ntheta0,naky))
  if (proc0) write (*,*) 'allocated'
  do i = 1,100
    call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
    if (proc0) write (*,*) 'Finished nonlinear_terms_unit_test_time_add_nl ', i
  end do
  !call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
  

  !stop

  !call finish_le_grids
  call finish_nonlinear_terms

  call close_module_test('time_ffts')
  !call finish_gs2

  call finish_mp
  write (*,*)  'called finish_mp'

end program time_ffts
