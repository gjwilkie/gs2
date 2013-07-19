
!> A program that runs unit tests on the le_grids module.
!! The test results were calculated using sage and are viewable at
!! http://www.sagenb.org/home/pub/5036
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_le_grids
  use unit_tests
  use le_grids
  use egrid
  use general_f0, only: init_general_f0
  use mp, only: init_mp, finish_mp, proc0
  use file_utils, only: init_file_utils
  use species, only: init_species, nspec, spec
  use constants, only: pi
  !use fields, only: init_fields
  !use fields_arrays, only: phi
  !use dist_fn, only: init_dist_fn
  !use dist_fn_arrays, only: g
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use theta_grid, only: ntgrid, init_theta_grid
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  implicit none
  real :: eps
  logical :: dummy
  integer, dimension(:), allocatable :: sizes
  real, dimension(:,:,:), allocatable :: energy_results
  real :: energy_min
  real :: vcut_local

  complex, dimension (:,:,:), allocatable :: integrate_species_results
  complex, dimension (:,:,:), allocatable :: g
  
  real, dimension(:), allocatable :: species_weights
  real :: ene
  integer :: iglo
  integer :: ie



  ! General config
  eps = 1.0e-7

  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy)
  call init_species
  call init_general_f0



  call announce_module_test('le_grids (adjust source)')

  !call announce_test('init_le_grids')
!
  !allocate(sizes(1))
  allocate(energy_results(8,3,2)) ! negrid, nspec, nresults 
  !sizes(1) = 24 ! Size of energy grid

  ! Resuls for energy grids
  energy_results(:,1,1) = (/4.0468821658005425E-003,  0.10438457502758539, 0.55159372522148142, 1.5625000000000000, 3.0881259213302132, 4.7389544850238012, 5.9359713343080411, 7.2500000000000000/) ! Quadrature grid for Maxwellian... this was not calculated inpendently but was taken from the code output
  energy_results(:,2,1) = energy_results(:,1,1)
  energy_results(:,3,1) = (/1.5104763817839494E-002, 4.6790434219629591E-002, 0.13496048316473824, 0.30250000000000005, 0.53674718302836155, 0.78090630796303029, 0.95472159703717052, 2.0000000000000000/) ! Quadrature grid for alphas... this was not calculated inpendently but was taken from the code output

  call init_le_grids(dummy, dummy)

  !deallocate(sizes)


  !call init_dist_fn
  call init_theta_grid
  call init_kt_grids
  call init_gs2_layouts
  !call init_fields

  allocate(sizes(3))
  sizes(1) = 1 !naky
  sizes(2) = 1 !ntheta0
  sizes(3) = 4 !ntgrid
  
  allocate(integrate_species_results(-ntgrid:ntgrid,naky,ntheta0))
  allocate(species_weights(nspec))
  allocate(g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))

  species_weights = 0.0
  g = cmplx(1.0, 1.0)
  !species_weights = 0.0
  species_weights(3) = 1.0
  integrate_species_results = cmplx(1.0,1.0)
  call announce_test('integrate species alphas only g = 1')
  call process_test(&
    le_grids_unit_test_integrate_species(&
      g, &
      species_weights, &
      sizes, &
      integrate_species_results, &
      0.01),&
    'integrate species alphas only  g = 1')


  deallocate(energy_results)
  deallocate(integrate_species_results)
  deallocate(species_weights)
  call finish_le_grids

  call close_module_test('le_grids')
  call finish_mp

end program test_le_grids
