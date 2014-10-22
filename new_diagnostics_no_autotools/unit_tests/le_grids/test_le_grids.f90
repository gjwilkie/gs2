
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
  !use general_f0, only: init_general_f0
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

  if (precision(eps).lt. 11) eps = eps * 100.0
  !if (precision(eps).lt. 11) write (*,*) 'low precision'

  ! Set up depenencies
  call init_mp
  if (proc0) call init_file_utils(dummy)
  call init_species
  !call init_general_f0



  call announce_module_test('le_grids')

  call announce_test('init_le_grids')

  allocate(sizes(1))
  allocate(energy_results(8,2,2)) ! negrid, nspec, nresults 
  sizes(1) = 24 ! Size of energy grid
  sizes(1) = 8 ! Size of energy grid

  ! Resuls for energy grids
  energy_results(:,1,1) = (/4.0468821658005425E-003,  0.10438457502758539, &
    0.55159372522148142, 1.5625000000000000, 3.0881259213302132, 4.7389544850238012, &
    5.9359713343080411, 7.2500000000000000/) ! Quadrature grid for Maxwellian... this was not calculated inpendently but was taken from the code output
  energy_results(:,2,1) = energy_results(:,1,1)
  !energy_results(:,3,1) = (/1.5104763817839494E-002, &
  !4.6790434219629591E-002, 0.13496048316473824, 0.30250000000000005, &
  !0.53674718302836155, 0.78090630796303029, 0.95472159703717052, 2.0000000000000000/) 
  ! Quadrature grid for alphas... this was not calculated inpendently but was taken from the code output

  !energy_min = 0.1 ! (the default value)
  !vcut_local = 2.5 !(the default value)
  !scal =  (energy_min - 1.0) / (energy_results(1,1,1) - vcut_local)
  !shift = energy_min - energy_results(1,1,1) * scal
  !energy_results(:,3,1) = energy_results(:,1,1)*scal + shift
  !write (*,*) 'aegrid', energy_results(:,3,1)


  ! Energy integration weights
  energy_results(:,1,2) = (/4.1155680601883985E-003, 0.22931213960689217, &
  1.6541668127613292, 5.1291308630037440, 9.2609382215892992, 10.410539988051941, &
  6.0367198818202885, 22.993938099111912 /) ! Weights grid for Maxwellian ions... this was not calculated inpendently but was taken from the code output
  ! Add maxwellian factor
  energy_results(:,1,2) = energy_results(:,1,2) * exp(-energy_results(:,1,1)) /(2.0*pi**1.5)
  energy_results(:,2,2) = energy_results(:,1,2) !  Electron weights
  ! For alphas,  mult weights taken from quadrature by alpha f0 calculated in sage
  !energy_results(:,3,2) = (/5.5300068428637454E-003, 3.7004138295655242E-002, &
  !0.14570320717563151, 0.35747990462790902, 0.57947173999263779, &
  !0.61757847512538577, 0.34953323523071972, 12.077007956766620/) * &
    !(/11.4454095380492, 7.12710579325183, 6.70488115401418, 6.49924424567821, 6.13726132537038, 5.72552031567843, 5.43111188929011, 2.41424807372596e-156/) / 2.0


  call process_test(le_grids_unit_test_init_le_grids(sizes, energy_results, eps), &
    'init_le_grids')
  !call init_le_grids(dummy, dummy)

  deallocate(sizes)


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
  species_weights(1) = 1.0
  g = cmplx(1.0, 1.0)
  integrate_species_results = cmplx(1.0,1.0)
  call announce_test('integrate species ions only g = 1')
  call process_test(&
    le_grids_unit_test_integrate_species(&
      g, &
      species_weights, &
      sizes, &
      integrate_species_results, &
      eps*100.0),&
    'integrate species ions only g = 1')

  ! Integrate a polynomial: yes it really does only get it to 
  ! 10% accuracy
  call announce_test(&
    'integrate species ions only g = energy^4.0 * 3.0*energy^1.25')
  g = cmplx(0.0, 0.0)
  do iglo = g_lo%llim_proc,g_lo%ulim_proc
    ie = ie_idx(g_lo, iglo)
    ene = energy_results(ie,1,1)
    !ene = energy_grid(ie,1)
    g(:,:,iglo) = cmplx(ene**4.0 + 3.0*ene**1.25,1.0)
  end do
  !write (*,*) 'g is ', size(g)
  integrate_species_results = cmplx(64.5070177949108,1.0)
  call process_test(&
    le_grids_unit_test_integrate_species(&
      g, &
      species_weights, &
      sizes, &
      integrate_species_results, &
      0.1),&
    'integrate species ions only g = energy^4.0 * 3.0*energy^1.25')

  !g = cmplx(1.0, 1.0)
  !species_weights = 0.0
  !species_weights(3) = 1.0
  !integrate_species_results = cmplx(25.3229494395902,25.3229494395902)
  !call announce_test('integrate species alphas only g = 1')
  !call process_test(&
    !le_grids_unit_test_integrate_species(&
      !g, &
      !species_weights, &
      !sizes, &
      !integrate_species_results, &
      !0.01),&
    !'integrate species alphas only  g = 1')


  deallocate(energy_results)
  deallocate(integrate_species_results)
  deallocate(species_weights)
  call finish_le_grids

  call close_module_test('le_grids')
  call finish_mp

end program test_le_grids
