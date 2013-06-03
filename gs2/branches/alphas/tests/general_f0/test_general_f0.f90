
!> A program that runs unit tests on the general_f0 module.
!! The test results were calculated using sage and are viewable at
!! http://www.sagenb.org/home/pub/5036
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_general_f0
  use unit_tests
  use general_f0
  use mp, only: init_mp, finish_mp
  use file_utils, only: init_file_utils
  use species, only: init_species, nspec
  use constants, only: pi
  implicit none
  real :: eps
  logical :: dummy

  real, dimension(:,:), allocatable :: epoints
  real, dimension(:,:), allocatable :: weights
  real, dimension(:,:,:), allocatable :: rslts
  integer :: negrid

  real :: vcut

  ! General config
  eps = 1.0e-7

  ! Set up depenencies
  call init_mp
  call init_file_utils(dummy)
  call init_species


  call announce_module_test('general_f0')

  call announce_test('init_general_f0')
  call process_test(general_f0_unit_test_init_general_f0(), 'init_general_f0')

  negrid = 4
  vcut = 2.5

  allocate(epoints(negrid, nspec))
  allocate(weights(negrid, nspec))
  allocate(rslts(negrid, nspec, 3))
  epoints(:,1) = (/0.0, 1.0, 2.5, 4.0/)
  rslts(:,1,1) = exp(-epoints(:,1))/(2.0*pi**1.5) ! ion f0
  rslts(:,2,1) = exp(-epoints(:,1))/(2.0*pi**1.5) ! electron f0
  rslts(:,3,1) = (/6.74132595300733, 6.26167292325401, 5.35544272159621, 2.00423248388538e-84/)

  rslts(:,1,2) = 1.0 ! ion temp
  rslts(:,2,2) = 0.8 ! electron temp
  rslts(:,3,2) = (/2023.61241984888, 1415.26340068638, 1158.15661293921, 1.00000000000000/) ! alpha temp

  rslts(:,1,3) = -(3.0 + (epoints(:,1) - 1.5) * 4.0) ! ion f0prim
  rslts(:,2,3) = -(3.0 + (epoints(:,1) - 1.5) * 6.0) ! electron f0prim
  rslts(:,3,3) = (/1.32051044369054, 0.856175603678116, -0.0397244931703406, -777.639723666515/) ! alpha f0prim

  call announce_test('calculate_f0_arrays')
  call process_test(&
     general_f0_unit_test_calculate_f0_arrays(epoints, weights, vcut, rslts, eps),&
     'calculate_f0_arrays')


  call finish_general_f0

  call close_module_test('general_f0')
  call finish_mp

end program test_general_f0
