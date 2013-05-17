
!> A program that runs unit tests on the analytical_falpha module.
!! The test results were calculated using sage and are viewable at
!! http://www.sagenb.org/home/pub/5036
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_analytical_falpha
  use unit_tests
  use analytical_falpha
  implicit none
  real :: eps
  type(analytical_falpha_parameters_type) :: parameters
  real, dimension(:,:,:), allocatable :: array

  ! General config
  eps = 1.0e-8


  write (*,*) '*************************'
  write (*,*) 'Testing analytical_falpha'
  write (*,*) '*************************'

  call announce_test('is_converged')
  call process_test(unit_test_is_converged(), 'is_converged')
  call announce_test('chandrasekhar')
  call process_test(analytical_falpha_unit_test_chandrasekhar(), 'chandrasekhar')

  parameters%alpha_ion_collision_rate = 0.1/4.0**2.0
  parameters%alpha_electron_collision_rate = &
     parameters%alpha_ion_collision_rate / 1836.0**2.0
  parameters%ion_vth = 1.0
  parameters%electron_vth = 42.0
  parameters%alpha_vth = (3.6e2 / 4.0)**0.5
  parameters%ion_temp = 1.0
  parameters%alpha_mass = 4.0
  parameters%source = 0.01
  parameters%alpha_injection_energy = 3.6e2

  call announce_test('nu_parallel')
  call process_test(analytical_falpha_unit_test_nu_parallel( &
    parameters, 0.8, 0.000121314067605073, eps), 'nu_parallel')

  
  call announce_test('falpha_integrand')
  call process_test(analytical_falpha_unit_test_falpha_integrand(&
    parameters, 0.4,  0.9, 74269.3906980909, eps), 'falpha_integrand 0.9')
  call process_test(analytical_falpha_unit_test_falpha_integrand(&
    parameters, 0.4,  0.01, 3428.53572237518, eps), 'falpha_integrand 0.01')
  call process_test(analytical_falpha_unit_test_falpha_integrand(&
    parameters, 0.4,  2.5, 3.03673247901918e6, eps), 'falpha_integrand 2.5')

  call announce_test('falpha')
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 0.02, 0.01, 128, 264.9390892941125*0.01, eps), 'falpha 0.02') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 0.5, 0.01, 256, 4571.551171450769*0.01, eps), 'falpha 0.5') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 1.0, 0.01, 256, 6233.062047471938*0.01, eps), 'falpha 1.0') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 2.5, 0.01, 256, 310.32588629862335*0.01, eps), 'falpha 2.5') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 0.01, 0.01, 128, 0.0, eps), 'falpha 0.01') 


  allocate(array(6, 1,  7))

  array(:,1,1) = (/0.0100000000000000, 0.675000000000000, 1.34000000000000, &
  2.00500000000000, 2.67000000000000, 3.33500000000000/) ! egrid
  array(:,1,5) = (/0.000000000000000, 53.4776774832068, 31.5777514771834, &
  8.35159722869193, 2.20880756252364, 0.584179374874592/) ! f0_rslt

  parameters%negrid = 6
  parameters%alpha_is = 1
  call announce_test('calculate_arrays')
  call process_test(analytical_falpha_unit_test_calculate_arrays(&
    parameters, array(:,:,1), array(:,:,2), array(:,:,3), array(:,:,4),&
    array(:,:,5), eps), 'calculate_arrays') 



end program test_analytical_falpha
