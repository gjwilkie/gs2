
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
    parameters, 0.02, 0.01, 16, 63.562943357601526, eps), 'falpha 0.02') 






end program test_analytical_falpha
