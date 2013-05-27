

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



  call announce_module_test('analytical_falpha')

  call announce_test('is_converged')
  call process_test(unit_test_is_converged(), 'is_converged')
  call announce_test('chandrasekhar')
  call process_test(analytical_falpha_unit_test_chandrasekhar(), 'chandrasekhar')

  parameters%alpha_ion_collision_rate = 0.1
  parameters%alpha_electron_collision_rate = 0.1
  parameters%ion_vth = 1.0
  parameters%electron_vth = (1836.0)**0.5
  parameters%alpha_vth = (3.6e2 / 4.0)**0.5
  parameters%source = 0.01
  parameters%alpha_injection_energy = 3.6e2
  parameters%ion_temp = 1.0
  parameters%energy_0 = 0.01
  parameters%alpha_charge = 2.0
  parameters%electron_charge = -1.0
  parameters%ion_charge = 1.0
  parameters%alpha_mass = 4.0
  parameters%ion_mass = 1.0
  parameters%electron_mass = 1.0/1836.0
  parameters%source_prim = 1.6

  call announce_test('nu_parallel')
  call process_test(analytical_falpha_unit_test_nu_parallel( &
    parameters, 0.8, 7.07304892219606e-7, eps), 'nu_parallel')

  
  call announce_test('falpha_integrand')
  call process_test(analytical_falpha_unit_test_falpha_integrand(&
    parameters, 0.9,  0.4, 2.40051081062416e-72, eps), 'falpha_integrand 0.9')
  call process_test(analytical_falpha_unit_test_falpha_integrand(&
    parameters, 0.01, 0.4, 3.37485927208412e67, eps), 'falpha_integrand 0.01')
  call process_test(analytical_falpha_unit_test_falpha_integrand(&
    parameters, 2.5,  0.4, 0.0, eps), 'falpha_integrand 2.5')

  call announce_test('falpha')
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 0.02, 0.01, 128, 10.6555616529750, eps), 'falpha 0.02') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 0.5, 0.01, 2048*8, 6.06172663857817, eps), 'falpha 0.5') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 1.0, 0.01, 2048*16, 5.07680938596173, eps), 'falpha 1.0') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 2.5, 0.01, 256, 0.0, eps), 'falpha 2.5') 
  call process_test(analytical_falpha_unit_test_falpha(&
    parameters, 0.01, 0.01, 128, 0.0, eps), 'falpha 0.01') 


  allocate(array(6, 1,  7))

  array(:,1,1) = (/0.0200000000000000, 0.266666666666667, 0.513333333333333, 0.760000000000000, 1.00666666666667, 1.25333333333333/) ! egrid
  array(:,1,5) = (/10.6555616529750, 6.48908209538287, 6.03563607197710, 5.54447700644160, 0.460557756734939, 1.25295106064775e-39/)! f0_rslt
!  array(:,1,6) = (/76.1206676344528, 0.677300806262160, -1.99999999999978, -2.00000000000422, -2.00000000000422, -2.00000000000644/) ! f0_rslt
! The array below was generated using my own Mathematica worksheet, which also agrees perfectly with the above results (GW)
! Amended again after major changes to normalisations, EGH
  array(:,1,6) = (/12.5256589692896, 1416.58825193068, 1108.25154924379, 1002.67534577706, 1.00000000000000, 1.00000000000000/) 
  array(:,1,7) = -parameters%source_prim

  parameters%negrid = 6
  parameters%alpha_is = 1
  call announce_test('calculate_arrays')
  call process_test(analytical_falpha_unit_test_calculate_arrays(&
    parameters, array(:,:,1), array(:,:,2), array(:,:,3), array(:,:,4),&
    array(:,:,5), array(:,:,6), array(:,:,7), eps), 'calculate_arrays') 


  call close_module_test('analytical_falpha')

end program test_analytical_falpha
