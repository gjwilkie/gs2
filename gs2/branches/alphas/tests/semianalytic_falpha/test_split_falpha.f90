


!> A program that runs unit tests on the analytical_falpha module.
!! The test results were calculated using sage and are viewable at
!! http://www.sagenb.org/home/pub/5036
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_semianalytic_falpha
  use unit_tests
  use semianalytic_falpha
  use mp, only: init_mp, finish_mp
  implicit none
  real :: eps

  type(semianalytic_falpha_parameters_type) :: parameters
  logical :: dummy

  real, dimension(:,:,:), allocatable :: array

  call init_mp

  ! General config
  eps = 1.0e-8



  call announce_module_test('semianalytic_falpha')

  call announce_test('is_converged')
  call process_test(unit_test_is_converged(), 'is_converged')
  call announce_test('chandrasekhar')
  call process_test(semianalytic_falpha_unit_test_chandrasekhar(), 'chandrasekhar')
  call announce_test('chandrasekhar prime')
  call process_test(&
     semianalytic_falpha_unit_test_chandrasekhar_prime(0.1, 0.369396267397740, eps), &
     'chandrasekhar prime')

  call announce_test('simpson')
  call process_test(semianalytic_falpha_unit_test_simpson(eps), 'simpson')

  parameters%energy_0 = 0.01
  parameters%source = 0.01
  parameters%source_prim = 0.0 !1.6

  parameters%alpha_ion_collision_rate = 0.1
  parameters%alpha_electron_collision_rate = 0.1

  parameters%alpha_vth = (3.6e2 / 4.0)**0.5
  parameters%alpha_injection_energy = 3.6e2
  parameters%alpha_mass = 4.0
  parameters%alpha_charge = 2.0
  !parameters%alpha_density = -1.0

  parameters%ion_temp = 1.0
  parameters%ion_tprim = 4.0
  parameters%ion_fprim = 3.0
  parameters%ion_mass = 1.0
  parameters%ion_vth = 1.0
  parameters%ion_charge = 1.0

  !parameters%electron_temp = 1.0
  parameters%electron_charge = -1.0
  parameters%electron_tprim = 6.0
  parameters%electron_fprim = 3.0
  parameters%electron_mass = 1.0/1836.0
  parameters%electron_vth = (1836.0)**0.5


  call announce_test('nu_parallel')
  call process_test(semianalytic_falpha_unit_test_nu_parallel( &
    parameters, 0.8, 7.07304892219606e-7, eps), 'nu_parallel')

  call announce_test('nu_parallel_prime')
  call process_test(semianalytic_falpha_unit_test_nu_parallel_prime( &
    parameters, 0.2, -0.000127328159032943, eps*1.0), 'nu_parallel_prime 0.2')
  call process_test(semianalytic_falpha_unit_test_nu_parallel_prime( &
    parameters, 0.8, -3.99791538405283e-6, eps*1.0), 'nu_parallel_prime 0.8')
  
  call announce_test('falpha_integrand')
  call process_test(semianalytic_falpha_unit_test_falpha_integrand(&
    parameters, 0.9,  0.4, 2.40051081062416e-72, eps), 'falpha_integrand 0.9')
  call process_test(semianalytic_falpha_unit_test_falpha_integrand(&
    parameters, 0.01, 0.4, 3.37485927208412e67, eps), 'falpha_integrand 0.01')
  call process_test(semianalytic_falpha_unit_test_falpha_integrand(&
    parameters, 2.5,  0.4, 0.0, eps), 'falpha_integrand 2.5')

  call announce_test('falpha')
  call process_test(semianalytic_falpha_unit_test_falpha(&
    parameters, 0.02, 0.01, 128, 10.6555616529750, eps), 'falpha 0.02') 
  call process_test(semianalytic_falpha_unit_test_falpha(&
    parameters, 0.5, 0.01, 2048*8, 6.06172663857817, eps), 'falpha 0.5') 
  call process_test(semianalytic_falpha_unit_test_falpha(&
    parameters, 1.0, 0.01, 2048*16, 5.07680938596173, eps), 'falpha 1.0') 
  call process_test(semianalytic_falpha_unit_test_falpha(&
    parameters, 2.5, 0.01, 256, 0.0, eps), 'falpha 2.5') 
  call process_test(semianalytic_falpha_unit_test_falpha(&
    parameters, 0.01, 0.01, 128, 0.0, eps), 'falpha 0.01') 


  call announce_test('dfalpha_dti')
  call process_test(semianalytic_falpha_unit_test_dfalpha_dti(&
    parameters, 0.5,  2048*8, -4.00357641291620, eps), 'dfalpha_dti 0.5')
  call process_test(semianalytic_falpha_unit_test_dfalpha_dti(&
    parameters, 0.9,  2048*16, -4.00410087513562, eps), 'dfalpha_dti 0.9')

  call announce_test('dfalpha_dnupar')
  call process_test(semianalytic_falpha_unit_test_dfalpha_dnupar(&
    parameters, 0.1,  2048*32, 6.90802623886552, eps*1000.0), 'dfalpha_dnupar 0.1')
  call process_test(semianalytic_falpha_unit_test_dfalpha_dnupar(&
    parameters, 0.5,  2048*8, 6.25450124922097 , eps*1000.0), 'dfalpha_dnupar 0.5')

  allocate(array(6, 1,  7))

!  array(:,1,1) = (/0.0200000000000000, 0.266666666666667, 0.513333333333333, 0.760000000000000, 1.00666666666667, 1.25333333333333/) ! egrid
  array(:,1,1) = (/0.0200000000000000, 0.266666666666667, 0.513333333333333, 0.760000000000000, 1.00666666666667, 1.25333333333333/) ! egrid
!  array(:,1,5) = (/10.6555616529750, 6.48908209538287, 6.03563607197710, 5.54447700644160, 0.460557756734939, 1.25295106064775e-39/)! f0_rslt
  array(:,1,5) = (/22.27137676396786, 12.97816419065629, 12.071217214408448800, 11.0889540128830400, -4.5772152276611050e-148, -1.9892625505801470e-186/)! f0_rslt
! Amended again after major changes to normalisations, EGH
!  array(:,1,6) = (/12.5256589692896, 1416.58825193068, 1108.25154924379, 1002.67534577706, 1.00000000000000, 1.00000000000000/) 
  array(:,1,6) = (/12.5256589692896, 1416.58825193068, 1108.25154924379, 1002.67534577706, 1.00000000000000, 1.00000000000000/) 
  array(:,1,7) = (/0.0484936118028753, 2.68468051893912, 2.22466808876727, 1.73271163434981, -8.32950066620934, -363.529500481621/)
  parameters%negrid = 6
  parameters%alpha_is = 1

  !parameters%ion_tprim = 0.0 !1200.0
  !parameters%electron_tprim = 0.0! 80.0
  !parameters%source_prim = 0.0!20.8
  !parameters%electron_fprim = 40000.8
  !parameters%ion_fprim = parameters%electron_fprim

  call announce_test('calculate_arrays')
  call process_test(semianalytic_falpha_unit_test_calculate_arrays(&
    parameters, array(:,:,1), array(:,:,2), array(:,:,3), array(:,:,4),&
    array(:,:,5), array(:,:,6), array(:,:,7), eps), 'calculate_arrays') 


  call announce_test('alpha_density')
  call process_test(semianalytic_falpha_unit_test_alpha_density(&
    parameters, 24.6272836977552, 0.001), 'alpha_density')
  parameters%alpha_density = 0.5
  parameters%source = -1
  dummy = semianalytic_falpha_unit_test_calculate_arrays(&
    parameters, array(:,:,1), array(:,:,2), array(:,:,3), array(:,:,4),&
    array(:,:,5), array(:,:,6), array(:,:,7), eps)
  call announce_test('alpha_density adjusted')
  call process_test(semianalytic_falpha_unit_test_alpha_density(&
    parameters, 1.0, eps), 'alpha_density adjusted')



  call close_module_test('semianalytic_falpha')

  call finish_mp

end program test_semianalytic_falpha
