
!> A program that runs unit tests on the analytical_falpha module.
!! The test results were calculated using sage and are viewable at
!! http://www.sagenb.org/home/pub/5036
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@sourceforge.net)
program test_general_f0
  use unit_tests
  use general_f0
  implicit none
  real :: eps

  ! General config
  eps = 1.0e-8



  call announce_module_test('general_f0')

  !call announce_test('is_converged')
  !call process_test(unit_test_is_converged(), 'is_converged')
  call close_module_test('general_f0')

end program test_general_f0
