
!> This unit tests the gs2 main interface,
!! and by repeatedly initializing and finalizing
!! gs2, tests that gs2 is being properly tidied up,
!! variables deallocated etc.
!
program test_optimisation
  use gs2_optimisation
  use gs2_main, only: gs2_program_state_type
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_comm
  implicit none
  real :: eps
  type(gs2_program_state_type) :: state
  !type(optimisation_test_type) :: optest

  eps = 1.0e-7
  if (precision(eps).lt. 11) eps = eps * 1000.0
  
  call init_mp
  
  call announce_module_test("optimisation")

  call initialize_gs2_optimisation(state)

  !call measure_timestep(state)

  call announce_test('optimise gs2')
  call optimise_gs2(state)

  call close_module_test("optimisation")

  call finish_mp



contains
  


end program test_optimisation
