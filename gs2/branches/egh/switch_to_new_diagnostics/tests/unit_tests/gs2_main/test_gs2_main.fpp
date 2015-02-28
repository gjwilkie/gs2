
!> This unit tests the gs2 main interface,
!! and by repeatedly initializing and finalizing
!! gs2, tests that gs2 is being properly tidied up,
!! variables deallocated etc.
!
program test_gs2_main
  use gs2_main
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_comm
  implicit none
  real :: eps
  type(gs2_program_state_type) :: state

  eps = 1.0e-7
  if (precision(eps).lt. 11) eps = eps * 1000.0
  
  call init_mp
  
  call announce_module_test("gs2_main")

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call initialize_equations(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call initialize_equations(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  !!program gs2
    !!type(gs2_program_state_type) :: state
    call initialize_gs2(state)
    call initialize_equations(state)
    call initialize_diagnostics(state)
    !if (state%eigsolve) then 
      !call solve_eigenproblem
    !else
    call evolve_equations(state, state%nstep/2)
    call evolve_equations(state, state%nstep/2)
    !! This call should do nothing and print a warning
    call evolve_equations(state, state%nstep/2)
    call evolve_equations(state, state%nstep/2)

    call calculate_outputs(state)


    call finalize_diagnostics(state)
    call finalize_equations(state)
    call finalize_gs2(state)
  !!end program gs2


    !!call init_mp

    !!test_driver_flag = .true.
    !!functional_test_flag = .true.



  call finalize_overrides(state)
  call close_module_test("gs2_main")

  !call finish_gs2


  call finish_mp



contains
  


end program test_gs2_main
