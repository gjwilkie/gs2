
!> A program to test the replay feature
!! (see the documentation for replay in 
!! gs2_main).
!!
!! This is free software release under the 
!! MIT license. 
!! Written by:
!!     Edmund Highcock (edmundhighcock@users.sourceforge.net)

program cyclone_itg_replay
  use gs2_main
  implicit none
  type(gs2_program_state_type) :: state
  !type(optimisation_type) :: optim
  call initialize_wall_clock_timer

  ! This is what we are testing
  state%replay = .true.

  call initialize_gs2(state)
  call initialize_equations(state)
  call initialize_diagnostics(state)
  state%print_times = .false.
  call evolve_equations(state, state%nstep)
  call finalize_diagnostics(state)
  call finalize_equations(state)
  state%print_times = .true.
  state%print_full_timers = .true.
  call finalize_gs2(state)



end program cyclone_itg_replay
