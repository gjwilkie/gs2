
!> This unit runs a tiny grid with nstep=10000, and checks
!! that both the .stop and the avail_cpu_time mechanisms
!! are working
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


    !call init_mp

    !test_driver_flag = .true.
    !functional_test_flag = .true.



  call close_module_test("gs2_main")

  call finish_mp

contains
  


end program test_gs2_main
