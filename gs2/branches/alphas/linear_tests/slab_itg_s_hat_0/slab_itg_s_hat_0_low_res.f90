!> A test program to run a linear slab itg benchmark
!! at low resolution and check the growth rate

program slab_itg_s_hat_0_low_res

! make_lib is a compiler flag used if running with 
! trinity (coupled flux tube code)

  use gs2_main, only: run_gs2, functional_test_flag, ilast_step, finish_gs2
  use functional_tests
  use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
  use gs2_diagnostics, only: finish_gs2_diagnostics

  implicit none

  call init_mp

  test_driver_flag = .true.
  functional_test_flag = .true.

  if (proc0) call announce_functional_test('Linear slab ITG with zero magnetic shear')

  call run_gs2(mp_comm)

  ! We only expect this low res test to get the growth rate to 10%
  call test_growth_rate(ilast_step-1, (/0.164/), 0.1)

  call finish_gs2_diagnostics(ilast_step)
  call finish_gs2

  if (proc0) call close_functional_test('Linear slab ITG with zero magnetic shear')

  call finish_mp

end program slab_itg_s_hat_0_low_res
