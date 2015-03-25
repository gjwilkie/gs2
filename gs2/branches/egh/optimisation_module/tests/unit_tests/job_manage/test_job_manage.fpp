
!> This unit runs a tiny grid with nstep=10000, and checks
!! that both the .stop and the avail_cpu_time mechanisms
!! are working
!
program test_job_manage
  !use functional_tests
  !use checks_mod
  !call test_gs2('Linear CBC (unit test) to test new diagnostics', checks)
    use gs2_main, only: run_gs2, finish_gs2, trin_finish_gs2
    !use gs2_main, only: run_gs2, trin_finish_gs2
    use gs2_main, only: old_iface_state, finalize_gs2, finalize_diagnostics
    use gs2_main, only: reset_equations, initialize_diagnostics
    use gs2_main, only: finalize_equations
    use unit_tests
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
    use gs2_diagnostics, only: finish_gs2_diagnostics
    use job_manage, only: timer_local
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
#endif 
    implicit none
    real :: eps
    real :: start_time
    real :: avail_time = 2.0

    eps = 1.0e-7
    if (precision(eps).lt. 11) eps = eps * 1000.0


    call init_mp

    test_driver_flag = .true.
    functional_test_flag = .true.

   call announce_module_test("job_manage")

   if (proc0) call system("touch test_job_manage.stop")
   call announce_test("that gs2 doesn't run without limit when run_name.stop is present")
   start_time = timer_local()
    call run_gs2(mp_comm)
   call process_test((timer_local() - start_time < avail_time), &
     "that gs2 doesn't run without limit when run_name.stop is present")
    if (proc0) call system("rm test_job_manage.stop")
    !write (*,*) 'calling finish_gs2_diagnostics'
    !call finish_gs2_diagnostics(ilast_step)
#ifdef NEW_DIAG
    !call finish_gs2_diagnostics_new
#endif 
    call reset_equations(old_iface_state)

    !write (*,*) 'calling trin_finish_gs2'
    !call trin_finish_gs2

   call announce_test("that gs2 doesn't run without limit when avail_cpu_time is set")
   start_time = timer_local()
    call run_gs2(mp_comm)
    !write (*,*) 'start_time', start_time, timer_local()
   call process_test(((timer_local() - start_time .gt. avail_time - 1.0) ), &
     "that gs2 doesn't run without limit when avail_cpu_time is set")


!    call finish_gs2_diagnostics(ilast_step)
!#ifdef NEW_DIAG
!    call finish_gs2_diagnostics_new
!#endif 
!
!   call finish_gs2

    !call trin_finish_gs2

    call finalize_diagnostics(old_iface_state)
    call finalize_equations(old_iface_state)
    call finalize_gs2(old_iface_state)



   call close_module_test("job_manage")

    call finish_mp

contains
  


end program test_job_manage
