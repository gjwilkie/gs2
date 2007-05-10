module gs2_reinit
  implicit none

  public :: reset_time_step

contains

  subroutine reset_time_step (istep, dt_cfl)

    use collisions, only: c_reset => reset_init
    use dist_fn, only: d_reset => reset_init
    use fields, only: f_reset => reset_init, init_fields
    use fields_implicit, only: fi_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use init_g, only: g_reset => reset_init
    use run_parameters, only: delt
    use gs2_save, only: gs2_save_for_restart
    use dist_fn, only: tstart
    use dist_fn_arrays, only: gnew
    use fields_arrays, only: aminv
    use gs2_time, only: stime
    use mp, only: iproc

    real :: dt_cfl, simtime
    integer :: istep
    integer :: istatus

! save fields and distribution function

    simtime = stime()
    call gs2_save_for_restart (gnew, aminv, simtime, istatus)

    gnew = 0.
    aminv = 0.

! reduce timestep 

    delt = delt/2.
    
! prepare to reinitialize inversion matrix, etc.
    call c_reset
    call f_reset
    call fi_reset
    call ft_reset
    call g_reset

! reinitialize
    call init_fields

  end subroutine reset_time_step

end module gs2_reinit
