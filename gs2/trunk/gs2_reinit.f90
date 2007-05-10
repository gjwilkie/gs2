module gs2_reinit
  implicit none

  public :: reset_time_step
  public :: check_time_step

  real :: delt_adj, dt0
  real :: delt_cushion = 1.5
  real :: delt_minimum 
  real, save :: time_nc = 0., time_reinit=0.

contains

  subroutine reset_time_step (istep, delt_cfl, exit)

    use additional_linear_terms, only: alt_reset => reset_init
!    use additional_terms, only: at_reset => reset_init
    use collisions, only: c_reset => reset_init
    use dist_fn, only: d_reset => reset_init
    use fields, only: f_reset => reset_init, init_fields
    use fields_implicit, only: fi_reset => reset_init
    use fields_explicit, only: fe_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use init_g, only: g_reset => reset_init
    use run_parameters, only: delt, tnorm, fphi, fapar, faperp
    use gs2_save, only: gs2_save_for_restart
    use dist_fn_arrays, only: gnew
    use gs2_time, only: stime, simdt
    use nonlinear_terms, only: nl_reset => reset_init
    use mp, only: proc0
    use file_utils, only: error_unit
    use antenna, only: dump_ant_amp

    logical :: exit
    real :: delt_cfl, simtime, deltsave
    integer :: istep 
    integer, save :: istep_last = -1 ! allow adjustment on first time step
    integer :: istatus

! save fields and distribution function

! calls on consecutive time steps is probably an error
    if (istep_last + 1 == istep) then
       exit = .true.
       if (proc0) write(error_unit(), *) 'Time step changing rapidly.  Abort run.'
       return
    end if

    if (delt/delt_adj <= delt_minimum*tnorm) then
       delt = delt_minimum*tnorm  ! set it so restart is ok
       exit = .true.
       return
    end if

    simtime = stime()/tnorm
    deltsave = simdt()/tnorm

    if (proc0) call dump_ant_amp
    call gs2_save_for_restart (gnew, simtime, deltsave, istatus, &
         fphi, fapar, faperp)

    gnew = 0.

! change timestep 

! If timestep is too big, make it smaller
    if (delt > delt_cfl) then
       delt = delt/delt_adj
! If timestep is too small, make it bigger
    else if (delt < min(dt0, delt_cfl/delt_adj/delt_cushion)) then
       delt = min(delt*delt_adj, dt0)
    endif
    
    if (proc0) write(*,*) 'Changing time step to ',delt/tnorm
    
! prepare to reinitialize inversion matrix, etc.
    call alt_reset
!    call at_reset
    call d_reset
    call c_reset
    call f_reset
    call fi_reset
    call fe_reset
    call ft_reset
    call g_reset
    call nl_reset

! reinitialize
    call init_fields
    if (proc0) call time_message(.true.,.false.,time_reinit,' Re-initialize')

    istep_last = istep

  end subroutine reset_time_step

  subroutine check_time_step (istep, delt, delt_cfl, reset, exit)

    integer :: istep
    real :: delt, delt_cfl
    logical :: reset, exit
    logical :: first = .true.

    if (first) call init_reinit
    first = .false.
    reset = .false.

! nothing to do if exiting in this iteration
    if (exit) return

! If timestep is too big, make it smaller
    if (delt > delt_cfl) reset = .true.
       
! If timestep is too small, make it bigger
    if (delt < min(dt0, delt_cfl/delt_adj/delt_cushion)) reset = .true.

! other choices
!     if (mod(istep,200) == 0) reset = .true.
!     if (delt > delt_cfl) exit = .true.

  end subroutine check_time_step

  subroutine init_reinit

    use run_parameters, only: delt_max, tnorm
    use mp, only: proc0, broadcast
    use file_utils, only: input_unit, input_unit_exist
    integer in_file
    logical exist
    namelist /reinit_knobs/ delt_adj, delt_minimum
    
    if (proc0) then
       dt0 = delt_max*tnorm
       delt_adj = 2.0
       delt_minimum = 1.e-5
       in_file = input_unit_exist("reinit_knobs",exist)
       if(exist) read (unit=in_file, nml=reinit_knobs)
    endif

    call broadcast (dt0)
    call broadcast (delt_adj)
    call broadcast (delt_minimum)

  end subroutine init_reinit

  subroutine time_message(lprint,lavg,tsum,chmessage)

    character(len=*), intent(in) :: chmessage
    logical, intent(in) :: lprint,lavg
    real, intent(in out) :: tsum
    character (len=10) :: zdate, ztime, zzone
    integer, dimension(8) :: ival
    logical :: ilprint
    real, save :: told=0., tnew=0., tavg=0.

    call date_and_time (zdate, ztime, zzone, ival)
    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
    ilprint=lprint
! ignore the error in told from a change of day in the next test:
    if (lavg .and. tnew-told < 1.2*tavg) ilprint=.false.
    if (told <= 0.) then
       ilprint=.false. ! first pass through this routine.
    else
       if (told > tnew) told=told-86400. ! evidently changed day
       if (lavg) then
! a single running average is maintained:
          tavg=0.9*tavg+0.1*(tnew-told)
       endif
       tsum=tsum+tnew-told
    endif
    if (ilprint)  then
       print *, chmessage,': ',tnew-told,' seconds'
    end if
    told = tnew

  end subroutine time_message


end module gs2_reinit

