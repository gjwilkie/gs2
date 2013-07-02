module gs2_reinit
  implicit none

  public :: reset_time_step
  public :: check_time_step

  real :: delt_adj, dt0
!  real :: delt_cushion = 1.5
  real :: delt_cushion
  real :: delt_minimum 
  real, save :: time_reinit(2)=0.

contains

  subroutine reset_time_step (istep, exit)

    !+PJK  Added code to prevent problems if the explicit DG scheme's adaptive
    !      timestep algorithm is in use
    use fields, only: fieldopt_switch, fieldopt_implicit, fieldopt_explicit
    use dg_scheme, only: adaptive_dt
    !-PJK

    use collisions, only: c_reset => reset_init, vnmult
    use dist_fn, only: d_reset => reset_init
    use fields, only: f_reset => reset_init, init_fields
    use fields_implicit, only: fi_reset => reset_init
    use fields_explicit, only: fe_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use init_g, only: g_reset => reset_init
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt
    use gs2_save, only: gs2_save_for_restart
    use dist_fn_arrays, only: gnew
    use gs2_time, only: user_time, code_dt_min
    use nonlinear_terms, only: nl_reset => reset_init
    use mp, only: proc0
    use file_utils, only: error_unit
    use antenna, only: dump_ant_amp
    use job_manage, only: time_message

    logical :: exit
    integer :: istep 
    integer, save :: istep_last = -1 ! allow adjustment on first time step
    integer :: istatus
    integer, save :: nconsec=0

! save fields and distribution function

! calls on consecutive time steps is probably an error
    if (istep_last + 1 == istep) then
       nconsec=nconsec+1
    else
       nconsec=0
    endif

    !+PJK    if (nconsec .gt. 4) then
!!    if ( (.not.adaptive_dt).and.(nconsec .gt. 4) ) then
    if ( ((fieldopt_switch == fieldopt_implicit).or. &
         ((fieldopt_switch == fieldopt_explicit).and.(.not.adaptive_dt)) ).and. &
         (nconsec .gt. 4) ) then
       !-PJK
       exit = .true.
       if (proc0) write(error_unit(), *) 'Time step changing rapidly.  Abort run.'
       return
    end if

    if (code_dt/delt_adj <= code_dt_min) then
       code_dt = code_dt_min  ! set it so restart is ok
       exit = .true.
       return
    end if

    if (proc0) call time_message(.true.,time_reinit,' Re-initialize')

    if (proc0) call dump_ant_amp
    call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)

    gnew = 0.

! change timestep 

! If timestep is too big, make it smaller
    if (code_dt > code_dt_cfl) then
       code_dt = code_dt/delt_adj

! If timestep is too small, make it bigger
    !+PJK  else if (code_dt < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) then
!!    else if ( (.not.adaptive_dt).and. &
    else if ( ((fieldopt_switch == fieldopt_implicit).or. &
         ((fieldopt_switch == fieldopt_explicit).and.(.not.adaptive_dt)) ).and. &
         (code_dt < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) ) then
       !-PJK
       code_dt = min(code_dt*delt_adj, dt0)

    endif
    
    call save_dt (code_dt)

    if (proc0) write(*,*) 'Changing time step to ', user_dt
    
! prepare to reinitialize inversion matrix, etc.
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

    if (proc0) call time_message(.true.,time_reinit,' Re-initialize')

    istep_last = istep

  end subroutine reset_time_step

!  subroutine check_time_step (istep, reset, exit)
  subroutine check_time_step (reset, exit)

    use gs2_time, only: code_dt_cfl, code_dt

!    integer :: istep
    logical :: reset, exit
    logical :: first = .true.

    if (first) call init_reinit
    first = .false.
    reset = .false.

! nothing to do if exiting in this iteration
    if (exit) return

! If timestep is too big, make it smaller
    if (code_dt > code_dt_cfl) reset = .true.
       
! If timestep is too small, make it bigger
    if (code_dt < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) reset = .true.

! other choices
!     if (mod(istep,200) == 0) reset = .true.
!     if (code_dt > code_dt_cfl) exit = .true.

  end subroutine check_time_step

  subroutine init_reinit

    use run_parameters, only: code_delt_max
    use mp, only: proc0, broadcast
    use file_utils, only: input_unit, input_unit_exist
    use gs2_time, only: save_dt_min
    integer in_file
    logical exist
    namelist /reinit_knobs/ delt_adj, delt_minimum, delt_cushion
    
    if (proc0) then
       dt0 = code_delt_max
       delt_adj = 2.0
       delt_minimum = 1.e-5
       delt_cushion = 1.5
       in_file = input_unit_exist("reinit_knobs",exist)
       if(exist) read (unit=in_file, nml=reinit_knobs)
    endif

    call broadcast (dt0)
    call broadcast (delt_adj)
    call broadcast (delt_minimum)
    call broadcast (delt_cushion)

    call save_dt_min (delt_minimum)

  end subroutine init_reinit

end module gs2_reinit

