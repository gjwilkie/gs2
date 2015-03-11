module gs2_reinit
  implicit none

  public :: reset_time_step, delt_adj
  public :: check_time_step, time_reinit
  public :: init_reinit, wnml_gs2_reinit
  public :: reduce_time_step, increase_time_step



  !> This subroutine reinitializes the modules which are responsible for
  !! solving the linear, collisonal and nonlinear parts of the GK eqn,
  !! as well as the modules which solve the field equations. This is 
  !! typically necessary because some parameter (e.g. the timestep) which
  !! is used to calculate the value of cached arrays (e.g. the response
  !! matrix) in those modules has changed. Note that this will not cause
  !! dist_fn to reread its namelist.. thus you can change g_exb (and any
  !! other physics parameter in any other module, of course) 
  !! before calling this function 
  !! and expect your changes to be preserved. The function 
  !! save_fields_and_dist_fn must be called before this, or 
  !! otherwise the current field and dist_fn values will be lost. The
  !! logical flag in_memory must be given the same value that was 
  !! set in save_fields_and_dist_fn. 
  public :: reinit_gk_and_field_equations

  !> This function overrides the in_memory flag
  !! and should only be used if you know what
  !! you are doing.
  public :: gs2_reinit_unit_test_set_in_memory

  private

  real :: delt_adj, dt0
  real :: delt_cushion
  real :: delt_minimum 
  real :: time_reinit(2)=0.
  logical :: abort_rapid_time_step_change
  logical :: first=.true.
  logical :: in_memory

contains

  subroutine gs2_reinit_unit_test_set_in_memory(in_memory_in)
    use gs2_init, only: gs2_init_in_memory=>in_memory
    logical, intent(in) :: in_memory_in
    in_memory = in_memory_in
    gs2_init_in_memory = in_memory
  end subroutine gs2_reinit_unit_test_set_in_memory

  subroutine wnml_gs2_reinit(unit)
    implicit none
    integer :: unit
    write (unit, *)
    write (unit, fmt="(' &',a)") "reinit_knobs"
    write (unit, fmt="(' delt_adj = ',e17.10)") delt_adj
    write (unit, fmt="(' delt_minimum = ',e17.10)") delt_minimum
    write (unit, fmt="(' /')")       
  end subroutine wnml_gs2_reinit

  subroutine reduce_time_step
    use gs2_time, only: code_dt
    implicit none
    if (first) call init_reinit
    code_dt = code_dt/delt_adj
  end subroutine reduce_time_step

  subroutine increase_time_step
    use gs2_time, only: code_dt
    implicit none
    if (first) call init_reinit
    code_dt = min(code_dt*delt_adj, dt0)
  end subroutine increase_time_step


  subroutine reinit_gk_and_field_equations(reset_antenna)
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: g_restart_tmp
    use fields_arrays, only: phinew, aparnew, bparnew
    use collisions, only: c_reset => reset_init
    use dist_fn, only: d_reset => reset_init
    use fields, only: f_reset => finish_fields, init_fields
    use init_g, only: g_reset => reset_init
    use nonlinear_terms, only: nl_reset => reset_init
    use antenna, only: a_reset => reset_init
    use gs2_init, only: load_saved_field_values
    logical, intent(in) :: reset_antenna
! prepare to reinitialize inversion matrix, etc.
    call d_reset
    call c_reset
    call f_reset
    call g_reset(.not.in_memory)
    call nl_reset

    if (reset_antenna) call a_reset

    write (*,*) 'EEEETTTT'

! reinitialize
    call init_fields

    write (*,*) 'EEEEFFF'

!Update fields if done in memory
!Don't need/want to update if force_maxwell_reinit
    call load_saved_field_values
  end subroutine reinit_gk_and_field_equations



  subroutine reset_time_step (current_init, istep, my_exit, job_id)
    use run_parameters, only: reset
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt_min
    use gs2_init, only: save_fields_and_dist_fn, init_level_type, init
    use gs2_init, only: init_level_list
    use mp, only: proc0
    use file_utils, only: error_unit
    use job_manage, only: time_message
    implicit none
    logical, intent(inout) :: my_exit
    logical :: reset_in
    integer :: istep 
    integer, save :: istep_last = -1 ! allow adjustment on first time step
    integer, save :: nconsec=0
    integer, intent (in), optional :: job_id
    type(init_level_type), intent(inout) :: current_init


    if (first) call init_reinit
    first = .false.

! save fields and distribution function

! calls on consecutive time steps is probably an error
    if (istep_last + 1 == istep) then
       nconsec=nconsec+1
    else
       nconsec=0
    endif

    if (nconsec .gt. 4 .and. abort_rapid_time_step_change) then
       my_exit = .true.
       if (proc0) write(error_unit(), *) 'Time step changing rapidly.  Abort run.'
       return
    end if

    if (code_dt/delt_adj <= code_dt_min) then
       code_dt = code_dt_min  ! set it so restart is ok
       my_exit = .true.
       if (proc0) write(error_unit(), *) 'Time step wants to fall below delt_min.  Abort run.'
       return
    end if

    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    !First disable the reset flag so we can call 
    !routines needed in reinit
    reset_in=reset
    reset=.false.

    call save_fields_and_dist_fn

    gnew = 0.

    ! Move to the correct init level
    call init(current_init, init_level_list%override_timestep)
! change timestep 

! If timestep is too big, make it smaller
    if (code_dt > code_dt_cfl) then
       call reduce_time_step
! If timestep is too small, make it bigger
    else if (code_dt < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) then
       call increase_time_step
    endif
    
    call save_dt (code_dt)

    if (proc0 .and. .not. present(job_id)) write(*,*) 'Changing time step to ', user_dt

    ! Don't reset antenna here because species parameters
    ! have not changed so resetting antenna would cause
    ! an unnecessary discontinuity
    !call reinit_gk_and_field_equations(reset_antenna=.false.)
    call init(current_init, init_level_list%fields)
    
    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    istep_last = istep

    !Now re-enable reset so we leave it in the same state as on entering
    reset=reset_in

  end subroutine reset_time_step

  subroutine check_time_step (reset, exit)
    use gs2_time, only: code_dt_cfl, code_dt
    implicit none
    logical, intent(in) :: exit
    logical, intent(out) :: reset

    if (first) call init_reinit
    first = .false.
    reset = .false.

! nothing to do if exiting in this iteration
    if (exit) return

! If timestep is too big, make it smaller
    if (code_dt > code_dt_cfl) reset = .true. !Note this logic is repeated in gs2_time::check_time_step_too_large
       
! If timestep is too small, make it bigger
    if (code_dt < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) reset = .true.

  end subroutine check_time_step

  subroutine init_reinit
    use run_parameters, only: code_delt_max
    use mp, only: proc0, broadcast
    use file_utils, only: input_unit, input_unit_exist
    use gs2_time, only: save_dt_min
    use gs2_init, only: gs2_init_in_memory=>in_memory
    implicit none
    integer :: in_file
    logical :: exist

    namelist /reinit_knobs/ delt_adj, delt_minimum, delt_cushion, &
                            abort_rapid_time_step_change, in_memory
    if(.not.first)return
    first=.false.

    if (proc0) then
       dt0 = code_delt_max
       delt_adj = 2.0
       delt_minimum = 1.e-5
       delt_cushion = 1.5
       abort_rapid_time_step_change = .true.
       in_memory=.false.
       in_file = input_unit_exist("reinit_knobs",exist)
       if(exist) read (unit=in_file, nml=reinit_knobs)
    endif

    call broadcast (dt0)
    call broadcast (delt_adj)
    call broadcast (delt_minimum)
    call broadcast (delt_cushion)
    call broadcast (abort_rapid_time_step_change)
    call broadcast (in_memory)
    call save_dt_min (delt_minimum)
    ! Override with the value from gs2 init
    in_memory = gs2_init_in_memory
  end subroutine init_reinit
end module gs2_reinit

