module gs2_reinit
  implicit none

  private

  public :: reset_time_step, delt_adj
  public :: check_time_step, time_reinit
  public :: init_reinit, wnml_gs2_reinit
  public :: reduce_time_step, increase_time_step
  public :: init_gs2_reinit, finish_gs2_reinit



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
  !public :: reinit_gk_and_field_equations

  !> This function overrides the in_memory flag
  !! and should only be used if you know what
  !! you are doing.
  public :: gs2_reinit_unit_test_set_in_memory


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
    integer, intent(in) :: unit
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





  subroutine reset_time_step (current_init, istep, my_exit, job_id)
    use run_parameters, only: reset
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt_min
    use gs2_init, only: init_type, init
    use gs2_init, only: init_level_list
    use mp, only: proc0
    use file_utils, only: error_unit
    use job_manage, only: time_message
    use nonlinear_terms, only: gryfx_zonal
    implicit none
    integer, intent(in) :: istep 
    logical, intent(inout) :: my_exit
    integer, intent (in), optional :: job_id
    logical :: reset_in
    integer, save :: istep_last = -1 ! allow adjustment on first time step
    integer, save :: nconsec=0
    type(init_type), intent(inout) :: current_init

    real :: fac = 1.0

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

    !call save_fields_and_dist_fn

    gnew = 0.

    ! Move to the correct init level
    call init(current_init, init_level_list%override_timestep)
! change timestep 

    if(gryfx_zonal%on) then
!both code_dt = dt_gs2 and code_dt_cfl = dt_cfl_gryfx are in gs2 units
!we want to check if dt_gryfx = 2*dt_gs2 is too big/small when compared to
!dt_cfl_gryfx
      fac = 2.0
    endif

! If timestep is too big, make it smaller
    if (code_dt*fac > code_dt_cfl) then
       call reduce_time_step
! If timestep is too small, make it bigger
    else if (code_dt*fac <= min(dt0, code_dt_cfl/delt_adj/delt_cushion)) then
       call increase_time_step
    endif
    
    call save_dt (code_dt)

    if (proc0 .and. .not. present(job_id)) write(*,*) 'Changing time step to ', user_dt

    ! Don't reset antenna here because species parameters
    ! have not changed so resetting antenna would cause
    ! an unnecessary discontinuity
    !call reinit_gk_and_field_equations(reset_antenna=.false.)
    call init(current_init, init_level_list%full)
    
    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    istep_last = istep

    !Now re-enable reset so we leave it in the same state as on entering
    reset=reset_in

  end subroutine reset_time_step

  subroutine check_time_step (reset, exit)
    use gs2_time, only: code_dt_cfl, code_dt
    use nonlinear_terms, only: gryfx_zonal
    use mp, only: broadcast
    implicit none
    logical, intent(in) :: exit
    logical, intent(out) :: reset

    real :: fac = 1.0

    if (first) call init_reinit
    first = .false.
    reset = .false.

    if(gryfx_zonal%on) then
      !code_dt_cfl is only set on proc 0 by gryfx in nlps.cu
      call broadcast(code_dt_cfl)
!both code_dt = dt_gs2 and code_dt_cfl = dt_cfl_gryfx are in gs2 units
!we want to check if dt_gryfx = 2*dt_gs2 is too big/small when compared to
!dt_cfl_gryfx
      fac = 2.0
    endif

! nothing to do if exiting in this iteration
    if (exit) return

! If timestep is too big, make it smaller
    if (code_dt*fac > code_dt_cfl) reset = .true. !Note this logic is repeated in gs2_time::check_time_step_too_large
       
! If timestep is too small, make it bigger
    if (code_dt*fac <= min(dt0, code_dt_cfl/delt_adj/delt_cushion)) reset = .true.

  end subroutine check_time_step


  subroutine init_gs2_reinit
    call init_reinit
  end subroutine init_gs2_reinit

  subroutine finish_gs2_reinit
    first = .true.
  end subroutine finish_gs2_reinit

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

