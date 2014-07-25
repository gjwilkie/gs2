module gs2_reinit
  implicit none

  public :: reset_time_step
  public :: check_time_step
  public :: init_reinit, wnml_gs2_reinit

  real :: delt_adj, dt0
!  real :: delt_cushion = 1.5
  real :: delt_cushion
  real :: delt_minimum 
  real, save :: time_reinit(2)=0.
  logical :: abort_rapid_time_step_change
  logical :: first=.true.
  logical :: in_memory
contains
  subroutine wnml_gs2_reinit(unit)
    implicit none
    integer :: unit
          write (unit, *)
          write (unit, fmt="(' &',a)") "reinit_knobs"
          write (unit, fmt="(' delt_adj = ',e16.10)") delt_adj
          write (unit, fmt="(' delt_minimum = ',e16.10)") delt_minimum
          write (unit, fmt="(' /')")       
  end subroutine wnml_gs2_reinit

  subroutine reset_time_step (istep, exit, job_id)

    use collisions, only: c_reset => reset_init, vnmult
    use dist_fn, only: d_reset => reset_init
    use fields, only: f_reset => reset_init, init_fields
    use init_g, only: g_reset => reset_init
    use run_parameters, only: fphi, fapar, fbpar, reset
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt
    use gs2_save, only: gs2_save_for_restart
    use dist_fn_arrays, only: gnew, g_restart_tmp
    use gs2_time, only: user_time, code_dt_min
    use nonlinear_terms, only: nl_reset => reset_init
    use mp, only: proc0
    use file_utils, only: error_unit
    use antenna, only: dump_ant_amp
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use job_manage, only: time_message
    logical :: exit, reset_in
    integer :: istep 
    integer, save :: istep_last = -1 ! allow adjustment on first time step
    integer :: istatus
    integer, save :: nconsec=0
    integer, intent (in), optional :: job_id
    integer :: iostat
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
       exit = .true.
       if (proc0) write(error_unit(), *) 'Time step changing rapidly.  Abort run.'
       return
    end if

    if (code_dt/delt_adj <= code_dt_min) then
       code_dt = code_dt_min  ! set it so restart is ok
       exit = .true.
       if (proc0) write(error_unit(), *) 'Time step wants to fall below delt_min.  Abort run.'
       return
    end if

    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    !First disable the reset flag so we can call 
    !routines needed in reinit
    reset_in=reset
    reset=.false.

    !If we want to do restarts in memory then try to allocate storage
    if(in_memory)then
       !Try to allocate storage to hold g
       allocate(g_restart_tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc),stat=iostat)

       !Should really make tmp copy of fields as well.

       !If allocate failed
       if(iostat.ne.0)then
          !Disable in_memory flag
          in_memory=.false.
          !Print error message
          if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for g --> Reverting to file based restart"
       else
          g_restart_tmp=gnew
       endif
    endif

    if(.not.in_memory)then

       !Should really do this with in_memory=.true. as well but
       !not sure that we really need to as we never read in the dumped data.
       if (proc0) call dump_ant_amp

       call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    endif

    gnew = 0.

! change timestep 

! If timestep is too big, make it smaller
    if (code_dt > code_dt_cfl) then
       code_dt = code_dt/delt_adj

! If timestep is too small, make it bigger
    else if (code_dt < min(dt0, code_dt_cfl/delt_adj/delt_cushion)) then
       code_dt = min(code_dt*delt_adj, dt0)

    endif
    
    call save_dt (code_dt)

    if (proc0 .and. .not. present(job_id)) write(*,*) 'Changing time step to ', user_dt
    
! prepare to reinitialize inversion matrix, etc.
    call d_reset
    call c_reset
    call f_reset
    call g_reset(.not.in_memory)
    call nl_reset

! reinitialize
    call init_fields

    if (proc0 .and. .not. present(job_id)) call time_message(.true.,time_reinit,' Re-initialize')

    istep_last = istep

    !Now re-enable reset so we leave it in the same state as on entering
    reset=reset_in

    !Deallocate tmp memory
    if(allocated(g_restart_tmp)) deallocate(g_restart_tmp)
  end subroutine reset_time_step

!  subroutine check_time_step (istep, reset, exit)
  subroutine check_time_step (reset, exit)

    use gs2_time, only: code_dt_cfl, code_dt

!    integer :: istep
    logical :: reset, exit

    if (first) call init_reinit
    first = .false.
    reset = .false.

! nothing to do if exiting in this iteration
    if (exit) return

! If timestep is too big, make it smaller
    if (code_dt > code_dt_cfl) reset = .true. !Note this logic is repeated in gs2_time::check_time_step_too_large
       
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

  end subroutine init_reinit

end module gs2_reinit

