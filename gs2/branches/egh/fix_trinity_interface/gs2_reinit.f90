module gs2_reinit
  implicit none

  public :: reset_time_step, delt_adj
  public :: check_time_step, time_reinit
  public :: init_reinit, wnml_gs2_reinit
  public :: reduce_time_step, increase_time_step

  !> Save the current state of the fields and 
  !! distribution function, either to a file 
  !! or to a temporary array, depending on the
  !! value of in_memory. (NB if there is 
  !! sufficient memory to allocate a temporary
  !! array, in_memory will be overriden to 
  !! false).
  public :: save_fields_and_dist_fn


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

    complex, dimension(:,:,:), allocatable :: phi_tmp, apar_tmp, bpar_tmp
contains

  subroutine gs2_reinit_unit_test_set_in_memory(in_memory_in)
    logical, intent(in) :: in_memory_in
    in_memory = in_memory_in
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

  subroutine save_fields_and_dist_fn
    use dist_fn_arrays, only: gnew, g_restart_tmp
    use gs2_save, only: gs2_save_for_restart
    use mp, only: proc0
    use collisions, only: vnmult
    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use file_utils, only: error_unit
    use kt_grids, only: ntheta0, naky
    use run_parameters, only: fphi, fapar, fbpar
    use fields, only:  force_maxwell_reinit
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_time, only: user_time, user_dt
    use antenna, only: dump_ant_amp
    integer :: iostat, istatus
    !If we want to do restarts in memory then try to allocate storage
    if(in_memory)then
       !Try to allocate storage to hold g
       allocate(g_restart_tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc),stat=iostat)

       !If allocate failed
       if(iostat.ne.0)then
          !Disable in_memory flag
          in_memory=.false.
          !Print error message
          if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for g --> Reverting to file based restart"
       else
          !Copy into temporary
          g_restart_tmp=gnew
       endif

       !!!!
       !! NOW WE MAKE COPIES OF THE FIELDS
       !! --> Don't bother if force_maxwell_reinit as we're going to recalculate
       !!!!

       !Try to allocate storage to hold phi
       if(fphi.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
          allocate(phi_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

          !If allocate failed
          if(iostat.ne.0)then
             !Disable in_memory flag
             in_memory=.false.
             !Print error message
             if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for phi --> Reverting to file based restart"
          else
             !Copy into temporary
             phi_tmp=phinew
          endif
       endif

       !Try to allocate storage to hold apar
       if(fapar.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
          allocate(apar_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

          !If allocate failed
          if(iostat.ne.0)then
             !Disable in_memory flag
             in_memory=.false.
             !Print error message
             if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for apar --> Reverting to file based restart"
          else
             !Copy into temporary
             apar_tmp=aparnew
          endif
       endif

       !Try to allocate storage to hold bpar
       if(fbpar.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
          allocate(bpar_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

          !If allocate failed
          if(iostat.ne.0)then
             !Disable in_memory flag
             in_memory=.false.
             !Print error message
             if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for bpar --> Reverting to file based restart"
          else
             !Copy into temporary
             bpar_tmp=bparnew
          endif
       endif

    endif

    if(.not.in_memory)then
       !Should really do this with in_memory=.true. as well but
       !not sure that we really need to as we never read in the dumped data.
       if (proc0) call dump_ant_amp

       call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    endif

  end subroutine save_fields_and_dist_fn

  subroutine reinit_gk_and_field_equations
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: g_restart_tmp
    use fields_arrays, only: phinew, aparnew, bparnew
    use collisions, only: c_reset => reset_init
    use dist_fn, only: d_reset => reset_init
    use fields, only: f_reset => reset_init, init_fields, force_maxwell_reinit
    use init_g, only: g_reset => reset_init
    use nonlinear_terms, only: nl_reset => reset_init
! prepare to reinitialize inversion matrix, etc.
    call d_reset
    call c_reset
    call f_reset
    call g_reset(.not.in_memory)
    call nl_reset

! reinitialize
    call init_fields

!Update fields if done in memory
!Don't need/want to update if force_maxwell_reinit
    if(in_memory.and.(.not.force_maxwell_reinit))then
       if(fphi.gt.0) phinew=phi_tmp
       if(fapar.gt.0) aparnew=apar_tmp
       if(fbpar.gt.0) bparnew=bpar_tmp
    endif
    !Deallocate tmp memory
    if(allocated(g_restart_tmp)) deallocate(g_restart_tmp)
    if(allocated(phi_tmp)) deallocate(phi_tmp)
    if(allocated(apar_tmp)) deallocate(apar_tmp)
    if(allocated(bpar_tmp)) deallocate(bpar_tmp)
  end subroutine reinit_gk_and_field_equations


  subroutine reset_time_step (istep, my_exit, job_id)
    use run_parameters, only: reset
    use gs2_time, only: code_dt, user_dt, code_dt_cfl, save_dt
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt_min
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

    call reinit_gk_and_field_equations
    
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
  end subroutine init_reinit
end module gs2_reinit

