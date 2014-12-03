!> A module which allows multiple values of certain parameters to be considered within a single simulation.
!!
!! In its simplest form, it starts with a given value of the parameter in question, runs for a given number of timesteps, and then changes the parameter by a given increment, until a lower limit is reached.
!!
!! In a more advanced scenario, the parameter scan continues until a different condition is satisfied: for example, zero heat flux.
!!
!! In a third scenario, the parameter is varied using a root finding algorithm.
!!
!! In addition, the condition for changing the parameter may be changed from a simple number of time steps to, for example, reaching a saturated state.

module parameter_scan

  !These should really be moved to specific routines where required
  use parameter_scan_arrays, only: current_scan_parameter_value, write_scan_parameter, run_scan
  use parameter_scan_arrays, only: scan_parameter_switch, scan_parameter_tprim, scan_parameter_g_exb
  use parameter_scan_arrays, only: scan_spec, hflux_tot, momflux_tot, phi2_tot, nout

  implicit none

  private

  public :: init_parameter_scan
  public :: finish_parameter_scan
  public :: update_scan_parameter_value
  public :: allocate_target_arrays
  public :: target_parameter_switch, scan_type_switch, scan_type_none
  public :: target_parameter_hflux_tot, target_parameter_momflux_tot, target_parameter_phi2_tot
  public :: scan_restarted

  logical :: scan_restarted

  integer :: target_parameter_switch, &
             scan_type_switch, &
             increment_condition_switch


  integer, parameter :: scan_type_none = 1, &
                        scan_type_range = 2, &
                        scan_type_target =3, &
                        scan_type_root_finding = 4

  integer, parameter :: target_parameter_hflux_tot = 1, &
                        target_parameter_momflux_tot = 2, &
                        target_parameter_phi2_tot = 3 

  integer, parameter :: increment_condition_n_timesteps = 1, &
                        increment_condition_delta_t = 2, &
                        increment_condition_saturated = 3

  real :: par_start, par_end, par_inc
  integer :: nstep_init, nstep_inc
  real :: delta_t_init, delta_t_inc
  real :: target_val
  integer :: scan_output_file
  real :: current_target_value = 0.0

contains
  subroutine init_parameter_scan
    use gs2_save, only: restore_current_scan_parameter_value 
    use init_g, only: init_init_g
    use file_utils, only: open_output_file
    logical, save :: initialized = .false.
    !write (*,*) "initializing parameter_scan"

    if (initialized) return
    initialized = .true.


    call read_parameters

    select case (scan_type_switch) 
    case (scan_type_none)
      write_scan_parameter = .false.
      run_scan = .false.
      return
    case (scan_type_target)
      !call allocate_target_arrays
    case (scan_type_root_finding)
      !call allocate_target_arrays
    end select


    call open_output_file(scan_output_file, ".par_scan")
      
    write_scan_parameter = .true.
    run_scan = .true.

    ! scan_restarted is set manually 
    if (scan_restarted) then
      call init_init_g ! To get restart file name
      call restore_current_scan_parameter_value(current_scan_parameter_value)
    end if
   
    !write (*,*) "initialized parameter_scan"

  end subroutine init_parameter_scan

  !subroutine set_restarted_scan_parameter_value(value)
    !restarted = .true.
    !current_scan_parameter_value = value
  !end subroutine set_restarted_scan_parameter_value 

  subroutine allocate_target_arrays(nwrite,write_nl_flux)
!    use gs2_diagnostics, only: nwrite, write_nl_flux
!    use gs2_diagnostics, only: write_flux_line
    use run_parameters, only: nstep
    use mp, only : mp_abort
    implicit none
    integer, intent(in) :: nwrite
    logical, intent(in) :: write_nl_flux
    !if (.not. write_flux_line) &
     !call mp_abort("write_flux_line must be set to true for target mode")

    !select case (scan_type_switch) 
    !case (scan_type_none)
      !return
    !case (scan_type_range)
      !return
    !end select

     allocate(hflux_tot(nstep/nwrite+1))
     allocate(momflux_tot(nstep/nwrite+1))
     allocate(phi2_tot(nstep/nwrite+1))
    if (scan_type_switch == scan_type_none) return

    select case (target_parameter_switch)
    case (target_parameter_hflux_tot)
      if (.not. write_nl_flux) &
       call mp_abort("write_nl_flux must be set to true for hflux target mode")
     case (target_parameter_momflux_tot)
      if (.not. write_nl_flux) &
       call mp_abort("write_nl_flux must be set to true for momflux target mode")
     end select 
     !write (*,*) "allocating target arrays ", nwrite, ",", nstep 
  end subroutine allocate_target_arrays


  subroutine finish_parameter_scan
    use file_utils, only: close_output_file
    if (allocated(hflux_tot)) deallocate(hflux_tot)
    if (allocated(momflux_tot)) deallocate(momflux_tot)
    if (allocated(phi2_tot)) deallocate(phi2_tot)
    call close_output_file(scan_output_file)
  end subroutine finish_parameter_scan
  
  subroutine update_scan_parameter_value(istep, reset, exit)
    use mp, only : mp_abort
    use gs2_time, only: user_time
    logical, intent (inout) :: exit, reset
    integer, intent (in) :: istep
    real :: increment
    logical :: increment_condition_satisfied, target_reached

    !write (*,*) "update_scan_parameter_value"
    increment_condition_satisfied = .false.
    target_reached = .false.

    select case (scan_type_switch)
    case (scan_type_none)
      return
    case (scan_type_range)
      write (scan_output_file, fmt="('time=  ',e11.4,'  parameter=  ',e11.4)") &
          user_time, current_scan_parameter_value
      call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      if (.not. increment_condition_satisfied) return
      increment = par_inc
      if ((increment .gt. 0.0 &
         .and. current_scan_parameter_value + increment .gt. par_end) &
         .or. (increment .lt. 0.0 &
         .and. current_scan_parameter_value + increment .lt. par_end)) &
         then
         exit = .true.
       else  
          call increment_scan_parameter(increment, reset)
       end if
     case (scan_type_target)
      write (scan_output_file, &
        fmt="('time=  ',e11.4,'  parameter=  ',e11.4,'  target=  ',e11.4)") &
          user_time, current_scan_parameter_value, current_target_value 
      call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      if (.not. increment_condition_satisfied) return
      increment = par_inc
      call check_target_reached(target_reached)
      if (target_reached) then
        exit = .true.
      else
        call increment_scan_parameter(increment, reset)
      end if
    case (scan_type_root_finding)
      call mp_abort("scan_type_root_finding not implemented yet!")
      !write (scan_output_file, &
        !fmt="(time=  e11.4  parameter=  e11.4  target=  e11.4)") &
          !user_time, current_scan_parameter_value, current_target_value 
      !call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      !if (.not. increment_condition_satisfied) return
      !increment = par_inc
      !call check_target_reached(target_reached)
      !if (target_reached) then
        !exit = .true.
      !else
        !call get_root_finding_increment(increment)
        !call increment_scan_parameter(increment, reset)
      !end if
    end select
    if (exit) write (scan_output_file, *) "Parameter scan is complete...."
  end subroutine update_scan_parameter_value

  subroutine check_target_reached(target_reached)
    logical, intent (out) :: target_reached
    logical, save :: first = .true.
    real, save :: last_value = 0.0
   
    target_reached = .false.
    
    last_value = current_target_value
    select case (target_parameter_switch)
    case (target_parameter_hflux_tot)
      current_target_value = hflux_tot(nout)
    case (target_parameter_momflux_tot)
      current_target_value = momflux_tot(nout)
    case (target_parameter_phi2_tot)
      current_target_value = phi2_tot(nout)
    end select

    if (first) then
      last_value = current_target_value
      first = .false.
      return
    end if

    if ((last_value .lt. target_val .and. &
          current_target_value .gt. target_val) .or. &
          (last_value .gt. target_val .and. &
          current_target_value .lt. target_val)) &
             target_reached = .true.

  end subroutine check_target_reached

  subroutine increment_scan_parameter(increment, reset)
    real, intent (in) :: increment
    logical, intent (inout) :: reset
    current_scan_parameter_value = current_scan_parameter_value + increment
    call set_scan_parameter(reset)
  end subroutine increment_scan_parameter
     
  subroutine set_scan_parameter(reset)
    !use parameter_scan_arrays, only: current_scan_parameter_value
    !use parameter_scan_arrays, only: scan_parameter_switch
    !use parameter_scan_arrays, only: scan_parameter_tprim
    !use parameter_scan_arrays, only: scan_parameter_g_exb
    !use parameter_scan_arrays
    use species, only: spec 
    use dist_fn, only: g_exb
    use mp, only: proc0
    logical, intent (inout) :: reset
     
    select case (scan_parameter_switch)
    case (scan_parameter_tprim)
       spec(scan_spec)%tprim = current_scan_parameter_value
       if (proc0) write (*,*) &
         "Set scan parameter tprim_1 to ", spec(scan_spec)%tprim
       reset = .true.
    case (scan_parameter_g_exb)
       g_exb = current_scan_parameter_value
       if (proc0) write (*,*) &
         "Set scan parameter g_exb to ", g_exb
       reset = .false.
    end select
  end subroutine set_scan_parameter

  subroutine check_increment_condition_satisfied(istep, increment_condition_satisfied)
    use gs2_time, only: user_time
    use mp, only : mp_abort, proc0
    logical, intent (out) :: increment_condition_satisfied
    integer, intent (in) :: istep
    integer, save :: last_timestep = 0
    real, save :: last_time = 0.0

    select case (increment_condition_switch)
    case (increment_condition_n_timesteps)
      if (istep .lt. nstep_init) then
        last_timestep = istep
        increment_condition_satisfied = .false.
        return
      end if 
      if ((istep - last_timestep) .gt. nstep_inc) then
        last_timestep = istep
        increment_condition_satisfied = .true.
        if (proc0) write (*,*) 'Updating parameter at ', user_time 
        return
      end if
    case (increment_condition_delta_t)
      if (user_time .lt. delta_t_init) then
        last_time = user_time
        increment_condition_satisfied = .false.
        return
      end if 
      if ((user_time - last_time) .gt. delta_t_inc) then
        last_time = user_time
        increment_condition_satisfied = .true.
        if (proc0) write (*,*) 'Updating parameter at ', user_time 
        return
      end if
    case (increment_condition_saturated)
      call mp_abort("increment_condition_saturated not implemented yet!")
      !if (istep .lt. nstep_init) then
        !last_timestep = istep
        !increment_condition_satisfied = .false.
        !return
      !end if 
      !call get_target_value_array(target_value_array)
      !call check_saturated(target_value_array, &
                           !istep, &
                           !last_timestep, &
                           !increment_condition_satisfied)
      !if (increment_condition_satisfied) last_timestep = istep
    end select 
  end subroutine check_increment_condition_satisfied

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (2), parameter :: scan_parameter_opts = &
         (/ text_option('tprim', scan_parameter_tprim), &
            text_option('g_exb', scan_parameter_g_exb) /)
    character(20) :: scan_par
    type (text_option), dimension (4), parameter :: scan_type_opts = &
         (/ text_option('none', scan_type_none), &
            text_option('range', scan_type_range), &
            text_option('target', scan_type_target), &
            text_option('root_finding', scan_type_root_finding) /)
    character(20) :: scan_type
    type (text_option), dimension (3), parameter :: target_parameter_opts = &
         (/ text_option('hflux_tot', target_parameter_hflux_tot), &
            text_option('momflux_tot', target_parameter_momflux_tot), &
            text_option('phi2_tot', target_parameter_phi2_tot) /)
    character(20) :: target_par
    type (text_option), dimension (3), parameter :: increment_condition_opts = &
         (/ text_option('n_timesteps', increment_condition_n_timesteps), &
            text_option('delta_t', increment_condition_delta_t), &
            text_option('saturated', increment_condition_saturated) /)
    character(20) :: inc_con 
    namelist /parameter_scan_knobs/ &
            scan_type, &
            scan_par, &
            target_par,& 
            par_start, &
            par_end, &
            par_inc, &
            inc_con, &
            nstep_init, &
            nstep_inc, &
            delta_t_init, &
            delta_t_inc, &
            scan_spec, &
            scan_restarted, &
            target_val

    integer :: ierr, in_file
    logical :: exist

    if (proc0) then
       scan_par = 'tprim'
       scan_type = 'none'
       target_par = 'hflux_tot'
       inc_con = 'delta_t'
       par_start = 0.0
       par_end = 0.0
       par_inc = 0.0
       target_val = 0.0
       nstep_init = 0
       nstep_inc = 0
       delta_t_init = 0.0
       delta_t_inc = 0.0
       scan_spec = 1
       scan_restarted = .false.

       in_file = input_unit_exist ("parameter_scan_knobs", exist)
       if (exist) read (unit=in_file, nml=parameter_scan_knobs)

       ierr = error_unit()
       call get_option_value &
            (scan_par, scan_parameter_opts, scan_parameter_switch, &
            ierr, "scan_par in parameter_scan_knobs",.true.)
       call get_option_value &
            (scan_type, scan_type_opts, scan_type_switch, &
            ierr, "scan_type in parameter_scan_knobs",.true.)
       call get_option_value &
            (target_par, target_parameter_opts, target_parameter_switch, &
            ierr, "target_par in parameter_scan_knobs",.true.)
       call get_option_value &
            (inc_con, increment_condition_opts, &
            increment_condition_switch, &
            ierr, "inc_con in parameter_scan_knobs",.true.)

    end if

    call broadcast (scan_parameter_switch)
    call broadcast (scan_type_switch)
    call broadcast (target_parameter_switch)
    call broadcast (par_start)
    call broadcast (par_end)
    call broadcast (par_inc)
    call broadcast (increment_condition_switch)
    call broadcast (nstep_init)
    call broadcast (nstep_inc)
    call broadcast (delta_t_init)
    call broadcast (delta_t_inc)
    call broadcast (scan_spec)
    call broadcast (target_val)
    call broadcast (scan_restarted)

    if (.not. scan_restarted) current_scan_parameter_value = par_start

  end subroutine read_parameters
end module parameter_scan

