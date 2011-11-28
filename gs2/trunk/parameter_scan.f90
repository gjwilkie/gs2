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

  !use parameter_scan_arrays, only: current_scan_parameter_value
  use parameter_scan_arrays

  public :: init_parameter_scan
  public :: update_scan_parameter_value
  public :: scan_restarted




  private

  integer :: scan_parameter_switch, &
             target_parameter_switch, &
             scan_type_switch, &
             increment_condition_switch

  integer, parameter :: scan_parameter_tprim = 1, &
                        scan_parameter_g_exb = 2

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

  real :: parameter_start_value, parameter_end_value, parameter_increment
  integer :: n_timesteps_initial, n_timesteps_increment
  integer :: species_scan_index
  real :: delta_t_initial, delta_t_increment
  real :: target_value
  logical :: scan_restarted


contains
  subroutine init_parameter_scan
    use gs2_save, only: restore_current_scan_parameter_value 
    logical, save :: initialized = .false.
    write (*,*) "initializing parameter_scan"

    if (initialized) return
    initialized = .true.


    call read_parameters

    select case (scan_type_switch) 
    case (scan_type_none)
      write_scan_parameter = .false.
      return
    case (scan_type_target)
      call allocate_target_arrays
    case (scan_type_root_finding)
      call allocate_target_arrays
    end select
      
    write_scan_parameter = .true.

    ! scan_restarted is set manually 
    if (scan_restarted) call &
      restore_current_scan_parameter_value(current_scan_parameter_value)
   
    ! To set the initial value of the actual parameter being varied:
    call increment_scan_parameter(0.0)

    write (*,*) "initialized parameter_scan"

  end subroutine init_parameter_scan

  !subroutine set_restarted_scan_parameter_value(value)
    !restarted = .true.
    !current_scan_parameter_value = value
  !end subroutine set_restarted_scan_parameter_value 

  subroutine allocate_target_arrays
    use gs2_diagnostics, only: nwrite, write_nl_flux, write_flux_line
    use run_parameters, only: nstep
    use mp, only : mp_abort

    !if (.not. write_flux_line) &
     !call mp_abort("write_flux_line must be set to true for target mode")
    select case (target_parameter_switch)
    case (target_parameter_hflux_tot)
      if (.not. write_nl_flux) &
       call mp_abort("write_nl_flux must be set to true for hflux target mode")
     case (target_parameter_momflux_tot)
      if (.not. write_nl_flux) &
       call mp_abort("write_nl_flux must be set to true for momflux target mode")
     end select 
     write (*,*) "allocating target arrays ", nwrite, ",", nstep 
     allocate(hflux_tot(nstep/nwrite))
     allocate(momflux_tot(nstep/nwrite))
     allocate(phi2_tot(nstep/nwrite))
  end subroutine allocate_target_arrays


  subroutine finish_parameter_scan
    if (allocated(hflux_tot)) deallocate(hflux_tot)
    if (allocated(momflux_tot)) deallocate(momflux_tot)
    if (allocated(phi2_tot)) deallocate(phi2_tot)
  end subroutine finish_parameter_scan
  
  subroutine update_scan_parameter_value(istep, reset, exit)
    use mp, only : mp_abort
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
      call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      if (.not. increment_condition_satisfied) return
      increment = parameter_increment
      if ((increment .gt. 0.0 &
         .and. current_scan_parameter_value + increment .gt. parameter_end_value) &
         .or. (increment .lt. 0.0 &
         .and. current_scan_parameter_value + increment .lt. parameter_end_value)) &
         then
         exit = .true.
         return
       end if
       call increment_scan_parameter(increment)
       reset = .true.
     case (scan_type_target)
      !call mp_abort("scan_type_target not implemented yet!")
      call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      if (.not. increment_condition_satisfied) return
      increment = parameter_increment
      call check_target_reached(target_reached)
      if (target_reached) then
        exit = .true.
        return
      end if
      call increment_scan_parameter(increment)
      reset = .true.
    case (scan_type_root_finding)
      call mp_abort("scan_type_root_finding not implemented yet!")
      !call check_increment_condition_satisfied(istep, increment_condition_satisfied)
      !if (.not. increment_condition_satisfied) return
      !increment = parameter_increment
      !call check_target_reached(target_reached)
      !if (target_reached) then
        !exit = .true.
        !return
      !end if
      !call get_root_finding_increment(increment)
      !call increment_scan_parameter(increment)
      !reset = .true.
    end select
  end subroutine update_scan_parameter_value

  subroutine check_target_reached(target_reached)
    logical, intent (out) :: target_reached
    logical, save :: first = .true.
    real, save :: last_value = 0.0
    real, save :: current_value = 0.0
   
    target_reached = .false.
    
    last_value = current_value
    select case (target_parameter_switch)
    case (target_parameter_hflux_tot)
      current_value = hflux_tot(nout)
    case (target_parameter_momflux_tot)
      current_value = momflux_tot(nout)
    case (target_parameter_phi2_tot)
      current_value = phi2_tot(nout)
    end select

    if (first) then
      last_value = current_value
      first = .false.
      return
    end if

    if ((last_value .lt. target_value .and. &
          current_value .gt. target_value) .or. &
          (last_value .gt. target_value .and. &
          current_value .lt. target_value)) &
             target_reached = .true.

  end subroutine check_target_reached

  subroutine increment_scan_parameter(increment)
    use species, only: spec 
    use dist_fn, only: g_exb
    real, intent (in) :: increment
     
    current_scan_parameter_value = current_scan_parameter_value + increment
    select case (scan_parameter_switch)
    case (scan_parameter_tprim)
       spec(species_scan_index)%tprim = current_scan_parameter_value
       write (*,*) "tprim_1 is ", spec(species_scan_index)%tprim
    case (scan_parameter_g_exb)
       g_exb = current_scan_parameter_value
    end select
  end subroutine increment_scan_parameter
     



  subroutine check_increment_condition_satisfied(istep, increment_condition_satisfied)
    use gs2_time, only: user_time
    use mp, only : mp_abort
    logical, intent (out) :: increment_condition_satisfied
    integer, intent (in) :: istep
    logical, save :: first = .true.
    integer, save :: last_timestep = 0
    real, save :: last_time = 0.0

    select case (increment_condition_switch)
    case (increment_condition_n_timesteps)
      if (istep .lt. n_timesteps_initial) then
        last_timestep = istep
        increment_condition_satisfied = .false.
        return
      end if 
      if ((istep - last_timestep) .gt. n_timesteps_increment) then
        last_timestep = istep
        increment_condition_satisfied = .true.
        return
      end if
    case (increment_condition_delta_t)
      if (user_time .lt. delta_t_initial) then
        last_time = user_time
        increment_condition_satisfied = .false.
        return
      end if 
      if ((user_time - last_time) .gt. delta_t_increment) then
        write (*,*) 'Updating parameter at ', user_time 
        last_time = user_time
        increment_condition_satisfied = .true.
        return
      end if
    case (increment_condition_saturated)
      call mp_abort("increment_condition_saturated not implemented yet!")
      !if (istep .lt. n_timesteps_initial) then
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
    character(20) :: scan_parameter
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
    character(20) :: target_parameter
    type (text_option), dimension (3), parameter :: increment_condition_opts = &
         (/ text_option('n_timesteps', increment_condition_n_timesteps), &
            text_option('delta_t', increment_condition_delta_t), &
            text_option('saturated', increment_condition_saturated) /)
    character(20) :: increment_condition 
    namelist /parameter_scan_knobs/ &
            scan_type, &
            scan_parameter, &
            target_parameter,& 
            parameter_start_value, &
            parameter_end_value, &
            parameter_increment, &
            increment_condition, &
            n_timesteps_initial, &
            n_timesteps_increment, &
            delta_t_initial, &
            delta_t_increment, &
            species_scan_index, &
            scan_restarted, &
            target_value

    integer :: ierr, in_file
    logical :: exist

    if (proc0) then
       scan_parameter = 'tprim'
       scan_type = 'none'
       target_parameter = 'hflux_tot'
       parameter_start_value = 0.0
       parameter_end_value = 0.0
       parameter_increment = 0.0
       target_value = 0.0
       n_timesteps_initial = 0
       n_timesteps_increment = 0
       delta_t_initial = 0.0
       delta_t_increment = 0.0
       species_scan_index = 1
       scan_restarted = .false.

       in_file = input_unit_exist ("parameter_scan_knobs", exist)
       if (exist) read (unit=in_file, nml=parameter_scan_knobs)

       ierr = error_unit()
       call get_option_value &
            (scan_parameter, scan_parameter_opts, scan_parameter_switch, &
            ierr, "scan_parameter in parameter_scan_knobs")
       call get_option_value &
            (scan_type, scan_type_opts, scan_type_switch, &
            ierr, "scan_type in parameter_scan_knobs")
       call get_option_value &
            (target_parameter, target_parameter_opts, target_parameter_switch, &
            ierr, "target_parameter in parameter_scan_knobs")
       call get_option_value &
            (increment_condition, increment_condition_opts, &
            increment_condition_switch, &
            ierr, "increment_condition in parameter_scan_knobs")

    end if

    call broadcast (scan_parameter_switch)
    call broadcast (scan_type_switch)
    call broadcast (target_parameter_switch)
    call broadcast (parameter_start_value)
    call broadcast (parameter_end_value)
    call broadcast (parameter_increment)
    call broadcast (increment_condition_switch)
    call broadcast (n_timesteps_initial)
    call broadcast (n_timesteps_increment)
    call broadcast (delta_t_initial)
    call broadcast (delta_t_increment)
    call broadcast (species_scan_index)
    call broadcast (target_value)

    if (.not. scan_restarted) current_scan_parameter_value = parameter_start_value

  end subroutine read_parameters
end module parameter_scan

