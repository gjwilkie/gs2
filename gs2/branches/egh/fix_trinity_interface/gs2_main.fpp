# ifndef MAKE_LIB

!> This module contains the functions gs2_main::run_gs2, gs2_main::finish_gs2 and gs2_main::reset_gs2, whose
!! purpose should be reasonably obvious. These functions were originally part of
!! program GS2, but have now been moved to this module so that GS2 itself can be
!! called as a library by, for example, Trinity. All the program GS2 does is
!! include this module and call run_gs2.

module gs2_main
  use gs2_init, only: init_level_type, init_level_list
  implicit none
  public :: run_gs2, finish_gs2, reset_gs2, trin_finish_gs2

  public :: gs2_program_state_type
  public :: initialize_gs2, initialize_equations
  public :: initialize_diagnostics, evolve_equations, run_eigensolver
  public :: finalize_diagnostics, finalize_equations, finalize_gs2
  public :: calculate_outputs
  public :: override

  public :: old_iface_state
  !> Unit tests

  !> This function calls reset_gs2 using
  !! all the default parameters, multiplied
  !! by a factor. It is used
  !! in the gs2_reinit unit test
  public :: gs2_main_unit_test_reset_gs2

  !> Starts the global wall clock timer
  !! used by check time. This is useful
  !! if you have multiple copies of gs2 
  !! running but you don't want to start 
  !! them all at the same time, but you
  !! want them all to respect avail_cpu_time
  public :: initialize_wall_clock_timer


  type gs2_timers_type
    real :: init(2) 
    real :: advance(2) = 0.
    real :: finish(2) = 0.
    real :: total(2) = 0. 
    real :: diagnostics(2)=0.
    !real :: interval
    real :: main_loop(2)
  end type gs2_timers_type

  type gs2_outputs_type
    !> The gradient of the flux tube volume wrt
    !! the flux label: related to the surface
    !! area of the flux tube
    real :: dvdrho
    !> The average gradient of the flux tube
    !! label <|grad rho|>
    real :: grho
    !> Particle flux by species
    real, dimension (:), pointer :: pflux
    !> Heat flux by species
    real, dimension (:), pointer :: qflux
    !> Turbulent heating by species
    real, dimension (:), pointer :: heat
    !> Momentum flux
    real :: vflux
  end type gs2_outputs_type

  type gs2_program_state_type

    ! Flags indicating the current state of the 
    ! program (used for error checking)
    logical :: gs2_initialized = .false.
    logical :: equations_initialized = .false.
    logical :: diagnostics_initialized = .false.

    !> A type for keeping track of the current
    !! initialization level of gs2
    type(init_level_type) :: init



    ! Timers
    type(gs2_timers_type) :: timers

    !> The exit flag is set to true by any 
    !! part of the main timestep loop that 
    !! wants to cause the loop to exit
    logical :: exit = .false.

    !> Whether the run has converged to a
    !! stationary state
    logical :: converged = .false.

    !> Whether to run the eigenvalue solver
    !! or not. Set equal to the input value
    !! in intialize_equations. Typically
    !! only important for the standalone
    !! gs2 program
    logical :: do_eigsolve = .false.

    !> Number of timesteps in evolve_equations
    !! It is set equal to run_parameters::nstep in 
    !! initialize_equations and is the maximum
    !! number of timesteps that can be run 
    !! without reinitalising the diagnostics. 
    !! Do not modify!
    integer :: nstep

    !> Gets set to the final value of istep
    !! in evolve equations. Any future calls
    !! to evolve_equations will increment this
    !! further. A call to finalize_diagnostics
    !! will set it back to 1. Note that setting
    !! this manually is NOT advised. 
    integer :: istep_end = 1

    !> Whether to print out debug messages
    !logical :: debug = .false.
    integer :: verb = 3

    !> Parameters pertaining to cases when gs2
    !! is being used as library.
    !! external_job_id is not to confused with the parameter
    !! job in mp, which identifies the subjob if
    !! running in list mode or with nensembles > 1
    !integer :: trin_job = -1
    integer :: external_job_id = -1

    !> is_external_job should be set to true when GS2
    !! is being used as a library 
    logical :: is_external_job = .false.

    !> Set true if using trinity. This does several things: 
    !!  * it forces the calculation of the fluxes
    !!  * it causes the species and theta_grid name lists to use parameters
    !!     provided by trinity
    !! Setting this flag true automatically sets is_external_job 
    !! to true
    logical :: is_trinity_job = .false.


    !> Parameters to be used when passing in an external communicator
    logical :: mp_comm_external = .false.
    integer :: mp_comm

    logical :: run_name_external = .false.
    character(2000) :: run_name

    !> Whether this is a list mode run
    logical :: list
    !> The number of identical runs happening
    !! simultaneously (used for ensemble averaging).
    !! Cannot be used in conjunction with list mode
    integer :: nensembles = 1

    ! Outputs (e.g. for Trinity)
    type(gs2_outputs_type) :: outputs

   
  end type gs2_program_state_type

  private

  interface override
    module procedure override_parameter_spec
    module procedure override_parameter_nospec
  end interface override

  !> This object is used for implementing the old interface
  !! and should not be modified
  type(gs2_program_state_type) :: old_iface_state
 
contains
# endif

  subroutine initialize_wall_clock_timer
    use job_manage, only: init_checktime
    call init_checktime
  end subroutine initialize_wall_clock_timer

  subroutine initialize_gs2(state)
    use file_utils, only: init_file_utils
    use file_utils, only: run_name, run_name_target
    use gs2_init, only: init_gs2_init
    use job_manage, only: checktime, time_message
    use job_manage, only: init_checktime, checktime_initialized
    use job_manage, only: job_fork
    use mp, only: init_mp, broadcast
    use mp, only: iproc, nproc, proc0
    use mp, only: trin_flag, job
    use redistribute, only: using_measure_scatter
    use run_parameters, only: avail_cpu_time
    use runtime_tests, only: verbosity
    use unit_tests, only: debug_message, set_job_id
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (state%init%level .ge. init_level_list%gs2) then
      write  (*,*) "ERROR: Called initialize_gs2 twice &
        without calling finalize_gs2"
      stop 1
    end if



    call debug_message(4, 'gs2_main::initialize_gs2 starting initialization')

    if (state%mp_comm_external) then
       call init_mp (state%mp_comm)
    else
       call init_mp
    end if

    if (state%is_trinity_job) state%is_external_job = .true.
    if (state%is_external_job) then 
      call broadcast(state%external_job_id)
      call set_job_id(state%external_job_id)
    end if

    ! We now manually set trin_flag, since we may now be passing 
    ! in a communicator but not running trinity.
    trin_flag = state%is_trinity_job
    call debug_message(state%verb, 'gs2_main::initialize_gs2 initialized mp')
    ! I don't think these should be necessary
    !call broadcast(state%trin_job)
    !call set_job_id(trin_job)

    call reset_timers(state%timers)
    if (.not. checktime_initialized) call init_checktime ! <doc> Initialize timer </doc>
    call debug_message(state%verb, &
      'gs2_main::initialize_gs2 called init_checktime')
  
    if (proc0) then
       ! Report # of processors being used </doc>
       if (nproc == 1) then
         write(*,*) 'Running on ',nproc,' processor'
       else
         if(state%mp_comm_external) then
            write(*,*) 'Job ',state%external_job_id,'Running on ',nproc,' processors'
         else
            write(*,*) 'Running on ',nproc,' processors'
         endif
       end if

       write (*,*) 
       ! Call init_file_utils, ie. initialize the run_name, checking 
       !  whether we are doing a Trinity run or a list of runs.
       ! Until the logic of init_file_utils is fixed we set trin_run
       ! when an external file name has been provided... this prevents
       ! it overriding the name from the command line.
       if (state%run_name_external) then
          call init_file_utils (state%list, trin_run=state%run_name_external, &
             name=state%run_name, n_ensembles=state%nensembles)
       else
          call init_file_utils (state%list, name="gs")
       end if
       call debug_message(state%verb, &
         'gs2_main::initialize_gs2 initialized file_utils')
    end if
    
    call broadcast (state%list)
    call debug_message(state%verb, 'gs2_main::initialize_gs2 broadcasted list')

    if (state%list) then
       call job_fork
       call set_job_id(job)
       call debug_message(state%verb, 'gs2_main::initialize_gs2 called job fork')
    else if (state%nensembles > 1) then 
       call job_fork (n_ensembles=state%nensembles)
       call debug_message(state%verb, &
         'gs2_main::initialize_gs2 called job fork (nensembles>1)')
    end if

    if (proc0) call time_message(.false.,state%timers%total,' Total')


    ! Pass run name to other procs
    if (proc0) then
       run_name_target = trim(run_name)
    end if
    call broadcast (run_name_target)
    if (.not. proc0) run_name => run_name_target
    call debug_message(state%verb, 'gs2_main::initialize_gs2 run_name = '//run_name)

    !Set using_measure_scatter to indicate we want to use in "gather/scatter" timings
    using_measure_scatter=.false.

    !> Initialize the gs2 initialization system
    call init_gs2_init

    !state%gs2_initialized = .true.
    state%init%level = init_level_list%gs2

  end subroutine initialize_gs2

  subroutine initialize_equations(state)
    use fields, only: init_fields
    use geometry, only: surfarea, dvdrhon
    use gs2_time, only: init_tstart
    use gs2_init, only: init
    use init_g, only: tstart
    use job_manage, only: time_message
    use mp, only: proc0, broadcast
    use parameter_scan, only: init_parameter_scan
    use run_parameters, only: nstep, do_eigsolve
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
       !call time_message(.false., time_init,' Initialization')
    !call init_parameter_scan
    call init_parameter_scan

    if (proc0) call time_message(.false., state%timers%init,' Initialization')
    !Set using_measure_scatter to indicate we want to use in "gather/scatter" timings
    call debug_message(state%verb, 'gs2_main::initialize_equations calling init_fields')

    ! This triggers initializing of all the grids, all the physics parameters
    ! and all the modules which solve the equations
    call init(state%init, init_level_list%full)

    ! Set the initial simulation time (must be after init_fields
    ! because initial time may be read from a restart file)
    call init_tstart(tstart)

    ! Here we copy some geometric information required by 
    ! Trinity to state%outputs
    if (proc0) then
       state%outputs%dvdrho = dvdrhon
       state%outputs%grho = surfarea/dvdrhon
    end if
    call broadcast (state%outputs%dvdrho)
    call broadcast (state%outputs%grho)
    
    if (proc0) call time_message(.false.,state%timers%init,' Initialization')

    if (state%is_trinity_job) call write_trinity_parameters

    ! Set defaults. These are typically only important
    ! for the standalone gs2 program
    state%nstep = nstep
    state%do_eigsolve = do_eigsolve

    call debug_message(state%verb, 'gs2_main::initialize_equations finished')

    call allocate_outputs(state)

  end subroutine initialize_equations

  subroutine initialize_diagnostics(state)
    use gs2_diagnostics, only: init_gs2_diagnostics
    use gs2_diagnostics, only: nwrite, write_nl_flux
    use gs2_diagnostics, only: loop_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: diagnostics_init_options_type
    use gs2_diagnostics_new, only: init_gs2_diagnostics_new
    use gs2_diagnostics_new, only: run_diagnostics
#endif
    use parameter_scan, only: allocate_target_arrays
    use run_parameters, only: nstep
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
#ifdef NEW_DIAG
    ! Configuration for the new diagnostics module
    type(diagnostics_init_options_type) :: diagnostics_init_options
    real :: precision_test

#ifdef NETCDF_PARALLEL
    diagnostics_init_options%parallel_io_capable = .true.
#else
    diagnostics_init_options%parallel_io_capable = .false.
#endif

    ! Here we check if reals have been promoted to doubles
    diagnostics_init_options%default_double =  (precision(precision_test).gt.10)
    ! Check whether this is a Trinity run... enforces calculation of the
    ! fluxes
    diagnostics_init_options%is_trinity_run = state%is_trinity_job
    call init_gs2_diagnostics_new(diagnostics_init_options)
    ! Create variables and write constants
    call run_diagnostics(-1,state%exit)
    ! Write initial values
    call run_diagnostics(0,state%exit)
#endif

    call debug_message(state%verb, &
      'gs2_main::initialize_diagnostics calling init_gs2_diagnostics')
    call init_gs2_diagnostics (state%list, nstep)
    call allocate_target_arrays(nwrite,write_nl_flux) ! must be after init_gs2_diagnostics
    call loop_diagnostics(0,state%exit)

  end subroutine initialize_diagnostics

  subroutine evolve_equations(state, nstep_run)
    use collisions, only: vnmult
    use dist_fn_arrays, only: gnew
    use fields, only: advance
    use gs2_diagnostics, only: nsave, loop_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: run_diagnostics
#endif
    use gs2_reinit, only: reset_time_step
    use gs2_reinit, only: check_time_step
    use gs2_save, only: gs2_save_for_restart
    use gs2_time, only: user_time, user_dt, update_time, write_dt
    use job_manage, only: time_message, checkstop, checktime
    use mp, only: proc0
    use mp, only: scope, subprocs
    use parameter_scan, only: update_scan_parameter_value
    use run_parameters, only: reset, fphi, fapar, fbpar, nstep
    use run_parameters, only: avail_cpu_time, margin_cpu_time
    use unit_tests, only: ilast_step
    type(gs2_program_state_type), intent(inout) :: state
    integer :: istep, istatus
    integer, intent(in) :: nstep_run
    
    
    if (state%nensembles > 1) &
          call scope (subprocs)

    call time_message(.false.,state%timers%main_loop,' Main Loop')

    ! Make sure exit is false before entering
    ! timestep loop
    state%exit = .false.

    if (proc0) write (*,*) 'istep_end', state%istep_end

    ! We run for nstep_run iterations, starting from whatever istep we got
    ! to in previous calls to this function. Note that calling
    ! finalize_diagnostics resets state%istep_end
    do istep = state%istep_end, state%istep_end + nstep_run - 1

       if (istep .gt. nstep) then
         if (proc0) write (*,*) 'Reached maximum number of steps allowed &
           & (set by nstep) without restarting diagnostics.'
         exit
       end if


       if (proc0) call time_message(.false.,state%timers%advance,' Advance time step')

       !Initialise reset to true
       reset=.true.

       do while(reset)
          reset=.false. !So that we only do this once unless something triggers a reset
          call advance (istep)

          !If we've triggered a reset then actually reset
          if (reset) then
             ! if called within trinity, do not dump info to screen
             if (state%is_external_job) then
                call reset_time_step (state%init, istep, state%exit, state%external_job_id)
             else       
                call reset_time_step (state%init, istep, state%exit)
             end if
          end if
          if(state%exit) exit
       enddo
       
       if (nsave > 0 .and. mod(istep, nsave) == 0) &
            call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
       call update_time
       if(proc0) call time_message(.false.,state%timers%diagnostics,' Diagnostics')
       call loop_diagnostics (istep, state%exit)
       if(state%exit) state%converged = .true.
#ifdef NEW_DIAG
       call run_diagnostics (istep, state%exit)
#endif
       if(proc0) call time_message(.false.,state%timers%diagnostics,' Diagnostics')
       if (proc0) call time_message(.false.,state%timers%advance,' Advance time step')

       if(.not.state%exit)then

          !Note this should only trigger a reset for timesteps too small
          !as timesteps too large are already handled
          call check_time_step(reset,state%exit)
          call update_scan_parameter_value(istep, reset, state%exit)

          !If something has triggered a reset then reset here
          if (reset) then
             ! if called within trinity, do not dump info to screen
             if (state%is_external_job) then
                call reset_time_step (state%init, istep, state%exit, state%external_job_id)
             else       
                call reset_time_step (state%init, istep, state%exit)
             end if
          end if

          if ((mod(istep,5) == 0).and.(.not.state%exit)) call checkstop(state%exit)
          if (.not.state%exit) call checktime(avail_cpu_time,state%exit,margin_cpu_time)
       endif

       state%istep_end = istep

       if (state%exit) then
          exit
       end if
    end do

    call time_message(.false.,state%timers%main_loop,' Main Loop')

    if (proc0 .and. .not. state%is_external_job) call write_dt
    
    if (state%is_external_job) call print_times(state, state%timers)

    ilast_step = state%istep_end

  end subroutine evolve_equations

  subroutine run_eigensolver(state)
#ifdef WITH_EIG
    use eigval, only: init_eigval, finish_eigval, time_eigval
    use eigval, only: BasicSolve
#endif 
    use job_manage, only: time_message
    use mp, only: mp_abort
    type(gs2_program_state_type), intent(inout) :: state
#ifdef WITH_EIG
   !Start timer
   call time_message(.false.,time_eigval,' Eigensolver')

   !Initialise slepc and the default/input eigensolver parameters
   call init_eigval

   !Create a solver based on input paramters, use it to solve and
   !then destroy it.
   call BasicSolve

   !Tidy up
   call finish_eigval

   !Stop timer
   call time_message(.false.,time_eigval,' Eigensolver')
#else
   call mp_abort("Require slepc/petsc")
#endif

  end subroutine run_eigensolver

  subroutine finalize_diagnostics(state)
    use gs2_diagnostics, only: finish_gs2_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new 
#endif
    use job_manage, only: time_message
    use mp, only: proc0
    use parameter_scan, only: deallocate_target_arrays
    type(gs2_program_state_type), intent(inout) :: state

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

#ifdef NEW_DIAG
    call finish_gs2_diagnostics_new 
#endif
    call finish_gs2_diagnostics (state%istep_end)

    call deallocate_target_arrays

    state%istep_end = 1

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')
  end subroutine finalize_diagnostics

  subroutine finalize_equations(state)
    use gs2_init, only: init
    use job_manage, only: time_message
    use mp, only: proc0
    use parameter_scan, only: finish_parameter_scan
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

    call debug_message(state%verb, 'gs2_main::finalize_equations starting')

    call deallocate_outputs(state)
    call finish_parameter_scan
    call init(state%init, init_level_list%gs2)
    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

  end subroutine finalize_equations

  subroutine finalize_gs2(state)
    use file_utils, only: finish_file_utils
    use gs2_init, only: finish_gs2_init
    use job_manage, only: time_message
    use mp, only: finish_mp, proc0
    use mp, only: mp_abort
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (state%init%level .ne. init_level_list%gs2) then
      write  (*,*) "ERROR: Called finalize_gs2 at the &
        & wrong init_level (perhaps you have called finalize_gs2 &
        & without calling initialize_gs2, or without calling &
        & finalize_equations"
      stop 1
    end if

    !if ((.not. state%gs2_initialized) .or. &
        !state%equations_initialized .or. &
        !state%diagnostics_initialized) then
        !write (*,*) 'ERROR: initialize_gs2 can only be called when &
        !& gs2_initialized is true, and equations_initialized &
        !& and diagnostics_initialized, &
        !& are all false. '
       !stop 1
     !end if

    call finish_gs2_init

    if (proc0) call finish_file_utils

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

    if (proc0) call time_message(.false.,state%timers%total,' Total')
    
    call debug_message(state%verb, 'gs2_main::finalize_gs2 calling print_times')

    call print_times(state, state%timers)

    call debug_message(state%verb, 'gs2_main::finalize_gs2 calling finish_mp')

    if (.not. state%mp_comm_external) call finish_mp

    !state%gs2_initialized = .false.

    state%init%level = 0


    !if (.not. present(mpi_comm) .and. .not. nofin) call finish_mp
  end subroutine finalize_gs2

  subroutine override_geometry
  end subroutine override_geometry

  subroutine override_parameter_spec(state, parameter_label, species_index, val)
    use gs2_profile_overrides, only: otemp, odens, ofprim, otprim, ovnewk
    use gs2_init, only: init, init_level_list
    use species, only: op_spec => override_parameter
    use mp, only: mp_abort
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    integer, intent(in) :: parameter_label, species_index
    real, intent(in) :: val


    select case (parameter_label)
    case (otemp, odens, ofprim, otprim, ovnewk)
      call init(state%init, init_level_list%override_profiles)
      state%init%profile_overrides_set = .true.
      call op_spec(parameter_label, species_index, val)
    case default
      write (*,*) "Unknown parameter label: ", parameter_label
      call mp_abort("Unknown parameter label in override", .true.)
    end select



  end subroutine override_parameter_spec

  subroutine override_parameter_nospec(state, parameter_label, val)
    use gs2_profile_overrides, only: og_exb, omach
    use gs2_miller_geometry_overrides, only: orhoc,  oqval,  oshat,  orgeo_lcfs
    use gs2_miller_geometry_overrides, only: orgeo_local, okap, okappri, otri
    use gs2_miller_geometry_overrides, only: otripri, oshift, obetaprim
    use gs2_init, only: init, init_level_list
    use mp, only: mp_abort
    use dist_fn, only: op_dist_fn => override_parameter
    use theta_grid_params, only: op_theta_grid_params => override_parameter
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    integer, intent(in) :: parameter_label
    real, intent(in) :: val


    select case (parameter_label)
    case (og_exb, omach)
      call init(state%init, init_level_list%override_profiles)  
      state%init%profile_overrides_set = .true.
      call op_dist_fn(parameter_label, val)
    case (orhoc,  oqval,  oshat,  orgeo_lcfs, &
          orgeo_local, okap, okappri, otri, &
          otripri, oshift, obetaprim)
      call init(state%init, init_level_list%override_miller_geometry)  
      state%init%miller_geometry_overrides_set = .true.
      call op_theta_grid_params(parameter_label, val)
    case default
      write (*,*) "Unknown parameter label: ", parameter_label
      call mp_abort("Unknown parameter label in override", .true.)
    end select
  end subroutine override_parameter_nospec
  
  subroutine calculate_outputs(state)
    use gs2_diagnostics, only: start_time
    use gs2_diagnostics, only: pflux_avg, qflux_avg, heat_avg, vflux_avg !, start_time, nwrite, write_nl_flux
    use gs2_diagnostics, only: diffusivity
    use gs2_diagnostics, only: ensemble_average
    use gs2_time, only: user_time
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use run_parameters, only: trinity_linear_fluxes
    use species, only: nspec, spec
    use species, only: ions, electrons, impurity
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    real :: time_interval
    real :: diff
    integer :: is

    time_interval = user_time-start_time
      write (*,*) 'GETTING FLUXES'

    if (state%nensembles > 1) &
      call ensemble_average (state%nensembles, time_interval)

      
   ! Use simple gamma / k^2 estimates for transport in the 
   ! linear case. This is for testing Trinity
   if (trinity_linear_fluxes .and. &
         nonlinear_mode_switch .eq. nonlinear_mode_none) then
     diff = diffusivity()
     do is = 1,nspec
       ! Q = n chi grad T = n (gamma / k^2) dT / dr
       ! = dens  n_r (gamma_N v_thr / k_N**2 rho_r a) dT / drho drho/dr
       ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) T a / L_T drho/dr
       ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) temp T_r tprim drho/dr_N/a
       ! Q / (n_r  v_r T_r rho_r**2/a**2) 
       ! = dens (gamma_N / k_N**2) temp tprim grho
       !   
       ! grho factored to diffusivity in diagnostics
       qflux_avg(is) =  diff * spec(is)%dens * spec(is)%temp * spec(is)%tprim 
       pflux_avg(is) =  diff * spec(is)%dens**2.0 * spec(is)%fprim 
       heat_avg = 0.0
     end do
     !if (present(converged)) then
     !end if
   end if

   if (size(state%outputs%pflux) > 1) then
      state%outputs%pflux(1) = pflux_avg(ions)/time_interval
      state%outputs%qflux(1) = qflux_avg(ions)/time_interval
      state%outputs%heat(1) = heat_avg(ions)/time_interval
      state%outputs%pflux(2) = pflux_avg(electrons)/time_interval
      state%outputs%qflux(2) = qflux_avg(electrons)/time_interval
      state%outputs%heat(2) = heat_avg(electrons)/time_interval
      if (size(state%outputs%pflux) > 2) then
         state%outputs%pflux(3) = pflux_avg(impurity)/time_interval
         state%outputs%qflux(3) = qflux_avg(impurity)/time_interval
         state%outputs%heat(3) = heat_avg(impurity)/time_interval
      end if
   else
      state%outputs%pflux = pflux_avg/time_interval
      state%outputs%qflux = qflux_avg/time_interval
      state%outputs%heat = heat_avg/time_interval
   end if
   state%outputs%vflux = vflux_avg(1)/time_interval
  end subroutine calculate_outputs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine reset_timers(timers)
    type(gs2_timers_type), intent(inout) :: timers
    timers%init = 0.
    timers%advance = 0.
    timers%finish = 0.
    timers%total = 0. 
    timers%diagnostics=0.
    !timers%interval
    timers%main_loop = 0.
  end subroutine reset_timers

  subroutine print_times(state, timers)
    use mp, only: proc0
    use redistribute, only: time_redist
    use fields_arrays, only: time_field
    use gs2_reinit, only: time_reinit
    implicit none
    type(gs2_program_state_type), intent(in) :: state
    type(gs2_timers_type), intent(in) :: timers

    if (proc0) then
       if (state%is_external_job) then
          print '(/,'' Job ID:'', i4,'', total from timer is:'', 0pf9.2,'' min'',/)', &
               state%external_job_id, state%timers%total(1)/60.
       else
!    if (proc0 .and. .not. nofin) then

          print '(/,'' Initialization'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
#ifdef WITH_EIG
               &'' Eigensolver'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
#endif
               &'' Advance steps'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
               &''(redistribute'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %)'',/, &
               &''(field solve'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %)'',/, &
               &''(diagnostics'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %)'',/, &
               &'' Re-initialize'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
               &'' Finishing'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/,  &
               &'' total from timer is:'', 0pf9.2,'' min'',/)', &
               timers%init(1)/60.,timers%init(1)/timers%total(1), &
#ifdef WITH_EIG
               timers%eigval(1)/60.,timers%eigval(1)/timers%total(1), &
#endif
               timers%advance(1)/60.,timers%advance(1)/timers%total(1), &
               time_redist(1)/60.,time_redist(1)/timers%total(1), &
               time_field(1)/60.,time_field(1)/timers%total(1), &
               timers%diagnostics(1)/60.,timers%diagnostics(1)/timers%total(1), &
               time_reinit(1)/60.,time_reinit(1)/timers%total(1), &
               timers%finish(1)/60.,timers%finish(1)/timers%total(1),timers%total(1)/60.
       endif
    end if
  end subroutine print_times

  subroutine allocate_outputs(state)
    use species, only: nspec
    type(gs2_program_state_type), intent(inout) :: state
    allocate(state%outputs%pflux(nspec))
    allocate(state%outputs%qflux(nspec))
    allocate(state%outputs%heat(nspec))
  end subroutine allocate_outputs

  subroutine deallocate_outputs(state)
    type(gs2_program_state_type), intent(inout) :: state
    deallocate(state%outputs%pflux)
    deallocate(state%outputs%qflux)
  end subroutine deallocate_outputs







  !> This is the main subroutine in which gs2 is initialized, equations are advanced,
  !!   and the program is finalized.
  !! \section Basic Structure 
  !!  This subroutine broadly falls into 3 sections: 
  !! -# Initialisation -  allocate arrays, calculate the response matrix etc.
  !! -# Running -  a loop which runs for run_parameters::nstep time steps, unless the code
  !!  is prematurely halted either through an error, reaching the available 
  !!  time limit or manually  (by using the command $ touch run_name.stop).
  !! -# Finishing up:  writing out results, deallocating arrays.
  !! \section Details
  !! -# Initialisation
  !!  - Initialize message passing 
  !!  - Initialize timer 
  !!  - Report numer of processors being used 
  !!  - If it is a Trinity run then filename (the name of the input file) 
  !!  is passed to  init_file_utils
  !!  - Otherwise, figure out run name or get list of jobs 
  !!  - If given a list of jobs, fork  
  !! \section arguments Arguments
  !! All arguments are optional and are not used for gs2. 
  !! (EGH - used for Trinity?)



subroutine run_gs2 (mpi_comm, job_id, filename, nensembles, &
     pflux, qflux, vflux, heat, dvdrho, grho, trinity_reset, converged)
   use job_manage, only: trin_reset

   
   use job_manage, only: trin_restart, trin_reset
    use unit_tests, only: functional_test_flag, ilast_step
    implicit none

    integer, intent (in), optional :: mpi_comm, job_id, nensembles
    character (*), intent (in), optional :: filename
    real, dimension (:), intent (out), optional :: pflux, qflux, heat
    real, intent (out), optional :: vflux
    real, intent (out), optional :: dvdrho, grho
    logical, intent (in), optional :: trinity_reset
    logical, intent (out), optional :: converged

    logical :: first_time = .true.
    logical :: debug

    !debug = (verbosity() .gt. 2)


!
!CMR, 12/2/2010: 
!     add nofinish optional variable to avoid deallocations at end of simulation
!     as may want to do post-processing
!
!    if (present(nofinish)) nofin=nofinish
     

! HJL tests on Trinity optionals for load balancing
    old_iface_state%external_job_id = -1
    if(present(job_id)) then
      old_iface_state%external_job_id = job_id + 1
      old_iface_state%is_trinity_job = .true.
      old_iface_state%is_external_job = .true.
    end if
    if(present(trinity_reset)) then
       old_iface_state%converged = .false.
       trin_reset = trinity_reset
       first_time = .true.
       !if (debug) write (*,*) 'gs2_main::run_gs2 set first_time = .true., jid=', &
         !old_iface_state%external_job_id
    endif
    if (present(filename)) then
      old_iface_state%run_name_external = .true.
      old_iface_state%run_name = filename
    end if
    if (present(mpi_comm)) then
      old_iface_state%mp_comm_external = .true.
      old_iface_state%mp_comm = mpi_comm
    end if
    if (present(nensembles)) then
      old_iface_state%nensembles = nensembles
    end if
    if(present(trinity_reset)) trin_restart = .true. ! All trinity runs are restarted except the first

    if (first_time) then
      call initialize_wall_clock_timer
      call initialize_gs2(old_iface_state)
      call initialize_equations(old_iface_state)
      call initialize_diagnostics(old_iface_state)
      first_time = .false.
    endif !firstime


    if(old_iface_state%do_eigsolve)then
      call run_eigensolver(old_iface_state)
    else
      call evolve_equations(old_iface_state, old_iface_state%nstep)
    endif

    if (present(converged)) converged = old_iface_state%converged


    if (present(pflux) .or. present(dvdrho) .or. present(grho)) then
      call calculate_outputs(old_iface_state)
    end if
    if (present(dvdrho)) dvdrho = old_iface_state%outputs%dvdrho
    if (present(grho)) grho = old_iface_state%outputs%grho
    if (present(pflux)) then
      write (*,*) 'SETTING FLUXES'
      pflux = old_iface_state%outputs%pflux
      qflux = old_iface_state%outputs%qflux
      vflux = old_iface_state%outputs%vflux
      heat = old_iface_state%outputs%heat
    else
       if (.not. functional_test_flag ) &
         call finalize_diagnostics(old_iface_state) 
    end if

    if (.not. present(mpi_comm) .and. .not. functional_test_flag) then 
      call finalize_equations(old_iface_state)
      call finalize_gs2(old_iface_state)
    end if
    
  end subroutine run_gs2
  
! HJL <
  subroutine trin_finish_gs2

    call finalize_diagnostics(old_iface_state)
    call finalize_equations(old_iface_state)
    call finalize_gs2(old_iface_state)
    !use gs2_layouts, only: finish_layouts
    !use gs2_transforms, only: finish_transforms
    !use gs2_diagnostics, only: finish_gs2_diagnostics
    !use gs2_save, only: gs2_save_for_restart, finish_save
    !use theta_grid, only: finish_theta_grid
    !use unit_tests, only: ilast_step
!#ifdef NEW_DIAG
    !use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
!#endif
    
    !call finish_gs2_diagnostics (ilast_step)
!#ifdef NEW_DIAG
    !call finish_gs2_diagnostics_new
!#endif
    !call finish_gs2
!! HJL Species won't change during a run so shouldn't need this    
!!    call finish_trin_species

    !call finish_layouts
    !call finish_transforms
    !call finish_save
    !call finish_theta_grid

  end subroutine trin_finish_gs2
! > HJL


  subroutine finish_gs2
    
    use antenna, only: finish_antenna
    use collisions, only: finish_collisions
    use dist_fn, only: finish_dist_fn
    use fields, only: finish_fields
    use file_utils, only: finish_file_utils
    use hyper, only: finish_hyper
    use init_g, only: finish_init_g
    use kt_grids, only: finish_kt_grids
    use le_grids, only: finish_le_grids
    use parameter_scan, only: finish_parameter_scan
    use mp, only: proc0
    use nonlinear_terms, only: finish_nonlinear_terms
    use run_parameters, only: finish_run_parameters
    use species, only: finish_species
    use theta_grid, only: finish_theta_grid
    use gs2_transforms, only: finish_transforms
    use gs2_save, only: finish_gs2_save
    implicit none

    call finish_antenna
    call finish_collisions
    call finish_dist_fn
    call finish_fields
    call finish_hyper
    call finish_init_g
    call finish_theta_grid
    call finish_kt_grids
    call finish_le_grids
    call finish_nonlinear_terms
    call finish_run_parameters
    call finish_species
    call finish_parameter_scan
    call finish_transforms
    call finish_gs2_save
    if (proc0) call finish_file_utils

  end subroutine finish_gs2
  
  function gs2_main_unit_test_reset_gs2(fac)
    use species, only: spec, nspec
    use dist_fn, only: g_exb 
    use mp, only: mp_abort
    logical :: gs2_main_unit_test_reset_gs2
    real, intent(in) :: fac

    if (nspec.ne.2)  & 
      call mp_abort("gs2_main_unit_test_reset_gs2 only works with 2 species", .true.)

    write (*,*) 'HEREEEE'

    call reset_gs2(nspec, &
      ! Deliberately leave fac off the next 2 lines
      ! so QN is unaffected
      (/spec(1)%dens, spec(2)%dens/), &
      (/spec(1)%temp, spec(2)%temp/), &
      (/spec(1)%fprim, spec(2)%fprim/)*fac, &
      (/spec(1)%tprim, spec(2)%tprim/)*fac, &
      g_exb, 0.0, &
      (/spec(1)%vnewk, spec(2)%vnewk/), &
      1)
    write (*,*) 'HERAAAA'

    gs2_main_unit_test_reset_gs2 = .true.

  end function gs2_main_unit_test_reset_gs2

  subroutine reset_gs2 (ntspec, dens, temp, fprim, tprim, gexb, mach, nu, nensembles)

    use dist_fn, only: dist_fn_g_exb => g_exb
    use gs2_init, only: save_fields_and_dist_fn, init, init_level_list
    use fields, only: init_fields, f_reset => reset_init
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use gs2_diagnostics, only: gd_reset => reset_init
    use gs2_save, only: gs2_save_for_restart!, gs_reset => reset_init
    use species, only: reinit_species
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt, save_dt
    use mp, only: scope, subprocs, allprocs
    use mp, only: mp_abort
    use run_parameters, only: trinity_linear_fluxes
    use species, only: nspec, spec, impurity, ions, electrons
    use species, only: determine_species_order
    use gs2_profile_overrides, only: otprim, ofprim, otemp, odens, ovnewk, og_exb, omach

    implicit none

    integer, intent (in) :: ntspec, nensembles
    real, intent (in) :: gexb, mach
    real, dimension (:), intent (in) :: dens, fprim, temp, tprim, nu
    integer :: is
    integer :: refspec
    integer, dimension(ntspec) :: gs2_spec

    real :: dummy

    !integer :: istatus

    ! doing nothing with gexb or mach for now, but in future will need to when
    ! using GS2 to evolve rotation profiles in TRINITY
    ! EGH add a check for this.
    if (gexb .ne. dist_fn_g_exb) call mp_abort(&
      "ERROR: Changing g_exb in gs2_reset is not implemented yet.", .true.)
    if (ntspec .gt. nspec) call mp_abort(&
      "Cannot pass more species to reset_gs2 than nspec", .true.)
    ! To prevent compiler warnings
    dummy = mach
    dummy = gexb


    if (trinity_linear_fluxes.and.nonlinear_mode_switch.eq.nonlinear_mode_none) &
      call reset_linear_magnitude

    if (nensembles > 1) call scope (subprocs)

    call save_fields_and_dist_fn
    ! EGH is this line necessary?
    gnew = 0.
    call save_dt (code_dt)
    ! This call to gs2_diagnostics::reset_init sets some time averages
    ! and counters to zero... used mainly for trinity convergence checks.
    call gd_reset

    call determine_species_order
    if (nspec==1) then
      gs2_spec(1) = ions
      ! The species with the reference DENSITY
      ! NB ref temp is always ions (trin spec 1)
      refspec = 1
    else if (nspec==2) then 
      gs2_spec(1) = ions
      gs2_spec(2) = electrons
      ! For two species or more, ref dens is electrons (trin spec 2)
      refspec = 2
    else if (nspec==3) then
      gs2_spec(1) = ions
      gs2_spec(2) = electrons
      gs2_spec(3) = impurity
      refspec = electrons
    else
      call mp_abort("Can't handle more than 3 species in reset_gs2", .true.)
    end if

    do is = 1,nspec
      call override(old_iface_state, odens, gs2_spec(is), dens(is)/dens(refspec))
      ! Temp gs2_spec(is) always normalised to the main ions
      call override(old_iface_state, otemp, gs2_spec(is), temp(is)/temp(1))
      call override(old_iface_state, ofprim, gs2_spec(is), fprim(is))
      call override(old_iface_state, otprim, gs2_spec(is), tprim(is))
      call override(old_iface_state, ovnewk, gs2_spec(is), nu(is))
    end do
    !call override(old_iface_state, og_exb, gexb)
    ! mach doesn't do anything atm
    !call override(old_iface_state, omach, mach)

    call init(old_iface_state%init, init_level_list%full)
    old_iface_state%istep_end = 1
    if (nensembles > 1) call scope (allprocs)

  end subroutine reset_gs2

  subroutine reset_linear_magnitude
    use dist_fn_arrays, only: g, gnew
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew,  bparnew
    use run_parameters, only: fphi, fapar, fbpar
    use init_g, only: ginit
!    real :: norm
    logical :: dummy

!    if (fphi .gt. epsilon(0.0)) then
!      norm = maxval(real(conjg(phi)*phi))
!    elseif (fapar .gt. epsilon(0.0)) then
!      norm = maxval(real(conjg(apar)*apar))
!    elseif (fbpar .gt. epsilon(0.0)) then
!      norm = maxval(real(conjg(bpar)*bpar))
!    endif
!
!    norm = norm**0.5
!
!    phi = phi/norm
!    apar = apar/norm
!    bpar = bpar/norm
!    gnew = gnew/norm
!    g = g/norm
     
     if (fphi > 0.) then 
       phi = 0.0 
       phinew = 0.0
     end if
       
     if (fapar > 0.) then 
       apar = 0.0
       aparnew = 0.0
     end if
     if (fbpar > 0.) then 
       bpar = 0.0
       bparnew = 0.0
     end if
     g = 0.0; gnew = 0.0

     call ginit(dummy)
  end subroutine reset_linear_magnitude

  subroutine write_trinity_parameters
      use file_utils, only: open_output_file, close_output_file
      use theta_grid_params, only: write_theta_grid => write_trinity_parameters
      use species, only: write_species => write_trinity_parameters
      use run_parameters, only: write_run_parameters => write_trinity_parameters
      use mp, only: proc0
      integer :: trinpars_unit
      
      if (proc0) then
        call open_output_file(trinpars_unit, '.trinpars')
        call write_theta_grid(trinpars_unit)
        call write_species(trinpars_unit)
        call write_run_parameters(trinpars_unit)
        call close_output_file(trinpars_unit)
      end if
  end subroutine write_trinity_parameters

  subroutine gs2_trin_init (rhoc, qval, shat, rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift, &
       betaprim, ntspec, dens, temp, fprim, tprim, gexb, mach, nu, use_gs2_geo)

    use species, only: init_trin_species
    use theta_grid_params, only: init_trin_geo
    use dist_fn, only: dist_fn_g_exb => g_exb
    !use mp, only: broadcast

    implicit none

    integer, intent (inout) :: ntspec
    real, intent (inout) :: rhoc, qval, shat, rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift
    real, intent (inout) :: betaprim, gexb, mach
    real, dimension (:), intent (inout) :: dens, fprim, temp, tprim, nu
    logical, intent (inout) :: use_gs2_geo

    ! for now do nothing with gexb or mach, but need to include later if want to use GS2
    ! with TRINITY to evolve rotation profiles
    ! EGH add a check for this.
    if (gexb .ne. dist_fn_g_exb) then 
      write (*,*) "ERROR: Changing g_exb in gs2_reset is not implemented yet."
      ! Can't use mp_abort here as mp may not have been initialised
      stop 1
    end if

   !call broadcast(rhoc)
   !call broadcast(qval)
   !call broadcast(shat)
   !call broadcast(rgeo_lcfs)
   !call broadcast(rgeo_local)
   !call broadcast(kap)
   !call broadcast(kappri)
   !call broadcast(tri)
   !call broadcast(tripri)
   !call broadcast(shift)
 
   !call broadcast(betaprim)
   !call broadcast(ntspec)
   !call broadcast(dens)
   !call broadcast(temp)
   !call broadcast(fprim)
   !call broadcast(tprim)
   !call broadcast(gexb)
   !call broadcast(mach)
   !call broadcast(nu)
   !call broadcast(use_gs2_geo)
    call init_trin_species (ntspec, dens, temp, fprim, tprim, nu)
    if (.not. use_gs2_geo) call init_trin_geo (rhoc, qval, shat, &
         rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift, betaprim)
    
  end subroutine gs2_trin_init

# ifndef MAKE_LIB
end module gs2_main
# endif
