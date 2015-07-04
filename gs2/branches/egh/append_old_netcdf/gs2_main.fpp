module old_interface_store
    integer :: ntspec_store
    real :: rhoc_store, qval_store, shat_store, rgeo_lcfs_store, rgeo_local_store
    real :: kap_store, kappri_store, tri_store, tripri_store, shift_store
    real :: betaprim_store, gexb_store, mach_store
    real, dimension (:), allocatable :: dens_store, fprim_store, temp_store, tprim_store, nu_store
    logical :: use_gs2_geo_store
    logical :: override_profiles = .false.
    logical :: override_miller_geometry = .false.
end module old_interface_store

# ifndef MAKE_LIB

!> This module provides the external interface to gs2. It contains functions
!! to initialize gs2, functions to run gs2, functions to finalize gs2 and 
!! functions to override/tweak gs2 parameters. 
!!
!! gs2_main contains both the new and the old interface to gs2. The new interface
!! comprises the gs2_program_state_type object, and functions for manipulating
!! that object, which are documented individually. 
!!
!! The old interface
!! consists of the functions run_gs2, finish_gs2, trin_finish_gs2, reset_gs2
!! and gs2_trin_init. These functions have been reimplemented using the 
!! new interface. They are deprecated and are provided only for backwards
!! compatibility. 
!!
!! The basic flow of a gs2 program can be seen in gs2.fpp, reproduced here:
!!
!!     type(gs2_program_state_type) :: state
!!     call initialize_gs2(state)
!!     call initialize_equations(state)
!!     call initialize_diagnostics(state)
!!     if (state%do_eigsolve) then 
!!       call run_eigensolver(state)
!!     else
!!       call evolve_equations(state, state%nstep)
!!     end if
!!     call finalize_diagnostics(state)
!!     call finalize_equations(state)
!!     call finalize_gs2(state)
!!
!! You can manipulate the gs2_program_state before running gs2. A typical
!! manipulation is to pass in an external mpi communicator (which also 
!! means that MPI_Init is not called within initialize_gs2. 
!!
!!     state%mp_comm_external = .true
!!     state%mp_comm = my_mpi_communicator
!!     call initialize_gs2(state)
!!     ...
!! 
!! You can also override parameters: typically this is done either before
!! calling initialize_equations, or between two calls to evolve equations;
!! here we manipulate the temperature gradient of the first species:
!!
!!     call prepare_profiles_overrides(state)
!!     state%init%prof_ov%override_tprim(1) = .true.
!!     state%init%prof_ov%tprim(1) = 5.5
!!     call initialize_equations(state)
!!     call initialize_diagnostics(state)
!!     call evolve_equations(state, state%nstep/2)
!!     call prepare_profiles_overrides(state)
!!     state%init%prof_ov%tprim(1) = 6.0
!!     call initialize_equations(state)
!!     call evolve_equations(state, state%nstep/2)
!!     ...
!!   
!! It is very important to note that just because this interface is presented
!! in an object-oriented way, it does not, unfortunately, mean that the entire
!! program state, i.e. data arrays etc, is encapsulated in the gs2_program_state 
!! object. In fact most data is stored globally as module public variables. 
!! The gs2_program_state object is provided for convenience, as a way of
!! monitoring the execution state of gs2, and because it is hoped that in the 
!! future gs2 will, in fact, be thread safe and have proper data encapsulation.
!! A particular consequence of this is that only one instance of the
!! gs2_program_state object should exist on each process, i.e. in a given 
!! memory space.

module gs2_main
  use gs2_init, only: init_type, init_level_list
  use optimisation_config, only: optimisation_type
  implicit none
  public :: run_gs2, finish_gs2, reset_gs2, trin_finish_gs2

  public :: gs2_program_state_type

  public :: initialize_gs2
  public :: initialize_equations
  public :: initialize_diagnostics, evolve_equations, run_eigensolver
  public :: finalize_diagnostics, finalize_equations, finalize_gs2
  public :: gs2_trin_init
  public :: calculate_outputs
  public :: reset_equations

  public :: determine_gs2spec_from_trin
  public :: gs2spec_from_trin
  public :: densrefspec

  public :: prepare_miller_geometry_overrides
  public :: prepare_profiles_overrides
  public :: prepare_kt_grids_overrides
  public :: prepare_optimisations_overrides
  public :: prepare_initial_values_overrides
  public :: set_initval_overrides_to_current_vals
  public :: finalize_overrides

  public :: old_iface_state
  !> Unit tests

  public :: gs2_main_unit_test_reset_gs2

  public :: initialize_wall_clock_timer


  type gs2_timers_type
    real :: init(2) 
    real :: advance(2) = 0.
    real :: timestep(2) = 0.
    real :: finish(2) = 0.
    real :: total(2) = 0. 
    real :: diagnostics(2)=0.
    !real :: interval
    real :: main_loop(2)
  end type gs2_timers_type

  !> A type for storing outputs of gs2
  !! for access externally
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


  !> The object which specifies and records the gs2 program state.
  !! Some attributes of the object, like mp_comm_external, are
  !! designed to directly manipulated by the user, and some are
  !! designed to store program state information and be 
  !! manipulated internally, like init\%level.
  type gs2_program_state_type

    ! Flags indicating the current state of the 
    ! program (used for error checking)
    logical :: gs2_initialized = .false.
    logical :: equations_initialized = .false.
    logical :: diagnostics_initialized = .false.

    !> A type for keeping track of the current
    !! initialization level of gs2, as well
    !! as storing all the overrides. See 
    !! gs2_init::init_type for more information.
    type(init_type) :: init


    !> Do not set manually. If fewer than the 
    !! available number of processors are being
    !! used, this is true for the processors
    !! that are active and false for those that
    !! lie idle.
    logical :: included = .true.


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

    !> This parameter is set equal to run_parameters::nstep in 
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
    integer :: istep_end = 0

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

    !> If true, don't initialize or evolve the equations, 
    !! just read their values from an existing output file.
    !! This allows additional calculations to be done, and 
    !! also allows GS2 to cheaply reproduce its functionality,
    !! e.g. within Trinity.
    logical :: replay = .false.

    !> If true, the equations are fully initialized during
    !! replays. This is typically needed for linear diffusive
    !! flux estimates, and is set true automatically in that 
    !! case.
    logical :: replay_full_initialize = .false.

    !> If true, 
    !! print full timing breakdown. 
    logical :: print_full_timers = .false.

    !> If true, print run time or
    !! full timing breakdown, depending
    !! on the value of print_full_timers
    logical :: print_times = .true.


    !> Parameters to be used when passing in an external communicator
    logical :: mp_comm_external = .false.
    integer :: mp_comm

    !> This is set in initialize_gs2 to the number
    !! of procs actually used
    integer :: nproc_actual

    logical :: run_name_external = .false.
    character(2000) :: run_name

    logical :: skip_diagnostics = .false.
    logical :: dont_change_timestep = .false.

    !> Whether this is a list mode run
    logical :: list
    !> The number of identical runs happening
    !! simultaneously (used for ensemble averaging).
    !! Cannot be used in conjunction with list mode
    integer :: nensembles = 1

    ! Outputs (e.g. for Trinity)
    type(gs2_outputs_type) :: outputs

    ! Optimisation configuration
    type(optimisation_type) :: optim

   
  end type gs2_program_state_type

  private


  integer :: densrefspec
  integer, dimension(:), allocatable :: gs2spec_from_trin


  !> This object is used for implementing the old interface
  !! and should not be modified
  type(gs2_program_state_type) :: old_iface_state
 
contains
# endif

  !> Starts the global wall clock timer
  !! used by check time. This is useful
  !! if you have multiple copies of gs2 
  !! running but you don't want to start 
  !! them all at the same time, but you
  !! want them all to respect avail_cpu_time
  subroutine initialize_wall_clock_timer
    use job_manage, only: init_checktime
    call init_checktime
  end subroutine initialize_wall_clock_timer

  !> Initialise message passing, open the input file, split the 
  !! communicator if running in list mode or if nensembles > 1,
  !! initialize the timers. After calling this function, gs2
  !! reaches init\%level = basic. If it is desired to provide
  !! an external commuicator or set the input file name externally,
  !! these should be set before calling this function. 
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
    use mp, only: use_nproc, included, mp_comm
    use redistribute, only: using_measure_scatter
    use run_parameters, only: avail_cpu_time
    use runtime_tests, only: verbosity
    use unit_tests, only: debug_message, set_job_id
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    integer :: oldrank

    if (state%init%level .ge. init_level_list%basic) then
      write  (*,*) "ERROR: Called initialize_gs2 twice &
        & without calling finalize_gs2"
      stop 1
    end if



    call debug_message(4, 'gs2_main::initialize_gs2 starting initialization')

    if (state%mp_comm_external) then
       call init_mp (state%mp_comm)
    else
       call init_mp
    end if

    oldrank = iproc
    if (state%init%opt_ov%override_nproc) then
      state%init%opt_ov%old_comm = mp_comm
      call use_nproc(state%init%opt_ov%nproc)
    end if
    state%included = included
    state%nproc_actual = nproc
    if (.not. state%included) return
    !write (*,*) 'I am included and my rank is ', iproc, 'and my old rank is', oldrank

    if (state%is_trinity_job) state%is_external_job = .true.
    if (state%is_external_job) then 
      call broadcast(state%external_job_id)
      call set_job_id(state%external_job_id)
    end if

    ! All subroutines that use this are now no longer 
    ! needed and will be deleted
    !trin_flag = state%is_trinity_job
    trin_flag = .false. 

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

    ! Initialize the gs2 initialization system
    call init_gs2_init(state%init)

    !state%gs2_initialized = .true.
    state%init%level = init_level_list%basic

    if (proc0) call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, 'gs2_main::initialize_gs2 finished')

  end subroutine initialize_gs2

  !> Initialize all the modules which are used to evolve the 
  !! equations. After calling this function, gs2 reaches 
  !! init\%level = full. 
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
    use gs2_reinit, only: init_gs2_reinit
    use fields_implicit, only: skip_initialisation
    !use dist_fn, only: replay_d=>replay
    use fields, only: fields_allocate_arrays=>allocate_arrays
    use fields_local, only: fieldmat
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use run_parameters, only: trinity_linear_fluxes
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return
    call debug_message(state%verb, 'gs2_main::initialize_equations starting')
       !call time_message(.false., time_init,' Initialization')
    !call init_parameter_scan
    if (proc0) call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, &
      'gs2_main::initialize_equations calling init_parameter_scan')
    call init_parameter_scan

    if (proc0) call time_message(.false., state%timers%init,' Initialization')
    !Set using_measure_scatter to indicate we want to use in "gather/scatter" timings
    call debug_message(state%verb, 'gs2_main::initialize_equations calling init_fields')

    !> If we are doing a replay then we are not evolving the equations
    !! so we disable calculating the implicit response
    !! matrix
    if (state%replay) then
      skip_initialisation = .true.
      !replay_d = .true.
      fieldmat%no_prepare = .true.
      fieldmat%no_populate = .true.
      call init(state%init, init_level_list%le_grids)

      if (trinity_linear_fluxes) then
        state%replay_full_initialize = .true.
        call init(state%init, init_level_list%full)
      end if
        
      !call fields_allocate_arrays
      !nonlinear_mode_switch = nonlinear_mode_none
    else
      ! This triggers initializing of all the grids, all the physics parameters
      ! and all the modules which solve the equations
      call init(state%init, init_level_list%full)
    end if

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


    ! Set defaults. These are typically only important
    ! for the standalone gs2 program
    state%nstep = nstep
    state%do_eigsolve = do_eigsolve

    call debug_message(state%verb, 'gs2_main::initialize_equations finished')

    ! Since the number of species may have changed
    ! since the last call to initialize_equations 
    ! we reallocate outputs every time.
    if (associated(state%outputs%pflux)) then
      call deallocate_outputs(state)
    end if
    call allocate_outputs(state)

    call init_gs2_reinit

    if (proc0) call time_message(.false.,state%timers%total,' Total')

  end subroutine initialize_equations

  !> Initialize the diagostics modules. This function can only be
  !! called when the init level is 'full', i.e., after calling 
  !! initialize_equations.  However, gs2 can be partially
  !! uninitialized without finalizing the diagnostics, for example to
  !! change parameters such as the timestep.
  subroutine initialize_diagnostics(state)
    use gs2_diagnostics, only: init_gs2_diagnostics
    use gs2_diagnostics, only: nwrite, write_nl_flux
    use gs2_diagnostics, only: loop_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: diagnostics_init_options_type
    use gs2_diagnostics_new, only: init_gs2_diagnostics_new
    use gs2_diagnostics_new, only: run_diagnostics
#endif
    use job_manage, only: time_message
    use mp, only: proc0, mp_abort, broadcast
    use parameter_scan, only: allocate_target_arrays
    use run_parameters, only: nstep, use_old_diagnostics
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
#ifdef NEW_DIAG
    ! Configuration for the new diagnostics module
    type(diagnostics_init_options_type) :: diagnostics_init_options
    real :: precision_test
#endif
    if (.not. state%included) return

    state%exit = .false. 

    if (proc0) call time_message(.false.,state%timers%total,' Total')

    if (state%init%diagnostics_initialized) &
      call mp_abort('Calling initialize_diagnostics twice', .true.)

    if (.not. use_old_diagnostics) then
#ifdef NEW_DIAG
      call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics calling init_gs2_diagnostics_new')

#ifdef NETCDF_PARALLEL
      diagnostics_init_options%parallel_io_capable = .true.
      diagnostics_init_options%parallel_io_capable = .false.
#else
      diagnostics_init_options%parallel_io_capable = .false.
#endif

      ! Here we check if reals have been promoted to doubles
      diagnostics_init_options%default_double =  (precision(precision_test).gt.10)
      ! Check whether this is a Trinity run... enforces calculation of the
      ! fluxes
      diagnostics_init_options%is_trinity_run = state%is_trinity_job
      ! Specifiy whether we are doing a replay.
      diagnostics_init_options%replay = state%replay

      call init_gs2_diagnostics_new(diagnostics_init_options)

      call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics calling run_diagnostics')
      ! Create variables and write constants
      call run_diagnostics(-1,state%exit, state%replay)
      !call broadcast(state%exit)
      ! Write initial values
      call run_diagnostics(0,state%exit, state%replay)
#else
      call mp_abort("use_old_diagnostics is .false. but you have &
        & not built gs2 with new diagnostics enabled", .true.)
#endif
    else ! if use_old_diagnostics

      call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics calling init_gs2_diagnostics')
      call init_gs2_diagnostics (state%list, nstep)
      call allocate_target_arrays(nwrite,write_nl_flux) ! must be after init_gs2_diagnostics
      call loop_diagnostics(0,state%exit)
    end if

    state%init%diagnostics_initialized = .true.

    if (proc0) call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics finished')

  end subroutine initialize_diagnostics



  !> Run the initial value solver. nstep_run must
  !! be less than or equal to state\%nstep, which is 
  !! set from the input file. The cumulative number of
  !! steps that can be run cannot exceed state\%nstep, 
  !! without a call to reset_equations.
  !! Examples:
  !!
  !!      call evolve_equations(state, state%nstep) ! Basic run
  !!
  !! Note that these two calls: 
  !!  
  !!      call evolve_equations(state, state%nstep/2) 
  !!      call evolve_equations(state, state%nstep/2)
  !!
  !! will in general have the identical effect to the single call above.
  !! There are circumstances when this is not true,
  !! for example, if the first of the two calls to evolve_equations  
  !! exits without running nstep/2 steps (perhaps because the
  !! growth rate has converged). 
  !!
  !! This example will cause an error because the total number
  !! of steps exceeds state\%nstep:
  !!
  !!      call evolve_equations(state, state%nstep/2) 
  !!      call evolve_equations(state, state%nstep/2)
  !!      call evolve_equations(state, state%nstep/2)
  !!
  !! This is OK:
  !!
  !!      call evolve_equations(state, state%nstep/2) 
  !!      call evolve_equations(state, state%nstep/2)
  !!      call reset_equations(state)
  !!      call evolve_equations(state, state%nstep/2)

  subroutine evolve_equations(state, nstep_run)
    use collisions, only: vnmult
    use dist_fn_arrays, only: gnew
    use fields, only: advance
    use gs2_diagnostics, only: nsave, loop_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: run_diagnostics
    use gs2_diagnostics_new, only: gnostics
#endif
    use gs2_reinit, only: reset_time_step
    use gs2_reinit, only: check_time_step
    use gs2_save, only: gs2_save_for_restart
    use gs2_time, only: user_time, user_dt, update_time, write_dt
    use job_manage, only: time_message, checkstop, checktime
    use mp, only: proc0
    use mp, only: scope, subprocs, broadcast
    use parameter_scan, only: update_scan_parameter_value
    use run_parameters, only: reset, fphi, fapar, fbpar, nstep
    use run_parameters, only: avail_cpu_time, margin_cpu_time
    use run_parameters, only: use_old_diagnostics
    use unit_tests, only: ilast_step, debug_message
    use gs2_init, only: init
    type(gs2_program_state_type), intent(inout) :: state
    integer :: istep, istatus
    integer, intent(in) :: nstep_run
    logical :: temp_initval_override_store
    integer :: istep_loop_max

    if (.not. state%included) return

    call debug_message(state%verb, &
        'gs2_main::evolve_equations starting')
    
    temp_initval_override_store = state%init%initval_ov%override
    
    if (proc0) call time_message(.false.,state%timers%total,' Total')
    
    if (state%nensembles > 1) &
          call scope (subprocs)

    if (state%is_trinity_job) call write_trinity_parameters

    call time_message(.false.,state%timers%main_loop,' Main Loop')

    ! Make sure exit is false before entering
    ! timestep loop
    state%exit = .false.

    !if (proc0) write (*,*) 'istep_end', state%istep_end

    call debug_message(state%verb, &
        'gs2_main::evolve_equations starting loop')
    ! We run for nstep_run iterations, starting from whatever istep we got
    ! to in previous calls to this function. Note that calling
    ! finalize_diagnostics resets state%istep_end
    istep_loop_max = state%istep_end + nstep_run
    do istep = state%istep_end+1, istep_loop_max

       if (istep .gt. nstep) then
         if (proc0) write (*,*) 'Reached maximum number of steps allowed &
           & (set by nstep) without restarting diagnostics.'
         exit
       end if


       if (proc0) call time_message(.false.,state%timers%advance,' Advance time step')
       !if (proc0) write (*,*) 'advance 1', state%timers%advance

       ! If we are doing a replay we don't actually evolve 
       ! the equations
       if (.not. state%replay) then 
         !Initialise reset to true
         reset=.true.

         call debug_message(state%verb+1, &
            'gs2_main::evolve_equations calling advance')
         do while(reset)
            reset=.false. !So that we only do this once unless something triggers a reset
            if (proc0) call time_message(.false.,state%timers%timestep,' Timestep')
            call advance (istep)
            if (proc0) call time_message(.false.,state%timers%timestep,' Timestep')
            call debug_message(state%verb+1, &
              'gs2_main::evolve_equations called advance')

            if(state%dont_change_timestep) reset = .false.
            !If we've triggered a reset then actually reset
            if (reset) then
               call prepare_initial_values_overrides(state)
               call debug_message(state%verb+1, &
                  'gs2_main::evolve_equations resetting timestep')
               call set_initval_overrides_to_current_vals(state%init%initval_ov)
               state%init%initval_ov%override = .true.
               if (state%is_external_job) then
                  call reset_time_step (state%init, istep, state%exit, state%external_job_id)
               else       
                  call reset_time_step (state%init, istep, state%exit)
               end if
            end if
            if(state%exit) exit
         enddo
         call debug_message(state%verb+1, &
            'gs2_main::evolve_equations calling gs2_save_for_restart')
         
         if (use_old_diagnostics) then
           if (nsave > 0 .and. mod(istep, nsave) == 0) &
                call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
         else
#ifdef NEW_DIAG
           if (gnostics%nsave > 0 .and. mod(istep, gnostics%nsave) == 0) &
                call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
#endif
         end if
         call update_time
         if(proc0) call time_message(.false.,state%timers%diagnostics,' Diagnostics')

       end if ! if (.not. state%replay)
       call debug_message(state%verb+1, &
          'gs2_main::evolve_equations calling diagnostics')
       if (.not. state%skip_diagnostics) then 
         if (use_old_diagnostics) then 
           call loop_diagnostics (istep, state%exit)
#ifdef NEW_DIAG
         else 
           !if (state%replay) call broadcast(state%exit)
           call run_diagnostics (istep, state%exit, state%replay)
#endif
         end if
       end if
       !if (state%replay) then 
         !!call debug_message(state%verb+1, &
           !!'gs2_main::evolve_equations broadcasting exit')
         !! Necessary with replay because other procs just read garbage and 
         !! exit may be wrong
         !!call broadcast(state%exit)
!#ifdef NEW_DIAG
         !if (.not. state%exit .and. istep < istep_loop_max ) then
           !call debug_message(state%verb+1, &
             !'gs2_main::evolve_equations rereading variable values')
           !call run_diagnostics (istep, state%exit, .true.)
         !end if
!#endif          
         !call broadcast(state%exit)
       !end if

       if(state%exit) state%converged = .true.
       if (state%exit) call debug_message(state%verb-1, &
         'gs2_main::evolve_equations exit true after diagnostics')


       if(proc0) call time_message(.false.,state%timers%diagnostics,' Diagnostics')

       if (proc0) call time_message(.false.,state%timers%advance,' Advance time step')

       if(.not.state%exit)then
          if (.not.state%replay) then
            call debug_message(state%verb+1, &
              'gs2_main::evolve_equations checking time step')

            !Note this should only trigger a reset for timesteps too small
            !as timesteps too large are already handled
            call check_time_step(reset,state%exit)
            !call update_scan_parameter_value(istep, reset, state%exit)
            call debug_message(state%verb-1, &
              'gs2_main::evolve_equations checked time step')

            !If something has triggered a reset then reset here
            if(state%dont_change_timestep) reset = .false.
            if (reset) then
               call prepare_initial_values_overrides(state)
               call set_initval_overrides_to_current_vals(state%init%initval_ov)
               state%init%initval_ov%override = .true.
               ! if called within trinity, do not dump info to screen
               if (state%is_external_job) then
                  call reset_time_step (state%init, istep, state%exit, state%external_job_id)
               else       
                  call reset_time_step (state%init, istep, state%exit)
               end if
            end if
          end if ! if (.not. state%replay)

          if ((mod(istep,5) == 0).and.(.not.state%exit)) call checkstop(state%exit)
          if (state%exit) call debug_message(state%verb-1, &
            'gs2_main::evolve_equations exit true after checkstop')
          if (.not.state%exit) call checktime(avail_cpu_time,state%exit,margin_cpu_time)
          if (state%exit) call debug_message(state%verb-1, &
            'gs2_main::evolve_equations exit true after checktime')
       endif

       state%istep_end = istep
       !if (proc0) write (*,*) 'advance 2', state%timers%advance, state%exit


       if (state%exit) then
           call debug_message(state%verb, &
                'gs2_main::evolve_equations exiting loop')
          exit
       end if

    end do

    call time_message(.false.,state%timers%main_loop,' Main Loop')

    if (proc0 .and. .not. state%is_external_job) call write_dt
    
    if (proc0) call time_message(.false.,state%timers%total,' Total')

    if (state%print_times) call print_times(state, state%timers)

    ilast_step = state%istep_end

    ! At the end evolve_equations, we guarantee
    ! that if the initial value overrides have
    ! been initialized then they are 
    ! set to the current value (otherwise
    ! this might vary depending on whether the 
    ! timestep has been reset).
    !
    ! We also restore the value of the override
    ! switch to what it was when this function was
    ! called.
    ! 
    call debug_message(state%verb, &
      'gs2_main::evolve_equations restoring initial value override')
    state%init%initval_ov%override = temp_initval_override_store
    if (state%init%initval_ov%init) then
      call set_initval_overrides_to_current_vals(state%init%initval_ov)
    end if

    call debug_message(state%verb,'gs2_main::evolve_equations finished')


  end subroutine evolve_equations

  subroutine run_eigensolver(state)
#ifdef WITH_EIG
    use eigval, only: init_eigval, finish_eigval, time_eigval
    use eigval, only: BasicSolve
#endif 
    use job_manage, only: time_message
    use mp, only: mp_abort, proc0
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
#ifdef WITH_EIG
   if (proc0) call time_message(.false.,state%timers%total,' Total')

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
    if (proc0) call time_message(.false.,state%timers%total,' Total')

  end subroutine run_eigensolver

  !subroutine save_current(state)
    !use gs2_init, only: save_fields_and_dist_fn
    !use mp, only: scope, subprocs, allprocs
    !implicit none
    !type(gs2_program_state_type), intent(inout) :: state


    !if (state%nensembles > 1) call scope (subprocs)
    !call save_fields_and_dist_fn
    !if (state%nensembles > 1) call scope (allprocs)
  !end subroutine save_current

  subroutine reset_equations (state)

    use gs2_diagnostics, only: gd_reset => reset_init
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: reset_averages_and_counters
#endif
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use run_parameters, only: trinity_linear_fluxes, use_old_diagnostics
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt, save_dt
    use mp, only: scope, subprocs, allprocs
    use mp, only: mp_abort
    use unit_tests, only: debug_message

    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return


    if (state%nensembles > 1) call scope (subprocs)

    call debug_message(state%verb, 'gs2_main::reset_equations starting')

    ! EGH is this line necessary?
    !gnew = 0.
    call save_dt (code_dt)
    ! This call to gs2_diagnostics::reset_init sets some time averages
    ! and counters to zero... used mainly for trinity convergence checks.
    if (use_old_diagnostics) then
      call gd_reset
    else 
#ifdef NEW_DIAG
      call reset_averages_and_counters
#endif
    end if

    if (trinity_linear_fluxes .and. &
       nonlinear_mode_switch .eq. nonlinear_mode_none) &
        call reset_linear_magnitude

    state%istep_end = 0

    call debug_message(state%verb, 'gs2_main::reset_equations finished')

    !write (*,*) 'nensembles = ', state%nensembles

    if (state%nensembles > 1) call scope (allprocs)

    !write (*,*) 'here at the close....'

  end subroutine reset_equations

  subroutine finalize_diagnostics(state)
    use gs2_diagnostics, only: finish_gs2_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new 
#endif
    use job_manage, only: time_message
    use mp, only: proc0, mp_abort
    use parameter_scan, only: deallocate_target_arrays
    use run_parameters, only: use_old_diagnostics
    use unit_tests, only: debug_message
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return

    if (proc0) call time_message(.false.,state%timers%total,' Total')

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

    call debug_message(state%verb, 'gs2_main::finalize_diagnostics starting')

    if (.not. state%init%diagnostics_initialized) &
      call mp_abort('Calling finalize_diagnostics when &
      & diagnostics are not initialized', .true.)

    if (.not. use_old_diagnostics) then
#ifdef NEW_DIAG
      call finish_gs2_diagnostics_new 
#endif
    else
      call finish_gs2_diagnostics (state%istep_end)
    end if

    call deallocate_target_arrays

    state%istep_end = 0

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

    state%init%diagnostics_initialized = .false.

    if (proc0) call time_message(.false.,state%timers%total,' Total')
  end subroutine finalize_diagnostics

  subroutine finalize_equations(state)
    use gs2_init, only: init
    use job_manage, only: time_message
    use mp, only: proc0
    use parameter_scan, only: finish_parameter_scan
    use unit_tests, only: debug_message
    use gs2_reinit, only: finish_gs2_reinit
    !use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')
    if (proc0) call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, 'gs2_main::finalize_equations starting')

    call finish_gs2_reinit

    call deallocate_outputs(state)
    call finish_parameter_scan
    !if (state%replay .and. allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew)

    call init(state%init, init_level_list%basic)
    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')
    if (proc0) call time_message(.false.,state%timers%total,' Total')

  end subroutine finalize_equations

  subroutine finalize_gs2(state)
    use file_utils, only: finish_file_utils
    use gs2_init, only: finish_gs2_init
    use job_manage, only: time_message
    use mp, only: finish_mp, proc0, barrier
    use mp, only: mp_abort, unsplit_all, included, nproc, iproc
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (state%included) then

      if (state%init%level .ne. init_level_list%basic) then
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

      if (proc0) call time_message(.false.,state%timers%total,' Total')

      call finish_gs2_init(state%init)

      if (proc0) call finish_file_utils

      if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

      if (proc0) call time_message(.false.,state%timers%total,' Total')
      
      call debug_message(state%verb, 'gs2_main::finalize_gs2 calling print_times')

      if (state%print_times) call print_times(state, state%timers)

    end if

    !write (*,*) 'I AM iproc ', iproc, 'and nproc = ', nproc
    if (state%init%opt_ov%override_nproc) &
      call unsplit_all(state%init%opt_ov%old_comm)
    call barrier
    state%included = included
    !write (*,*) 'I am iproc ', iproc, 'and nproc = ', nproc

    call debug_message(state%verb, 'gs2_main::finalize_gs2 calling finish_mp')

    if (.not. state%mp_comm_external) call finish_mp

    !state%gs2_initialized = .false.

    state%init%level = 0


    !if (.not. present(mpi_comm) .and. .not. nofin) call finish_mp
  end subroutine finalize_gs2

  subroutine prepare_optimisations_overrides(state)
    use overrides, only: init_optimisations_overrides
    use gs2_init, only: init, init_level_list
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    !call init(state%init, init_level_list%override_optimisations-1)
    call init_optimisations_overrides(state%init%opt_ov)
  end subroutine prepare_optimisations_overrides

  subroutine prepare_miller_geometry_overrides(state)
    use overrides, only: init_miller_geometry_overrides
    use gs2_init, only: init, init_level_list
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_miller_geometry-1)
    call init_miller_geometry_overrides(state%init%mgeo_ov)
  end subroutine prepare_miller_geometry_overrides

  subroutine prepare_profiles_overrides(state)
    use overrides, only: init_profiles_overrides
    use gs2_init, only: init, init_level_list
    use species, only: nspec
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_profiles-1)
    call init_profiles_overrides(state%init%prof_ov, nspec)
  end subroutine prepare_profiles_overrides

  subroutine prepare_kt_grids_overrides(state)
    use overrides, only: init_kt_grids_overrides
    use gs2_init, only: init, init_level_list
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_kt_grids-1)
    call init_kt_grids_overrides(state%init%kt_ov)
  end subroutine prepare_kt_grids_overrides

  subroutine prepare_initial_values_overrides(state)
    use overrides, only: init_initial_values_overrides
    use fields, only: force_maxwell_reinit
    use gs2_init, only: in_memory
    use gs2_init, only: init, init_level_list
    use gs2_layouts, only: g_lo
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use species, only: nspec
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_initial_values-1)
    call init_initial_values_overrides(state%init%initval_ov,&
      ntgrid, ntheta0, naky, g_lo%llim_proc, g_lo%ulim_alloc, &
      force_maxwell_reinit=force_maxwell_reinit, &
      in_memory=in_memory)
  end subroutine prepare_initial_values_overrides

  subroutine set_initval_overrides_to_current_vals(initval_ov)
    use dist_fn_arrays, only: gnew, g_restart_tmp
    use gs2_save, only: gs2_save_for_restart
    use mp, only: proc0, broadcast, mp_abort
    use collisions, only: vnmult
    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use file_utils, only: error_unit
    use kt_grids, only: ntheta0, naky
    use run_parameters, only: fphi, fapar, fbpar
    use fields, only:  force_maxwell_reinit
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_time, only: user_time, user_dt
    use overrides, only: initial_values_overrides_type
    use antenna, only: dump_ant_amp
    implicit none
    type(initial_values_overrides_type), intent(inout) :: initval_ov
    integer :: iostat
    integer :: istatus
    !if (.not. state%included) return

    if (.not.initval_ov%init) &
      call mp_abort("Trying to set initial value overrides &
      & before they are initialized... have you called &
      & prepare_initial_values_overrides ? ", .true.)

    if(initval_ov%in_memory)then
      initval_ov%g=gnew
      if (.not. initval_ov%force_maxwell_reinit) then
        if(fphi.gt.0) then
          initval_ov%phi=phinew
        end if
        if(fapar.gt.0) then 
          initval_ov%apar=aparnew
        end if
        if(fbpar.gt.0) then 
          initval_ov%bpar=bparnew
        end if
      end if
    else ! if(.not.in_memory)then
      !Should really do this with in_memory=.true. as well but
      !not sure that we really need to as we never read in the dumped data.
      if (proc0) call dump_ant_amp

      call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    endif

  end subroutine set_initval_overrides_to_current_vals

  subroutine finalize_overrides(state)
    use overrides, only: finish_profiles_overrides
    use overrides, only: finish_miller_geometry_overrides
    use overrides, only: finish_initial_values_overrides
    type(gs2_program_state_type), intent(inout) :: state

    if (state%init%mgeo_ov%init) call &
      finish_miller_geometry_overrides(state%init%mgeo_ov)
    if (state%init%prof_ov%init) call &
      finish_profiles_overrides(state%init%prof_ov)
    if (state%init%initval_ov%init) call &
      finish_initial_values_overrides(state%init%initval_ov)
  end subroutine finalize_overrides

  
  subroutine calculate_outputs(state)
    use gs2_diagnostics, only: start_time
    use gs2_diagnostics, only: pflux_avg, qflux_avg, heat_avg, vflux_avg !, start_time, nwrite, write_nl_flux
    use gs2_diagnostics, only: diffusivity
    use gs2_diagnostics, only: ensemble_average
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: gnostics
#endif 
    use gs2_time, only: user_time
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use run_parameters, only: trinity_linear_fluxes
    use species, only: nspec, spec
    use species, only: ions, electrons, impurity
    use run_parameters, only: use_old_diagnostics
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    real :: time_interval
    real :: diff
    real, dimension(nspec) :: qf, pf, ht, vf
    integer :: is

    if (.not. state%included) return

    !write (*,*) 'GETTING FLUXES', 'user_time', user_time, start_time, 'DIFF', time_interval

    if (state%nensembles > 1) &
      call ensemble_average (state%nensembles, time_interval)

      
   ! Use simple gamma / k^2 estimates for transport in the 
   ! linear case. This is for testing Trinity
   if (trinity_linear_fluxes .and. &
         nonlinear_mode_switch .eq. nonlinear_mode_none) then
     if (use_old_diagnostics) then
       diff = diffusivity()
     else
#ifdef NEW_DIAG
       diff = gnostics%current_results%diffusivity
#endif 
     endif
     !write(*,*) 'difff is', diff
     do is = 1,nspec
       ! Q = n chi grad T = n (gamma / k^2) dT / dr
       ! = dens  n_r (gamma_N v_thr / k_N**2 rho_r a) dT / drho drho/dr
       ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) T a / L_T drho/dr
       ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) temp T_r tprim drho/dr_N/a
       ! Q / (n_r  v_r T_r rho_r**2/a**2) 
       ! = dens (gamma_N / k_N**2) temp tprim grho
       !   
       ! grho factored to diffusivity in diagnostics
       qf(is) =  diff * spec(is)%dens * spec(is)%temp * spec(is)%tprim 
       pf(is) =  diff * spec(is)%dens**2.0 * spec(is)%fprim 
       ht = 0.0
       vf = 0.0
     end do
     time_interval = 1.0
   else 
     if (use_old_diagnostics) then
       time_interval = user_time-start_time
       qf = qflux_avg
       pf = pflux_avg
       ht = heat_avg
       vf = vflux_avg
     else
#ifdef NEW_DIAG
       time_interval = user_time - gnostics%start_time
       qf = gnostics%current_results%species_heat_flux_avg
       pf = gnostics%current_results%species_particle_flux_avg
       ht = gnostics%current_results%species_heating_avg
       vf = gnostics%current_results%species_momentum_flux_avg
#endif
     endif
   end if

   state%outputs%pflux = pf/time_interval
   state%outputs%qflux = qf/time_interval
   state%outputs%heat = ht/time_interval
   state%outputs%vflux = vf(1)/time_interval
  end subroutine calculate_outputs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine reset_timers(timers)
    type(gs2_timers_type), intent(inout) :: timers
    timers%init = 0.
    timers%advance = 0.
    timers%timestep = 0.
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
       if (.not. state%print_full_timers) then
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
    if (.not. associated(state%outputs%pflux)) then
      allocate(state%outputs%pflux(nspec))
      allocate(state%outputs%qflux(nspec))
      allocate(state%outputs%heat(nspec))
    end if
  end subroutine allocate_outputs

  subroutine deallocate_outputs(state)
    type(gs2_program_state_type), intent(inout) :: state
    if (associated(state%outputs%pflux)) then
      deallocate(state%outputs%pflux)
      deallocate(state%outputs%qflux)
      deallocate(state%outputs%heat)
    end if
  end subroutine deallocate_outputs







  !> Deprecated. This is the main subroutine in which gs2 is initialized, equations are advanced,
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

   
   use job_manage, only: trin_restart, trin_reset, time_message
    use unit_tests, only: functional_test_flag, ilast_step
    use mp, only: proc0
    use old_interface_store, only: override_profiles, override_miller_geometry
    implicit none

    integer, intent (in), optional :: mpi_comm, job_id, nensembles
    character (*), intent (in), optional :: filename
    real, dimension (:), intent (out), optional :: pflux, qflux, heat
    real, intent (out), optional :: vflux
    real, intent (out), optional :: dvdrho, grho
    logical, intent (in), optional :: trinity_reset
    logical, intent (out), optional :: converged

    logical, save :: first_time = .true.
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
      if (override_profiles .or. override_miller_geometry) then 
        if (proc0) call time_message(.false., old_iface_state%timers%init,' Initialization')
        call old_interface_set_overrides
        if (proc0) call time_message(.false., old_iface_state%timers%init,' Initialization')
      end if
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
      !write (*,*) 'SETTING FLUXES'
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
    call finalize_overrides(old_iface_state)

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
    use normalisations, only: finish_normalisations
    implicit none

    call finish_normalisations
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
  
  !> This function calls reset_gs2 using
  !! all the default parameters, multiplied
  !! by a factor. It is used
  !! in the gs2_reinit unit test
  function gs2_main_unit_test_reset_gs2(fac)
    use species, only: spec, nspec
    use dist_fn, only: g_exb 
    use mp, only: mp_abort
    logical :: gs2_main_unit_test_reset_gs2
    real, intent(in) :: fac

    if (nspec.ne.2)  & 
      call mp_abort("gs2_main_unit_test_reset_gs2 only works with 2 species", .true.)

    !write (*,*) 'HEREEEE'

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
    !write (*,*) 'HERAAAA'

    gs2_main_unit_test_reset_gs2 = .true.

  end function gs2_main_unit_test_reset_gs2

  subroutine reset_gs2 (ntspec, dens, temp, fprim, tprim, gexb, mach, nu, nensembles)

    use dist_fn, only: dist_fn_g_exb => g_exb
    use gs2_init, only: init, init_level_list
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use gs2_diagnostics, only: gd_reset => reset_init
    use gs2_save, only: gs2_save_for_restart!, gs_reset => reset_init
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt, save_dt
    use mp, only: scope, subprocs, allprocs
    use mp, only: mp_abort
    use run_parameters, only: trinity_linear_fluxes
    use species, only: nspec, spec, impurity, ions, electrons
    use species, only: determine_species_order
    use run_parameters, only: use_old_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: reset_averages_and_counters
#endif

    implicit none

    integer, intent (in) :: ntspec, nensembles
    real, intent (in) :: gexb, mach
    real, dimension (:), intent (in) :: dens, fprim, temp, tprim, nu
    integer :: is
    integer :: isg
    logical :: temp_initval_override_store

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


    if (nensembles > 1) call scope (subprocs)

    if (trinity_linear_fluxes.and.nonlinear_mode_switch.eq.nonlinear_mode_none) &
      call reset_linear_magnitude


    call prepare_initial_values_overrides(old_iface_state)
    call set_initval_overrides_to_current_vals(old_iface_state%init%initval_ov)
    temp_initval_override_store = old_iface_state%init%initval_ov%override 
    old_iface_state%init%initval_ov%override = .true.

    ! EGH is this line necessary?
    gnew = 0.
    call save_dt (code_dt)
    ! This call to gs2_diagnostics::reset_init sets some time averages
    ! and counters to zero... used mainly for trinity convergence checks.
    if (use_old_diagnostics) then
      call gd_reset
    else
#ifdef NEW_DIAG
      call reset_averages_and_counters
#endif
    end if

    call prepare_profiles_overrides(old_iface_state)


    call determine_gs2spec_from_trin(ntspec)
    do is = 1,nspec
      isg = gs2spec_from_trin(is)
      old_iface_state%init%prof_ov%override_dens(isg) = .true.
      old_iface_state%init%prof_ov%dens(isg) = dens(is)/dens(densrefspec)

      old_iface_state%init%prof_ov%override_temp(isg) = .true.
      old_iface_state%init%prof_ov%temp(isg) = temp(is)/temp(1)

      old_iface_state%init%prof_ov%override_fprim(isg) = .true.
      old_iface_state%init%prof_ov%fprim(isg) = fprim(is)

      old_iface_state%init%prof_ov%override_tprim(isg) = .true.
      old_iface_state%init%prof_ov%tprim(isg) = tprim(is)

      old_iface_state%init%prof_ov%override_vnewk(isg) = .true.
      old_iface_state%init%prof_ov%vnewk(isg) = nu(is)

      !call override(old_iface_state, odens, gs2spec_from_trin(is), dens(is)/dens(densrefspec))
      !! Temp is always normalised to the main ions
      !call override(old_iface_state, otemp, gs2spec_from_trin(is), temp(is)/temp(1))
      !call override(old_iface_state, ofprim, gs2spec_from_trin(is), fprim(is))
      !call override(old_iface_state, otprim, gs2spec_from_trin(is), tprim(is))
      !call override(old_iface_state, ovnewk, gs2spec_from_trin(is), nu(is))
    end do

    !old_iface_state%init%prof_ov%override_g_exb = .true.
    !old_iface_state%init%prof_ov%g_exb = gexb
    !old_iface_state%init%prof_ov%override_mach = .true.
    !old_iface_state%init%prof_ov%mach = mach

    !call override(old_iface_state, og_exb, gexb)
    ! mach doesn't do anything atm
    !call override(old_iface_state, omach, mach)

    call init(old_iface_state%init, init_level_list%full)
    old_iface_state%istep_end = 1
    old_iface_state%init%initval_ov%override = temp_initval_override_store
    if (nensembles > 1) call scope (allprocs)


  end subroutine reset_gs2

  subroutine old_interface_set_overrides
    use old_interface_store
    use gs2_init, only: init, init_level_list
    integer :: is, isg

    if (override_miller_geometry) then 
      call prepare_miller_geometry_overrides(old_iface_state)
      old_iface_state%init%mgeo_ov%override_rhoc = .true.
      old_iface_state%init%mgeo_ov%rhoc = rhoc_store
      old_iface_state%init%mgeo_ov%override_qinp = .true.
      old_iface_state%init%mgeo_ov%qinp = qval_store
      old_iface_state%init%mgeo_ov%override_shat = .true.
      old_iface_state%init%mgeo_ov%shat = shat_store
      old_iface_state%init%mgeo_ov%override_rgeo_lcfs = .true.
      old_iface_state%init%mgeo_ov%rgeo_lcfs = rgeo_lcfs_store
      old_iface_state%init%mgeo_ov%override_rgeo_local = .true.
      old_iface_state%init%mgeo_ov%rgeo_local = rgeo_local_store
      old_iface_state%init%mgeo_ov%override_akappa = .true.
      old_iface_state%init%mgeo_ov%akappa = kap_store
      old_iface_state%init%mgeo_ov%override_akappri = .true.
      old_iface_state%init%mgeo_ov%akappri = kappri_store
      old_iface_state%init%mgeo_ov%override_tri = .true.
      old_iface_state%init%mgeo_ov%tri = tri_store
      old_iface_state%init%mgeo_ov%override_tripri = .true.
      old_iface_state%init%mgeo_ov%tripri = tripri_store
      old_iface_state%init%mgeo_ov%override_shift = .true.
      old_iface_state%init%mgeo_ov%shift = shift_store
      old_iface_state%init%mgeo_ov%override_betaprim = .true.
      old_iface_state%init%mgeo_ov%betaprim = betaprim_store
    end if


    if (override_profiles) then
      call prepare_profiles_overrides(old_iface_state)
      call determine_gs2spec_from_trin(ntspec_store)
      do is = 1,ntspec_store
        isg = gs2spec_from_trin(is)
        old_iface_state%init%prof_ov%override_dens(isg) = .true.
        old_iface_state%init%prof_ov%dens(isg) = dens_store(is)/dens_store(densrefspec)

        old_iface_state%init%prof_ov%override_temp(isg) = .true.
        old_iface_state%init%prof_ov%temp(isg) = temp_store(is)/temp_store(1)

        old_iface_state%init%prof_ov%override_fprim(isg) = .true.
        old_iface_state%init%prof_ov%fprim(isg) = fprim_store(is)

        old_iface_state%init%prof_ov%override_tprim(isg) = .true.
        old_iface_state%init%prof_ov%tprim(isg) = tprim_store(is)

        old_iface_state%init%prof_ov%override_vnewk(isg) = .true.
        old_iface_state%init%prof_ov%vnewk(isg) = nu_store(is)
      end do
    end if
    override_profiles = .false.
    override_miller_geometry = .false.
  end subroutine old_interface_set_overrides

  subroutine determine_gs2spec_from_trin(ntspec)
    use species, only: determine_species_order, spec, nspec, impurity, ions
    use species, only: electrons
    use mp, only: mp_abort
    integer, intent(in) :: ntspec
    call determine_species_order
    if (.not. allocated(gs2spec_from_trin)) allocate(gs2spec_from_trin(ntspec))
    !write (*,*) 'IONS', ions, 'ELECTRONS', electrons, 'IMPURITY', impurity
    if (nspec==1) then
      gs2spec_from_trin(1) = ions
      ! ref temp is always ions (trin spec 1)
      ! For one species, reference dens is ions
      densrefspec = 1
    else if (nspec==2) then 
      gs2spec_from_trin(1) = ions
      gs2spec_from_trin(2) = electrons
      ! For two species or more, ref dens is electrons (trin spec 2)
      densrefspec = 2
    else if (nspec==3) then
      gs2spec_from_trin(1) = ions
      gs2spec_from_trin(2) = electrons
      gs2spec_from_trin(3) = impurity
      densrefspec = 2
    else
      call mp_abort("Can't handle more than 3 species in &
       & determine_gs2spec_from_trin", .true.)
    end if
  end subroutine determine_gs2spec_from_trin

  subroutine reset_linear_magnitude
    use dist_fn_arrays, only: g, gnew
    use fields, only: set_init_fields
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
     call set_init_fields
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

    use dist_fn, only: dist_fn_g_exb => g_exb
    use old_interface_store
    !use mp, only: broadcast

    implicit none

    integer, intent (inout) :: ntspec
    real, intent (inout) :: rhoc, qval, shat, rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift
    real, intent (inout) :: betaprim, gexb, mach
    real, dimension (:), intent (inout) :: dens, fprim, temp, tprim, nu
    logical, intent (inout) :: use_gs2_geo

    ! for now do nothing with gexb or mach, but need to include later if want to use GS2
    ! with TRINITY to evolve rotation profiles
     if (.not. allocated(dens_store)) then
       allocate(dens_store(ntspec))
       allocate(temp_store(ntspec))
       allocate(tprim_store(ntspec))
       allocate(fprim_store(ntspec))
       allocate(nu_store(ntspec))
     end if
     ntspec_store = ntspec
     dens_store = dens
     temp_store = temp
     tprim_store = tprim
     fprim_store = fprim
     nu_store = nu
     override_profiles = .true.
     if (.not. use_gs2_geo) then
       rhoc_store = rhoc
       qval_store = qval
       shat_store = shat
       rgeo_lcfs_store = rgeo_lcfs
       rgeo_local_store = rgeo_local
       kap_store = kap
       kappri_store = kappri
       tri_store = tri
       tripri_store = tripri
       shift_store = shift
       betaprim_store = betaprim
       override_miller_geometry = .true.
     end if
     
    
  end subroutine gs2_trin_init

# ifndef MAKE_LIB
end module gs2_main
# endif


