# ifndef MAKE_LIB

!> This module contains the functions gs2_main::run_gs2, gs2_main::finish_gs2 and gs2_main::reset_gs2, whose
!! purpose should be reasonably obvious. These functions were originally part of
!! program GS2, but have now been moved to this module so that GS2 itself can be
!! called as a library by, for example, Trinity. All the program GS2 does is
!! include this module and call run_gs2.

module gs2_main
  implicit none
  public :: run_gs2, finish_gs2, reset_gs2, trin_finish_gs2

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
    real :: dvdrho
    real :: grho
  end type gs2_outputs_type

  type gs2_program_state_type

    ! Flags indicating the current state of the 
    ! program (used for error checking)
    logical :: gs2_initialized = .false.
    logical :: equations_initialized = .false.
    logical :: diagnostics_initialized = .false.

    ! Timers
    type(gs2_timers_type) :: timers

    !> The exit flag is set to true by any 
    !! part of the main timestep loop that 
    !! wants to cause the loop to exit
    logical :: exit = .false.

    !> Whether to print out debug messages
    !logical :: debug = .false.
    integer :: verb = 3

    !> Parameters pertaining to Trinity runs
    !! external_job_id is not to confused with the parameter
    !! job in mp, which identifies the subjob if
    !! running in list mode or with nensembles > 1
    !integer :: trin_job = -1
    integer :: external_job_id = -1

    logical :: is_external_job = .false.


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
 
contains
# endif

  subroutine initialize_wall_clock_timer
    use job_manage, only: init_checktime
    call init_checktime
  end subroutine initialize_wall_clock_timer

  subroutine initialize_gs2(state)
    use file_utils, only: init_file_utils
    use file_utils, only: run_name, run_name_target
    use job_manage, only: checktime, time_message
    use job_manage, only: init_checktime, checktime_initialized
    use job_manage, only: job_fork
    use mp, only: init_mp, broadcast
    use mp, only: iproc, nproc, proc0
    use mp, only: trin_flag
    use redistribute, only: using_measure_scatter
    use run_parameters, only: avail_cpu_time
    use runtime_tests, only: verbosity
    use unit_tests, only: debug_message, set_job_id
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (state%gs2_initialized .or. &
        state%equations_initialized .or. &
        state%diagnostics_initialized) then
        write (*,*) 'ERROR: initialize_gs2 can only be called when &
        & gs2_initialized, equations_initialized and diagnostics_initialized, &
        & are all false. '
       stop 1
     end if

    if (state%is_external_job) then 
      call set_job_id(state%external_job_id)
    end if

    call debug_message(4, 'gs2_main::initialize_gs2 starting initialization')

    if (state%mp_comm_external) then
       call init_mp (state%mp_comm)
    else
       call init_mp
    end if
    ! We now manually set trin_flag, since we may now be passing 
    ! in a communicator but not running trinity.
    trin_flag = state%is_external_job
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

    state%gs2_initialized = .true.

  end subroutine initialize_gs2

  subroutine initialize_equations(state)
    use fields, only: init_fields
    use geometry, only: surfarea, dvdrhon
    use gs2_time, only: init_tstart
    use init_g, only: tstart
    use job_manage, only: time_message
    use mp, only: proc0, broadcast
    use parameter_scan, only: init_parameter_scan
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
    call init_fields

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

    call debug_message(state%verb, 'gs2_main::initialize_equations finished')

  end subroutine initialize_equations

  subroutine initialize_diagnostics
  end subroutine initialize_diagnostics

  subroutine evolve_equations(state)
    type(gs2_program_state_type), intent(inout) :: state

    ! Make sure exit is false before entering
    ! timestep loop
    state%exit = .false.

  end subroutine evolve_equations

  subroutine finalize_diagnostics
  end subroutine finalize_diagnostics

  subroutine finalize_equations(state)
    use gs2_layouts, only: finish_layouts
    use gs2_transforms, only: finish_transforms
    use gs2_save, only: gs2_save_for_restart, finish_save
    use theta_grid, only: finish_theta_grid
    use unit_tests, only: ilast_step
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
    use gs2_save, only: finish_save
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

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
    call finish_save
! HJL Species won't change during a run so shouldn't need this    
!    call finish_trin_species

    call finish_layouts
    call finish_transforms
    call finish_save
    call finish_theta_grid
  end subroutine finalize_equations

  subroutine finalize_gs2(state)
    use file_utils, only: finish_file_utils
    use job_manage, only: time_message
    use mp, only: finish_mp, proc0
    use mp, only: mp_abort
    use unit_tests, only: debug_message
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if ((.not. state%gs2_initialized) .or. &
        state%equations_initialized .or. &
        state%diagnostics_initialized) then
        write (*,*) 'ERROR: initialize_gs2 can only be called when &
        & gs2_initialized is true, and equations_initialized &
        & and diagnostics_initialized, &
        & are all false. '
       stop 1
     end if

    if (proc0) call finish_file_utils

    if (proc0) call time_message(.false.,state%timers%finish,' Finished run')

    if (proc0) call time_message(.false.,state%timers%total,' Total')
    
    call debug_message(state%verb, 'gs2_main::finalize_gs2 calling print_times')

    call print_times(state, state%timers)

    call debug_message(state%verb, 'gs2_main::finalize_gs2 calling finish_mp')

    if (.not. state%mp_comm_external) call finish_mp

    state%gs2_initialized = .false.



    !if (.not. present(mpi_comm) .and. .not. nofin) call finish_mp
  end subroutine finalize_gs2

  subroutine override_geometry
  end subroutine override_geometry

  subroutine override_profiles
  end subroutine override_profiles


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

    use job_manage, only: checkstop, job_fork, checktime, time_message, trin_reset, trin_restart, trin_job
    use mp, only: init_mp, finish_mp, proc0, nproc, broadcast, scope, subprocs
    use mp, only: max_reduce, min_reduce, sum_reduce, mp_abort
    use file_utils, only: init_file_utils, run_name!, finish_file_utils
    use fields, only: init_fields, advance
    use species, only: ions, electrons, impurity, spec, nspec
    use gs2_diagnostics, only: init_gs2_diagnostics, finish_gs2_diagnostics
    use parameter_scan, only: init_parameter_scan, allocate_target_arrays
    use gs2_diagnostics, only: nsave, pflux_avg, qflux_avg, heat_avg, vflux_avg, start_time, nwrite, write_nl_flux
    use run_parameters, only: nstep, fphi, fapar, fbpar, avail_cpu_time, margin_cpu_time
    use run_parameters, only: trinity_linear_fluxes, do_eigsolve, reset
    use dist_fn_arrays, only: gnew
    use gs2_save, only: gs2_save_for_restart
    use gs2_diagnostics, only: loop_diagnostics, ensemble_average
    use gs2_diagnostics, only: diffusivity
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: run_diagnostics, init_gs2_diagnostics_new
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new, diagnostics_init_options_type
#endif
    use gs2_reinit, only: reset_time_step, check_time_step, time_reinit
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use gs2_time, only: update_time, write_dt, init_tstart
    use gs2_time, only: user_time, user_dt
    use init_g, only: tstart
    use collisions, only: vnmult
    use geometry, only: surfarea, dvdrhon
    use redistribute, only: time_redist, using_measure_scatter
    use fields_arrays, only: time_field
    use parameter_scan, only: update_scan_parameter_value
    use unit_tests, only: functional_test_flag, ilast_step
    use unit_tests,only: should_print
    use unit_tests, only: set_job_id
    use runtime_tests, only: verbosity
    use mp, only: broadcast, iproc
#ifdef WITH_EIG
    use eigval, only: init_eigval, finish_eigval, time_eigval
    use eigval, only: BasicSolve
#endif
    implicit none

    integer, intent (in), optional :: mpi_comm, job_id, nensembles
    character (*), intent (in), optional :: filename
    real, dimension (:), intent (out), optional :: pflux, qflux, heat
    real, intent (out), optional :: vflux
    real, intent (out), optional :: dvdrho, grho
    logical, intent (in), optional :: trinity_reset
    logical, intent (out), optional :: converged

    real :: time_init(2) = 0., time_advance(2) = 0., time_finish(2) = 0.
    real :: time_total(2) = 0., time_diagnostics(2)=0.
    real :: time_interval
    real :: time_main_loop(2)
#ifdef NEW_DIAG
    real :: precision_test
#endif
    real :: diff
    integer :: istep = 0, istatus, istep_end
    integer :: is
    logical :: exit, list
    logical :: first_time = .true.
    logical :: nofin= .false.
    logical :: debug
!    logical, optional, intent(in) :: nofinish
    character (500), target :: cbuff

#ifdef NEW_DIAG
    type(diagnostics_init_options_type) :: diagnostics_init_options
#endif

    time_main_loop(1) = 0.
    time_main_loop(2) = 0.
    exit=.false.

    debug = (verbosity() .gt. 2)


!
!CMR, 12/2/2010: 
!     add nofinish optional variable to avoid deallocations at end of simulation
!     as may want to do post-processing
!
!    if (present(nofinish)) nofin=nofinish
     

! HJL tests on Trinity optionals for load balancing
    trin_job = -1
    if(present(job_id)) trin_job = job_id + 1
    if(present(trinity_reset)) then
       converged = .false.
       trin_reset = trinity_reset
       first_time = .true.
       if (debug) write (*,*) 'gs2_main::run_gs2 set first_time = .true., jid=', trin_job
    endif

#ifdef NEW_DIAG
    diagnostics_init_options%initialized = (.not. first_time)
#endif

    if (.not. first_time) then
      call broadcast(trin_job)
      call set_job_id(trin_job)
      debug = should_print(3)
    end if
    if (first_time) then
      if (debug) write (*,*) 'gs2_main::run_gs2 starting initialization', iproc, trin_job

       ! <doc> Initialize message passing </doc>
       if (present(mpi_comm)) then
          call init_mp (mpi_comm)
       else
          call init_mp
       end if
       call broadcast(trin_job)
       call set_job_id(trin_job)
       debug = should_print(3)
       write (*,*) 'gs2_main::run_gs2 initialized mp, nproc= ', nproc, iproc, trin_job
       call checktime(avail_cpu_time,exit) ! <doc> Initialize timer </doc>
       write (*,*) 'gs2_main::run_gs2 called checktime, avail_cpu_time = ', avail_cpu_time, iproc, trin_job
       
       ! <doc> Report # of processors being used </doc>
       if (proc0) then
          if (nproc == 1) then
             if (.not. nofin) then
                write(*,*) 'Running on ',nproc,' processor'
             end if
          else
             if (.not. nofin) then
                if(present(job_id)) then
                   write(*,*) 'Job ',trin_job,'Running on ',nproc,' processors'
                else
                   write(*,*) 'Running on ',nproc,' processors'
                endif
             end if
          end if

          write (*,*) 
          ! <doc> Call init_file_utils, ie. initialize the inputs and outputs, checking 
          !  whether we are doing a [[Trinity]] run or a list of runs. </doc>
          ! <doc>If it is a [[Trinity]] run then [[filename]] (the name of the input file?) is passed to  init_file_utils</doc>
          ! <doc> Figure out run name or get list of jobs </doc>
          if (present(filename)) then
             call init_file_utils (list, trin_run=.true., name=filename, n_ensembles=nensembles)
          else
             call init_file_utils (list, name="gs")
          end if
          if (debug) write (*,*) 'gs2_main::run_gs2 initialized file_utils', iproc, trin_job
       end if
       
       call broadcast (list)
       if (debug) write (*,*) 'gs2_main::run_gs2 broadcasted list', iproc, trin_job
       
       ! <doc> If given a list of jobs, fork </doc>
       if (list) then
          call job_fork
          if (debug) write (*,*) 'gs2_main::run_gs2 called job fork', trin_job
       else if (present(nensembles)) then
          if (nensembles > 1) call job_fork (n_ensembles=nensembles)
       end if
       if (proc0) call time_message(.false.,time_total,' Total')

       if (proc0) then
          call time_message(.false., time_init,' Initialization')
          cbuff = trim(run_name)
       end if
       
       !write (*,*) 'gs2_main::run_gs2 before broadcast, run_name = ', cbuff, iproc, trin_job
       call broadcast (cbuff)
       if (.not. proc0) run_name => cbuff
       if (debug) write (*,*) 'gs2_main::run_gs2 run_name = ', run_name, iproc, trin_job
       call init_parameter_scan
       !Set using_measure_scatter to indicate we want to use in "gather/scatter" timings
       using_measure_scatter=.false.
       if (debug) write (*,*) 'gs2_main::run_gs2 calling init_fields', trin_job
       call init_fields
       if (debug) write (*,*) 'gs2_main::run_gs2 calling init_gs2_diagnostics', trin_job
       call init_gs2_diagnostics (list, nstep)

#ifdef NEW_DIAG
#ifdef NETCDF_PARALLEL
       diagnostics_init_options%parallel_io_capable = .true.
#else
       diagnostics_init_options%parallel_io_capable = .false.
#endif

       ! Here we check if reals have been promoted to doubles
       diagnostics_init_options%default_double =  (precision(precision_test).gt.10)
       ! Check whether this is a Trinity run... enforces calculation of the
       ! fluxes
       diagnostics_init_options%is_trinity_run = present(mpi_comm)

       !if (first_time) diagnostics_init_options%initialized = .false.
       !write (*,*) 'diagnostics_init_options%initialized', diagnostics_init_options%initialized
       if (.not. diagnostics_init_options%initialized) then
         call init_gs2_diagnostics_new(diagnostics_init_options)
       end if
#endif

       call allocate_target_arrays(nwrite,write_nl_flux) ! must be after init_gs2_diagnostics
       call init_tstart (tstart)   ! tstart is in user units 
       
       if (present(dvdrho)) then
          if (proc0) then
             dvdrho = dvdrhon
             grho = surfarea/dvdrhon
          end if
          call broadcast (dvdrho)
          call broadcast (grho)
       end if
       
       if (proc0) call time_message(.false.,time_init,' Initialization')
       
#ifdef NEW_DIAG
       ! Create variables
       if (.not. diagnostics_init_options%initialized) then
         call run_diagnostics(-1,exit)
         !diagnostics_init_options%initialized = .true.
       end if
       ! Write initial values
       call run_diagnostics(0,exit)
#endif
       first_time = .false.

       ! When trinity starts a new step it needs to reset after initialisation
       if(present(trinity_reset) .and. trin_reset) return

    else if (present(nensembles)) then
       if (nensembles > 1) then
          call scope (subprocs)
       end if
    endif !firstime

    istep_end = nstep
    ilast_step = nstep
    
    call loop_diagnostics(0,exit)

    if(present(trinity_reset)) trin_restart = .true. ! All trinity runs are restarted except the first


    if (present(pflux)) call write_trinity_parameters

    call time_message(.false.,time_main_loop,' Main Loop')
if(do_eigsolve)then
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
else
    do istep = 1, nstep

       if (proc0) call time_message(.false.,time_advance,' Advance time step')

       !Initialise reset to true
       reset=.true.

       do while(reset)
          reset=.false. !So that we only do this once unless something triggers a reset
          call advance (istep)

          !If we've triggered a reset then actually reset
          if (reset) then
             ! if called within trinity, do not dump info to screen
             if (present(job_id)) then
                call reset_time_step (istep, exit, job_id)
             else       
                call reset_time_step (istep, exit)
             end if
          end if
          if(exit) exit
       enddo
       
       if (nsave > 0 .and. mod(istep, nsave) == 0) &
            call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
       call update_time
       if(proc0) call time_message(.false.,time_diagnostics,' Diagnostics')
       call loop_diagnostics (istep, exit)
       if(exit .and. present(converged)) converged = .true.
#ifdef NEW_DIAG
       call run_diagnostics (istep, exit)
#endif
       if(proc0) call time_message(.false.,time_diagnostics,' Diagnostics')
       if (proc0) call time_message(.false.,time_advance,' Advance time step')

       if(.not.exit)then

          !Note this should only trigger a reset for timesteps too small
          !as timesteps too large are already handled
          call check_time_step(reset,exit)
          call update_scan_parameter_value(istep, reset, exit)

          !If something has triggered a reset then reset here
          if (reset) then
             ! if called within trinity, do not dump info to screen
             if (present(job_id)) then
                call reset_time_step (istep, exit, job_id)
             else       
                call reset_time_step (istep, exit)
             end if
          end if

          if ((mod(istep,5) == 0).and.(.not.exit)) call checkstop(exit)
          if (.not.exit) call checktime(avail_cpu_time,exit,margin_cpu_time)
       endif

       if (exit) then
          istep_end = istep
          ilast_step = istep
          exit
       end if
    end do
endif
    call time_message(.false.,time_main_loop,' Main Loop')

    if (proc0) call time_message(.false.,time_finish,' Finished run')

    if (proc0 .and. .not. present(job_id)) call write_dt

    time_interval = user_time-start_time

    if (present(nensembles)) then
       if (nensembles > 1) call ensemble_average (nensembles, time_interval)
    end if

    if (present(pflux)) then
      
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
         if (present(converged)) then
         end if
       end if
       if (size(pflux) > 1) then
          pflux(1) = pflux_avg(ions)/time_interval
          qflux(1) = qflux_avg(ions)/time_interval
          heat(1) = heat_avg(ions)/time_interval
          pflux(2) = pflux_avg(electrons)/time_interval
          qflux(2) = qflux_avg(electrons)/time_interval
          heat(2) = heat_avg(electrons)/time_interval
          if (size(pflux) > 2) then
             pflux(3) = pflux_avg(impurity)/time_interval
             qflux(3) = qflux_avg(impurity)/time_interval
             heat(3) = heat_avg(impurity)/time_interval
          end if
       else
          pflux = pflux_avg/time_interval
          qflux = qflux_avg/time_interval
          heat = heat_avg/time_interval
       end if
       vflux = vflux_avg(1)/time_interval
    else
#ifdef NEW_DIAG
       if (.not.nofin .and. .not. functional_test_flag ) call finish_gs2_diagnostics_new 
#endif
       if (.not.nofin .and. .not. functional_test_flag ) call finish_gs2_diagnostics (istep_end)
       if (.not.nofin .and. .not. functional_test_flag) call finish_gs2
    end if

    if (proc0) call time_message(.false.,time_finish,' Finished run')

    if (proc0) call time_message(.false.,time_total,' Total')

    if (proc0) then
       if (present(job_id)) then
          print '(/,'' Job ID:'', i4,'', total from timer is:'', 0pf9.2,'' min'',/)', &
               job_id+1, time_total(1)/60.
       else if (.not. nofin) then
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
               time_init(1)/60.,time_init(1)/time_total(1), &
#ifdef WITH_EIG
               time_eigval(1)/60.,time_eigval(1)/time_total(1), &
#endif
               time_advance(1)/60.,time_advance(1)/time_total(1), &
               time_redist(1)/60.,time_redist(1)/time_total(1), &
               time_field(1)/60.,time_field(1)/time_total(1), &
               time_diagnostics(1)/60.,time_diagnostics(1)/time_total(1), &
               time_reinit(1)/60.,time_reinit(1)/time_total(1), &
               time_finish(1)/60.,time_finish(1)/time_total(1),time_total(1)/60.
       endif
    end if

    if (.not. present(mpi_comm) .and. .not. nofin) call finish_mp
    
  end subroutine run_gs2
  
! HJL <
  subroutine trin_finish_gs2

    use gs2_layouts, only: finish_layouts
    use gs2_transforms, only: finish_transforms
    use gs2_diagnostics, only: finish_gs2_diagnostics
    use gs2_save, only: gs2_save_for_restart, finish_save
    use theta_grid, only: finish_theta_grid
    use unit_tests, only: ilast_step
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
#endif
    
    call finish_gs2_diagnostics (ilast_step)
#ifdef NEW_DIAG
    call finish_gs2_diagnostics_new
#endif
    call finish_gs2
! HJL Species won't change during a run so shouldn't need this    
!    call finish_trin_species

    call finish_layouts
    call finish_transforms
    call finish_save
    call finish_theta_grid

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
    use gs2_save, only: finish_save
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
    call finish_save
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

    call reset_gs2(nspec, &
      (/spec(1)%dens, spec(2)%dens/)*fac, &
      ! Deliberately leave fac off this line
      (/spec(1)%temp, spec(2)%temp/), &
      (/spec(1)%fprim, spec(2)%fprim/)*fac, &
      (/spec(1)%tprim, spec(2)%tprim/)*fac, &
      g_exb, 0.0, &
      (/spec(1)%vnewk, spec(2)%vnewk/), &
      1)

    gs2_main_unit_test_reset_gs2 = .true.

  end function gs2_main_unit_test_reset_gs2

  subroutine reset_gs2 (ntspec, dens, temp, fprim, tprim, gexb, mach, nu, nensembles)

    use dist_fn, only: dist_fn_g_exb => g_exb
    use gs2_reinit, only: save_fields_and_dist_fn
    use gs2_reinit, only: reinit_gk_and_field_equations
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

    implicit none

    integer, intent (in) :: ntspec, nensembles
    real, intent (in) :: gexb, mach
    real, dimension (:), intent (in) :: dens, fprim, temp, tprim, nu

    real :: dummy

    !integer :: istatus

    ! doing nothing with gexb or mach for now, but in future will need to when
    ! using GS2 to evolve rotation profiles in TRINITY
    ! EGH add a check for this.
    if (gexb .ne. dist_fn_g_exb) call mp_abort(&
      "ERROR: Changing g_exb in gs2_reset is not implemented yet.", .true.)

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


    call reinit_species (ntspec, dens, temp, fprim, tprim, nu)

    ! Antenna should be reset because it depends on species paramters
    call reinit_gk_and_field_equations(reset_antenna=.true.)

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
