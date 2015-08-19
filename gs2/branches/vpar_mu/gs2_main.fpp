# ifndef MAKE_LIB

!> This module contains the functions gs2_main::run_gs2, gs2_main::finish_gs2 and gs2_main::reset_gs2, whose
!! purpose should be reasonably obvious. These functions were originally part of
!! program GS2, but have now been moved to this module so that GS2 itself can be
!! called as a library by, for example, Trinity. All the program GS2 does is
!! include this module and call run_gs2.

module gs2_main
  
  public :: run_gs2, finish_gs2, reset_gs2
  
contains
# endif


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
     pflux, qflux, vflux, heat, dvdrho, grho)

    use job_manage, only: checkstop, job_fork, checktime, time_message
    use mp, only: init_mp, finish_mp, proc0, nproc, broadcast, scope, subprocs
    use mp, only: max_reduce, min_reduce, sum_reduce
    use file_utils, only: init_file_utils, run_name
    use fields, only: init_fields, advance
    use species, only: ions, electrons, impurity
    use gs2_diagnostics, only: init_gs2_diagnostics, finish_gs2_diagnostics
    use parameter_scan, only: init_parameter_scan, allocate_target_arrays
    use gs2_diagnostics, only: nsave, pflux_avg, qflux_avg, heat_avg, vflux_avg, start_time
    use run_parameters, only: nstep, fphi, fapar, fbpar, avail_cpu_time
    use dist_fn_arrays, only: gnew
    use gs2_save, only: gs2_save_for_restart
    use gs2_diagnostics, only: loop_diagnostics, ensemble_average
    use gs2_reinit, only: reset_time_step, check_time_step, time_reinit
    use gs2_time, only: update_time, write_dt, init_tstart
    use gs2_time, only: user_time, user_dt
    use init_g, only: tstart
    use geometry, only: surfarea, dvdrhon
    use redistribute, only: time_redist
    use fields_implicit, only: time_field
    use gs2_layouts, only: layout
    use parameter_scan, only: update_scan_parameter_value
    implicit none

    integer, intent (in), optional :: mpi_comm, job_id, nensembles
    character (*), intent (in), optional :: filename
    real, dimension (:), intent (out), optional :: pflux, qflux, heat
    real, intent (out), optional :: vflux
    real, intent (out), optional :: dvdrho, grho

    real :: time_init(2) = 0., time_advance(2) = 0., time_finish(2) = 0.
    real :: time_total(2) = 0.
    real :: time_interval
    real :: time_main_loop(2)

!    real :: t1, t2, t3, t4

    integer :: istep = 0, istatus, istep_end
    logical :: exit, reset, list
    logical :: first_time = .true.
    logical :: nofin= .false.
!    logical, optional, intent(in) :: nofinish
    character (500), target :: cbuff

    time_main_loop(1) = 0.
    time_main_loop(2) = 0.

!
!CMR, 12/2/2010: 
!     add nofinish optional variable to avoid deallocations at end of simulation
!     as may want to do post-processing
!
!    if (present(nofinish)) nofin=nofinish
    
    if (first_time) then

       ! <doc> Initialize message passing </doc>
       if (present(mpi_comm)) then
          call init_mp (mpi_comm)
       else
          call init_mp
       end if
       call checktime(avail_cpu_time,exit) ! <doc> Initialize timer </doc>
       
       ! <doc> Report # of processors being used </doc>
       if (proc0) then
          if (nproc == 1) then
             if (.not. nofin) then
                write(*,*) 'Running on ',nproc,' processor'
             end if
          else
             if (.not. nofin) then
                write(*,*) 'Running on ',nproc,' processors'
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
       end if
       
       call broadcast (list)
       
       ! <doc> If given a list of jobs, fork </doc>
       if (list) then
          call job_fork
       else if (present(nensembles)) then
          if (nensembles > 1) call job_fork (n_ensembles=nensembles)
       end if
       if (proc0) call time_message(.false.,time_total,' Total')

       if (proc0) then
          call time_message(.false., time_init,' Initialization')
          cbuff = trim(run_name)
       end if
       
       call broadcast (cbuff)
       if (.not. proc0) run_name => cbuff

       call init_parameter_scan

       call init_fields

       call init_gs2_diagnostics (list, nstep)
       call allocate_target_arrays ! must be after init_gs2_diagnostics
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
       
       first_time = .false.

    else if (present(nensembles)) then
       if (nensembles > 1) then
          call scope (subprocs)
       end if
    end if
    
    istep_end = nstep
    
    call loop_diagnostics(0,exit)
    
    if (proc0) write(*,*) 'layout ',layout

    call time_message(.false.,time_main_loop,' Main Loop')

    do istep = 1, nstep

       if (proc0) call time_message(.false.,time_advance,' Advance time step')
       call advance (istep)

       if (nsave > 0 .and. mod(istep, nsave) == 0) &
            call gs2_save_for_restart (gnew, user_time, user_dt, istatus, fphi, fapar, fbpar)
       call update_time
       call loop_diagnostics (istep, exit)
       call check_time_step (reset, exit)
!       call update_scan_parameter_value(istep, reset, exit)
       if (proc0) call time_message(.false.,time_advance,' Advance time step')
       if (reset) then
          ! if called within trinity, do not dump info to screen
          if (present(job_id)) then
             call reset_time_step (istep, exit, job_id)
          else       
             call reset_time_step (istep, exit)
          end if
       end if

       if (mod(istep,5) == 0) call checkstop(exit)
       
       call checktime(avail_cpu_time,exit)

       if (exit) then
          istep_end = istep
          exit
       end if
    end do

    call time_message(.false.,time_main_loop,' Main Loop')

    if (proc0) call time_message(.false.,time_finish,' Finished run')

    if (proc0 .and. .not. present(job_id)) call write_dt

    time_interval = user_time-start_time

    if (present(nensembles)) then
       if (nensembles > 1) call ensemble_average (nensembles, time_interval)
    end if

    if (present(pflux)) then
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
       if (.not.nofin ) call finish_gs2_diagnostics
       if (.not.nofin) call finish_gs2
    end if
    
    if (proc0) call time_message(.false.,time_finish,' Finished run')

    if (proc0) call time_message(.false.,time_total,' Total')

    if (proc0) then
       if (present(job_id)) then
          print '(/,'' Job ID:'', i4,'', total from timer is:'', 0pf9.2,'' min'',/)', &
               job_id+1, time_total(1)/60.
       else if (.not. nofin) then
!    if (proc0 .and. .not. nofin) then

          write (*,*)
          write (*,fmt=101) 'Initialization', time_init(1)/60., 'min', time_init(1)/time_total(1), '%'
          write (*,fmt=101) 'Advance steps', time_advance(1)/60., 'min', time_advance(1)/time_total(1), '%'
          write (*,fmt=101) '(redistribute)', time_redist(1)/60., 'min', time_redist(1)/time_total(1), '%'
          write (*,fmt=101) '(field solve)', time_field(1)/60., 'min', time_field(1)/time_total(1), '%'
          write (*,fmt=101) 'Re-initialize', time_reinit(1)/60., 'min', time_reinit(1)/time_total(1), '%'
          write (*,fmt=101) 'Finishing', time_finish(1)/60., 'min', time_finish(1)/time_total(1), '%'
          write (*,fmt=102) 'total from timer is:', time_total(1)/60., 'min'
       endif
    end if

101 format (a25,0pf8.2,a4,T40,2pf5.1,a2)
102 format (a25,0pf8.2,a4)
    
    if (.not. present(mpi_comm) .and. .not. nofin) call finish_mp
    
  end subroutine run_gs2
  
  subroutine finish_gs2
    
    use antenna, only: finish_antenna
    use dist_fn, only: finish_dist_fn
    use fields, only: finish_fields
    use file_utils, only: finish_file_utils
    use hyper, only: finish_hyper
    use init_g, only: finish_init_g
    use kt_grids, only: finish_kt_grids
    use vpamu_grids, only: finish_vpamu_grids
    use parameter_scan, only: finish_parameter_scan
    use mp, only: proc0
    use nonlinear_terms, only: finish_nonlinear_terms
    use run_parameters, only: finish_run_parameters
    use species, only: finish_species

    implicit none

    call finish_antenna
    call finish_dist_fn
    call finish_fields
    call finish_hyper
    call finish_init_g
    call finish_kt_grids
    call finish_vpamu_grids
    call finish_nonlinear_terms
    call finish_run_parameters
    call finish_species
    call finish_parameter_scan
    if (proc0) call finish_file_utils

  end subroutine finish_gs2

!  subroutine reset_gs2 (ntspec, dens, temp, fprim, tprim, gexb, mach, nu, nensembles)
  subroutine reset_gs2 (ntspec, dens, temp, fprim, tprim, nu, nensembles)

    use dist_fn, only: d_reset => reset_init
    use fields, only: init_fields, f_reset => reset_init
    use fields_implicit, only: fi_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use init_g, only: g_reset => reset_init
    use nonlinear_terms, only: nl_reset => reset_init
    use gs2_diagnostics, only: gd_reset => reset_init
    use gs2_save, only: gs2_save_for_restart!, gs_reset => reset_init
    use species, only: reinit_species
    use dist_fn_arrays, only: gnew
    use gs2_time, only: code_dt, user_dt, save_dt, user_time
    use run_parameters, only: fphi, fapar, fbpar
    use antenna, only: a_reset => reset_init
    use mp, only: scope, subprocs, allprocs

    implicit none

    integer, intent (in) :: ntspec, nensembles
!    real, intent (in) :: gexb, mach
    real, dimension (:), intent (in) :: dens, fprim, temp, tprim, nu

    integer :: istatus

    ! doing nothing with gexb or mach for now, but in future will need to when
    ! using GS2 to evolve rotation profiles in TRINITY

    if (nensembles > 1) call scope (subprocs)

    call gs2_save_for_restart (gnew, user_time, user_dt, istatus, fphi, fapar, fbpar)
    gnew = 0.

    call save_dt (code_dt)

    call d_reset
    call f_reset
    call fi_reset
    call ft_reset
!    if (.not. ql_flag) call g_reset
    call g_reset
    call nl_reset
    call gd_reset
    call a_reset
!!    call gs_reset

    call reinit_species (ntspec, dens, temp, fprim, tprim, nu)
    call init_fields

    if (nensembles > 1) call scope (allprocs)

  end subroutine reset_gs2

  subroutine gs2_trin_init (rhoc, qval, shat, rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift, &
!       betaprim, ntspec, dens, temp, fprim, tprim, gexb, mach, nu, use_gs2_geo)
       betaprim, ntspec, dens, temp, fprim, tprim, nu, use_gs2_geo)

    use species, only: init_trin_species
    use theta_grid_params, only: init_trin_geo

    implicit none

    integer, intent (in) :: ntspec
    real, intent (in) :: rhoc, qval, shat, rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift
!    real, intent (in) :: betaprim, gexb, mach
    real, intent (in) :: betaprim
    real, dimension (:), intent (in) :: dens, fprim, temp, tprim, nu
    logical, intent (in) :: use_gs2_geo

    ! for now do nothing with gexb or mach, but need to include later if want to use GS2
    ! with TRINITY to evolve rotation profiles

    call init_trin_species (ntspec, dens, temp, fprim, tprim, nu)
    if (.not. use_gs2_geo) call init_trin_geo (rhoc, qval, shat, &
         rgeo_lcfs, rgeo_local, kap, kappri, tri, tripri, shift, betaprim)
    
  end subroutine gs2_trin_init

# ifndef MAKE_LIB
end module gs2_main
# endif
