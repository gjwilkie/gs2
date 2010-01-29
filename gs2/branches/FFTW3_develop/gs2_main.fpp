# ifndef MAKE_LIB
module gs2_main
  
  public :: run_gs2, finish_gs2, reset_gs2
  
contains
# endif

  subroutine run_gs2 (mpi_comm, filename, nensembles, pflux, qflux, heat, dvdrho, grho)

    ! <doc> Main subroutine in which gs2 is initialized, equations are advanced,
    ! and the program is finalized </doc>

    use job_manage, only: checkstop, job_fork, checktime
    use mp, only: init_mp, finish_mp, proc0, nproc, broadcast, scope, subprocs
    use file_utils, only: init_file_utils, run_name, list_name!, finish_file_utils
    use fields, only: init_fields
    use gs2_diagnostics, only: init_gs2_diagnostics, finish_gs2_diagnostics
    use gs2_diagnostics, only: nsave, pflux_avg, qflux_avg, heat_avg, start_time
    use run_parameters, only: nstep
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: avail_cpu_time
    use fields, only: advance
    use dist_fn_arrays, only: gnew
    use gs2_save, only: gs2_save_for_restart
    use gs2_diagnostics, only: loop_diagnostics, ensemble_average
    use gs2_reinit, only: reset_time_step, check_time_step
    use gs2_reinit, only: time_message, time_nc, time_reinit
    use gs2_time, only: update_time
    use gs2_time, only: write_dt, init_tstart
    use gs2_time, only: user_time, user_dt
    use gs2_time, only: code_time
    use init_g, only: tstart
    use collisions, only: vnmult
    use geometry, only: surfarea, dvdrhon

    implicit none

    integer, intent (in), optional :: mpi_comm, nensembles
    character (*), intent (in), optional :: filename
    real, dimension (:), intent (out), optional :: pflux, qflux, heat
    real, intent (out), optional :: dvdrho, grho

    real :: time_init = 0., time_advance = 0., time_finish = 0., time_total
    real :: time_interval
    integer :: istep = 0, istatus, istep_end
    logical :: exit, reset, list
    logical :: first_time = .true.
    character (500), target :: cbuff
    
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
             write(*,*) 'Running on ',nproc,' processor'
          else
             write(*,*) 'Running on ',nproc,' processors'
          end if
          write (*,*) 
          ! <doc> Call init_file_utils, ie. initialize the inputs and outputs, checking 
          !  whether we are doing a [[Trinity]] run or a list of runs </doc>
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

       if (proc0) then
          call time_message(.false., .false., time_init,' Initialization')
          cbuff = trim(run_name)
       end if
       
       call broadcast (cbuff)
       if (.not. proc0) run_name => cbuff
       
       call init_fields
       call init_gs2_diagnostics (list, nstep)
       call init_tstart (tstart)   ! tstart is in user units 
       
       if (present(dvdrho)) then
          if (proc0) then
             dvdrho = dvdrhon
             grho = surfarea/dvdrhon
          end if
          call broadcast (dvdrho)
          call broadcast (grho)
       end if
       
       if (proc0) call time_message(.false.,.false.,time_init,' Initialization')
       
       first_time = .false.

    else if (present(nensembles)) then
       if (nensembles > 1) then
          call scope (subprocs)
       end if
    end if
    
    istep_end = nstep
    
    call loop_diagnostics(0,exit)
    
    do istep = 1, nstep
       
       call advance (istep)
       
       if (nsave > 0 .and. mod(istep, nsave) == 0) &
            call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
       
       call update_time
       call loop_diagnostics (istep, exit)
       call check_time_step (reset, exit)
       if (proc0) call time_message(.false.,.true.,time_advance,' Advance time step')
       if (reset) call reset_time_step (istep, exit)
       
       if (mod(istep,5) == 0) call checkstop(exit)
       
       call checktime(avail_cpu_time,exit)
       
       if (exit) then
          istep_end = istep
          exit
       end if
    end do
    if (proc0) call write_dt

    time_interval = user_time-start_time

    if (present(nensembles)) then
       if (nensembles > 1) call ensemble_average (nensembles, time_interval)
    end if

    if (present(pflux)) then
       pflux = pflux_avg/time_interval
       qflux = qflux_avg/time_interval
       heat = heat_avg/time_interval
    else
       call finish_gs2_diagnostics (istep_end)
       call finish_gs2
    end if
    
    if (proc0) then
       call time_message(.false., .false., time_finish,'Finished run')
       time_total=time_init+time_advance+time_nc+time_reinit+time_finish
       print '(/,'' Initialization'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
            &'' Advance steps'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
            &'' Write restart'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
            &'' Re-initialize'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/, &
            &'' Finishing'',T25,0pf8.2,'' min'',T40,2pf5.1,'' %'',/,  &
            &'' total from timer is:'', 0pf9.2,'' min'',/)', &
            time_init/60.,time_init/time_total, &
            time_advance/60.,time_advance/time_total, &
            time_nc/60.,time_nc/time_total, &
            time_reinit/60.,time_reinit/time_total, &
            time_finish/60.,time_finish/time_total,time_total/60.
    endif
    
    if (.not. present(mpi_comm)) call finish_mp
    
  end subroutine run_gs2
  
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
    use mp, only: proc0
    use nonlinear_terms, only: finish_nonlinear_terms
    use run_parameters, only: finish_run_parameters
    use species, only: finish_species

    implicit none

    call finish_antenna
    call finish_collisions
    call finish_dist_fn
    call finish_fields
    call finish_hyper
    call finish_init_g
    call finish_kt_grids
    call finish_le_grids
    call finish_nonlinear_terms
    call finish_run_parameters
    call finish_species
    if (proc0) call finish_file_utils

  end subroutine finish_gs2

  subroutine reset_gs2 (ntspec, dens, temp, fprim, tprim, nu, nensembles)

    use dist_fn, only: d_reset => reset_init
    use collisions, only: vnmult, c_reset => reset_init
    use fields, only: init_fields, f_reset => reset_init
    use fields_implicit, only: fi_reset => reset_init
    use fields_explicit, only: fe_reset => reset_init
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
    use mp, only: proc0, scope, subprocs, allprocs

    implicit none

    integer, intent (in) :: ntspec, nensembles
    real, intent (in) :: dens, fprim
    real, dimension (:), intent (in) :: temp, tprim, nu

    integer :: istatus

    if (nensembles > 1) call scope (subprocs)

    call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    gnew = 0.

    call save_dt (code_dt)

    call d_reset
    call c_reset
    call f_reset
    call fi_reset
    call fe_reset
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

  subroutine gs2_trin_init (rhoc, qval, shat, aspr, kap, kappri, tri, tripri, shift, &
       betaprim, ntspec, dens, temp, fprim, tprim, nu)

    use species, only: init_trin_species
    use theta_grid_params, only: init_trin_geo

    implicit none

    integer, intent (in) :: ntspec
    real, intent (in) :: rhoc, qval, shat, aspr, kap, kappri, tri, tripri, dens, fprim, shift
    real, intent (in) :: betaprim
    real, dimension (:), intent (in) :: temp, tprim, nu

    call init_trin_species (ntspec, dens, temp, fprim, tprim, nu)
    call init_trin_geo (rhoc, qval, shat, aspr, kap, kappri, tri, tripri, shift, betaprim)
    
  end subroutine gs2_trin_init

# ifndef MAKE_LIB
end module gs2_main
# endif
