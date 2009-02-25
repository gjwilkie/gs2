program gs2
  use job_manage, only: checkstop, job_fork, checktime
  use mp, only: init_mp, finish_mp, proc0, nproc, broadcast
!  use mp, only: init_mp, finish_mp, proc0, iproc, nproc, broadcast, init_jobs
!  use mp, only: scope, allprocs, subprocs, send, receive, barrier, job
!  use mp, only: send, receive, barrier, job
  use file_utils, only: init_file_utils, finish_file_utils, run_name, list_name
  use fields, only: init_fields
  use gs2_diagnostics, only: init_gs2_diagnostics, finish_gs2_diagnostics
  use gs2_diagnostics, only: nsave
  use run_parameters, only: nstep
  use run_parameters, only: fphi, fapar, fbpar
  use run_parameters, only: avail_cpu_time
  use fields, only: advance
  use dist_fn_arrays, only: gnew
  use gs2_save, only: gs2_save_for_restart
  use gs2_diagnostics, only: loop_diagnostics
  use gs2_reinit, only: reset_time_step, check_time_step
  use gs2_reinit, only: time_message, time_nc, time_reinit
  use gs2_time, only: update_time
  use gs2_time, only: write_dt, init_tstart
  use gs2_time, only: user_time, user_dt
  use gs2_time, only: code_time
  use init_g, only: tstart
 ! use check, only: checkstop
  use collisions, only: vnmult

  implicit none
  real :: time_init = 0., time_advance = 0., time_finish = 0., time_total
  integer :: istep = 0, istep_end, unit, istatus
  logical :: exit, reset, list
  character (500), target :: cbuff

! initialize message passing 
  call init_mp
  call checktime(avail_cpu_time,exit) ! initialize timer

! report # of processors being used
  if (proc0) then
     if (nproc == 1) then
        write(*,*) 'Running on ',nproc,' processor'
     else
        write(*,*) 'Running on ',nproc,' processors'
     end if
     write (*,*) 
! figure out run name or get list of jobs
     call init_file_utils (list, name="gs")
  end if

  call broadcast (list)

! if given a list of jobs, fork
  if (list) call job_fork

  if (proc0) then
     call time_message(.false., .false., time_init,' Initialization')
     cbuff = trim(run_name)
  end if

  call broadcast (cbuff)
  if (.not. proc0) run_name => cbuff

  call init_fields
  call init_gs2_diagnostics (list, nstep)
  call init_tstart (tstart)   ! tstart is in user units 
  if (proc0) call time_message(.false.,.false.,time_init,' Initialization')
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

  call finish_gs2_diagnostics (istep_end)
  if (proc0) call finish_file_utils
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

  call finish_mp

!!$contains
!!$
!!$  subroutine timer (i)
!!$    
!!$    character (len=10) :: zdate, ztime, zzone
!!$    integer, dimension(8) :: ival
!!$    real, dimension (:), allocatable, save :: tsave
!!$    real, save :: told=0., tnew=0., tavg
!!$    integer :: i
!!$    integer :: navg=10, j=0
!!$    logical :: first_time = .true.
!!$    
!!$    if (first_time) then
!!$       allocate (tsave(0:navg-1))
!!$       tsave = 0.
!!$       first_time = .false.
!!$    end if
!!$
!!$    call date_and_time (zdate, ztime, zzone, ival)
!!$    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
!!$    if (i == 1) then
!!$       tsave(mod(j,navg)) = tnew-told
!!$       j=j+1
!!$       if (j>=navg) then
!!$          tavg = sum(tsave)/real(navg)
!!$          if (abs((tnew-told)/tavg-1.) > 2.0) then
!!$             print *, 'Avg time = ',tavg, &
!!$                  & '  Time since last called: ',tnew-told,' seconds'
!!$          end if
!!$       else if (told > 0.) then
!!$          print *, 'Time since last called: ',tnew-told,' seconds'
!!$       end if
!!$    end if
!!$    told = tnew
!!$  end subroutine timer
!!$
!!$  subroutine job_fork
!!$
!!$    use file_utils
!!$    implicit none
!!$    integer, dimension(:), allocatable :: group0
!!$    integer :: i, njobs
!!$
!!$    character (len=500), dimension(:), allocatable :: job_list
!!$
!!$    integer :: list_unit, ierr
!!$
!!$    if (proc0) then
!!$       call get_unused_unit(list_unit)
!!$       open (unit=list_unit, file=trim(list_name))
!!$       read (list_unit,*) njobs
!!$    end if
!!$    call broadcast (njobs)
!!$    
!!$    if (nproc < njobs) then
!!$       if (proc0) then
!!$          write (*,*) 
!!$          write (*,*) 'Number of jobs = ',njobs,' and number of processors = ',nproc
!!$          write (*,*) 'Number of processors must not be less than the number of jobs'
!!$          write (*,*) 'Stopping'
!!$          write (*,*) 
!!$       end if
!!$       call finish_mp
!!$       stop
!!$    end if
!!$
!!$    if (mod(nproc, njobs) /= 0) then
!!$       if (proc0) then
!!$          write (*,*) 
!!$          write (*,*) 'Number of jobs = ',njobs,' and number of processors = ',nproc
!!$          write (*,*) 'Number of jobs must evenly divide the number of processors.'
!!$          write (*,*) 'Stopping'
!!$          write (*,*) 
!!$       end if
!!$       call finish_mp
!!$       stop
!!$    end if
!!$
!!$    allocate (job_list(0:njobs-1))
!!$
!!$    if (proc0) then
!!$       do i=0,njobs-1
!!$          read (list_unit, fmt="(a)") job_list(i)
!!$       end do
!!$       close (list_unit)
!!$    end if
!!$
!!$    do i=0,njobs-1
!!$       call broadcast (job_list(i))
!!$    end do
!!$
!!$    allocate (group0(0:njobs-1))
!!$
!!$    call init_jobs (njobs, group0, ierr)
!!$!    call init_job_name (njobs, group0, job_list)
!!$    call init_job_name (job_list(job))
!!$
!!$    if (nproc > 1 .and. proc0) &
!!$         & write(*,*) 'Job ',job,' is called ',trim(run_name),&
!!$         & ' and is running on ',nproc,' processors'
!!$    if (nproc == 1) write(*,*) 'Job ',job,' is called ',trim(run_name),&
!!$         & ' and is running on ',nproc,' processor'
!!$
!!$  end subroutine job_fork

end program gs2
