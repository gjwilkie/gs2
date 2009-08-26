# include "define.inc"

module job_manage

  implicit none
  private
  public :: timer
  public :: timer_local
  public :: job_fork
  public :: checkstop
  public :: checktime
! MAB> -- needed for Trinity
  public :: njobs

  integer :: njobs
! <MAB
  logical :: debug=.false.

contains

!!! returns CPU time in second
  function timer_local()
# ifdef MPI
# ifndef MPIINC
    use mpi, only: mpi_wtime
# else
    include "mpif.h" ! CMR following Michele Weiland's advice
# endif

# endif
    real :: timer_local

    timer_local=0.

# ifdef MPI    
    timer_local=mpi_wtime()
# else
    ! ther compilers are not tested, but may work
    ! other standard timer routines are clock, etime, etc
# if ( FCOMPILER == _XL_ || FCOMPILER == _GFORTRAN_ || FCOMPILER == _G95_ \
    || FCOMPILER == _INTEL_ || FCOMPILER == _PGI_ \
    || FCOMPILER == _PATHSCALE_ )
    ! this routine is F95 starndard
    call cpu_time(timer_local)
# endif
# endif

  end function timer_local
    
  subroutine timer (i)

    integer, intent(in) :: i
    character (len=10) :: zdate, ztime, zzone
    integer, dimension(8) :: ival
    real, dimension (:), allocatable, save :: tsave
    real, save :: told=0., tnew=0., tavg
    integer :: navg=10, j=0
    logical :: first_time = .true.

    if (first_time) then
       allocate (tsave(0:navg-1))
       tsave = 0.
       first_time = .false.
    end if

    call date_and_time (zdate, ztime, zzone, ival)
    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
    if (i == 1) then
       tsave(mod(j,navg)) = tnew-told
       j=j+1
       if (j>=navg) then
          tavg = sum(tsave)/real(navg)
          if (abs((tnew-told)/tavg-1.) > 2.0) then
             print *, 'Avg time = ',tavg, &
                  & '  Time since last called: ',tnew-told,' seconds'
          end if
       else if (told > 0.) then
          print *, 'Time since last called: ',tnew-told,' seconds'
       end if
    end if
    told = tnew
  end subroutine timer

  subroutine job_fork (nrad, ngrads, fork_flag)

    use file_utils, only: get_unused_unit, list_name, run_name, init_job_name
! MAB> -- moved init_error_unit and init_input_unit calls here from file_utils
! because they were being called there on all procs when they should be called
! only on proc0
!    use file_utils, only: init_error_unit, init_input_unit, init_nt_name, trin
    use file_utils, only: init_error_unit, init_input_unit, trin, list_name
! <MAB
    use mp, only: job, scope, allprocs
    use mp, only: proc0, nproc
    use mp, only: init_jobs, broadcast, finish_mp
    implicit none
    integer, dimension(:), allocatable :: group0
    integer :: i, l

    character (len=500) :: ext, hack ! needed for trinity
    character (len=500), dimension(:), allocatable :: job_list

    integer :: list_unit, ierr
    logical :: err = .true., inp = .true., fork = .true. ! MAB

    integer, intent (in), optional :: nrad, ngrads
    logical, intent (in), optional :: fork_flag

    if (trin) then
       ! number of files to run computed from # radial grid points x # evolved grads
       ! given by trinity transport input file
       njobs = nrad*(ngrads+1)
    else
       ! open file containing list of input files to run and read total
       ! number of input files (i.e. # flux tubes x # evolved grads) from first line
       if (proc0) then
          call get_unused_unit(list_unit)
          open (unit=list_unit, file=trim(list_name))
          read (list_unit,*) njobs
       end if
       call broadcast (njobs)
    end if

    if (present(fork_flag)) then
       if (.not. fork_flag) then
          fork = .false.
       end if
    end if

    if (fork) then
       
       if (nproc < njobs) then
          if (proc0) then
             write (*,*) 
             write (*,*) 'Number of jobs = ',njobs,' and number of processors = ',nproc
             write (*,*) 'Number of processors must not be less than the number of jobs'
             write (*,*) 'Stopping'
             write (*,*) 
          end if
          call finish_mp
          stop
       end if
       
       if (mod(nproc, njobs) /= 0) then
          if (proc0) then
             write (*,*) 
             write (*,*) 'Number of jobs = ',njobs,' and number of processors = ',nproc
             write (*,*) 'Number of jobs must evenly divide the number of processors.'
             write (*,*) 'Stopping'
             write (*,*) 
          end if
          call finish_mp
          stop
       end if
       
       allocate (job_list(0:njobs-1))
       
       if (trin) then
          if (proc0) then
!          l = len_trim(trin_name)-5
             l = len_trim(list_name)-5
             do i=0,njobs-1
                write (ext,*) i+1
                ext = adjustl(ext)
                job_list(i) = list_name(:l)//trim(ext)
             end do
          end if
       else
          if (proc0) then
             do i=0,njobs-1
                read (list_unit, fmt="(a)") job_list(i)
             end do
          
          ! read in name of input file for profile evolution
!          read (list_unit, fmt="(a)") nt_name
             close (list_unit)
          end if
       end if

       do i=0,njobs-1
          call broadcast (job_list(i))
       end do
!       call broadcast (nt_name)
!       call init_nt_name (nt_name)
       
       allocate (group0(0:njobs-1))

       call init_jobs (njobs, group0, ierr)
       ! TT> brought up one line [call scope(subprocs)] from file_utils.fpp
       !     to init_jobs
       !    call init_job_name (njobs, group0, job_list)
       call init_job_name (job_list(job))
       ! <TT

       ! MAB> moved from file_utils because had to be within proc0, 
       ! which is undefined there
       if (proc0) then
          call init_error_unit (err)
          call init_input_unit (inp)
       end if
       ! <MAB

       if (nproc > 1 .and. proc0) &
            & write(*,*) 'Job ',job,' is called ',trim(run_name),&
            & ' and is running on ',nproc,' processors'
       if (nproc == 1) write(*,*) 'Job ',job,' is called ',trim(run_name),&
            & ' and is running on ',nproc,' processor'
       
       deallocate (group0, job_list) ! MAB
       
       if (trin) call scope (allprocs)
       
    end if

  end subroutine job_fork

  subroutine checkstop(exit,list)

    use mp, only: proc0, broadcast
    use file_utils, only: run_name, list_name
    logical, intent (in), optional :: list
    logical, intent (in out) :: exit
    character (len=300) :: filename
    
    logical :: exit_local

    ! If .stop file has appeared, set exit flag
    filename=trim(run_name)//".stop"
    if(present(list)) then
       if(list) filename=list_name(:len_trim(list_name)-5)//".stop"
    endif
    
    if (proc0) then
       inquire(file=filename,exist=exit_local)
       exit = exit .or. exit_local
    end if

    call broadcast (exit)

  end subroutine checkstop

  subroutine checktime(avail_time,exit)
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit

    ! available time in second
    real, intent(in) :: avail_time
    ! true if elapse time exceed available time
    logical, intent(in out) :: exit
    logical, save :: initialized=.false.
    real :: elapse_time=0.
    real :: initial_time=0.
    real :: margin=300. ! 5 minutes

    if(.not.initialized) then
       initial_time=timer_local()  ! timer_local() returns #seconds from fixed time in past
       initialized=.true.
       return
    endif

    elapse_time=timer_local()-initial_time

    if(proc0) then
       if(elapse_time >= avail_time-margin) then
          write(error_unit(),'(a,f12.4,a,f12.4)') &
               & 'Elapse time ',elapse_time, &
               & ' exceeds available time',avail_time-margin
          write(error_unit(),'(a,f12.4,a,f12.4,a)') &
               & '  (Given CPU time: ',avail_time, &
               & '  Margin: ',margin,')'
          exit=.true.
       endif
    endif

    call broadcast(exit)

  end subroutine checktime

end module job_manage
