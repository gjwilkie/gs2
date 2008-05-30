module job_manage

  implicit none
  private
  public :: timer
  public :: job_fork
  public :: checkstop

contains
  
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

  subroutine job_fork

    use file_utils, only: get_unused_unit, list_name, run_name, init_job_name
    use mp, only: job
    use mp, only: proc0, nproc
    use mp, only: init_jobs, broadcast, finish_mp
    implicit none
    integer, dimension(:), allocatable :: group0
    integer :: i, njobs

    character (len=500), dimension(:), allocatable :: job_list

    integer :: list_unit, ierr

    if (proc0) then
       call get_unused_unit(list_unit)
       open (unit=list_unit, file=trim(list_name))
       read (list_unit,*) njobs
    end if
    call broadcast (njobs)
    
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

    if (proc0) then
       do i=0,njobs-1
          read (list_unit, fmt="(a)") job_list(i)
       end do
       close (list_unit)
    end if

    do i=0,njobs-1
       call broadcast (job_list(i))
    end do

    allocate (group0(0:njobs-1))

    call init_jobs (njobs, group0, ierr)
!    call init_job_name (njobs, group0, job_list)
    call init_job_name (job_list(job))

    if (nproc > 1 .and. proc0) &
         & write(*,*) 'Job ',job,' is called ',trim(run_name),&
         & ' and is running on ',nproc,' processors'
    if (nproc == 1) write(*,*) 'Job ',job,' is called ',trim(run_name),&
         & ' and is running on ',nproc,' processor'

  end subroutine job_fork

  subroutine checkstop(exit)

    use mp, only: proc0, broadcast
    use file_utils, only: get_unused_unit, run_name
    logical, intent (in out) :: exit
    integer :: unit

! If .stop file has appeared, set exit flag
    if (proc0) then
       call get_unused_unit (unit)
       open (unit=unit, file=trim(run_name)//".stop", status="old", err=100)
       exit = .true.
       close (unit=unit)
100    continue
    end if

    call broadcast (exit)

  end subroutine checkstop

end module job_manage
