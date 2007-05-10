module check

  implicit none
  public :: checkstop
  private

contains

  subroutine checkstop(exit)

! in this version, margin refers to the number of minutes rather 
! than the fraction of run time.

    use mp, only: proc0, nproc, iproc, sum_allreduce, broadcast
    use file_utils, only: get_unused_unit, run_name
    use run_parameters, only: margin
    logical, intent (in out) :: exit
    integer :: unit
    real :: left
!    integer, dimension (0:nproc-1) :: exit_array

!    exit_array = 0
!!! just look at one PE for now
    if (proc0) then
!       call timeleft (trim(run_name), left)
!       if (left < margin) exit_array(iproc) = 1
!!!       if (left < margin) exit = .true.
    endif

! If .stop file has appeared, set exit flag
    if (proc0) then
       call get_unused_unit (unit)
       open (unit=unit, file=trim(run_name)//".stop", status="old", err=100)
!       exit_array(0) = 1
       exit = .true.
       close (unit=unit)
100    continue
    end if

!    call sum_allreduce(exit_array)
!    if (any(exit_array /= 0)) exit = .true.
   
    call broadcast (exit)

  end subroutine checkstop

  subroutine timeleft (runname, time_left)

    character (*), intent (in) :: runname
    real, intent (out) :: time_left

    integer, external :: system
    character (500) :: cmd,file,str
    real :: hr, min, sec
    integer :: i
    
    file=trim(runname)//'.time'
    cmd = 'llqs | grep '//trim(runname)//'>'//trim(runname)//'.time'
    i = system(trim(cmd))

    open (unit=12,file=trim(file), err=100)
    read (12,*) str,str,str,str,str,str,str
    close(12)

    hr = (ichar(str(1:1))-48)*10 + (ichar(str(2:2))-48)
    min = (ichar(str(4:4))-48)*10 + (ichar(str(5:5))-48)
    sec = (ichar(str(7:7))-48)*10 + (ichar(str(8:8))-48)

    time_left = ((hr*60+min)*60+sec)/60
    return

100 time_left = 1.e8

  end subroutine timeleft

end module check
