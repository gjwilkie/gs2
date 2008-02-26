module check

  implicit none
  public :: checkstop
  private

contains

  subroutine checkstop(exit)

    use mp, only: proc0, nproc, iproc, sum_allreduce
    use file_utils, only: get_unused_unit, run_name
    use run_parameters, only: margin
    logical, intent (in out) :: exit

    real :: timelimit, timeleft, timeused
    integer :: iret, unit
    integer, dimension(0:nproc-1) :: exit_array

! If any PE is running out of time, set exit flag
    exit_array = 0
!    iret=mpptime(timelimit, timeleft, timeused)
!    if (timeleft/timelimit < margin) exit_array(iproc) = 1

! If .stop file has appeared, set exit flag
    if (proc0) then
       call get_unused_unit (unit)
       open (unit=unit, file=trim(run_name)//".stop", status="old", err=100)
       exit_array(0) = 1
       close (unit=unit)
100    continue
    end if

    call sum_allreduce(exit_array)
    if (any(exit_array /= 0)) exit = .true.

  end subroutine checkstop

end module check
