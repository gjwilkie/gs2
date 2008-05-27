module check

  implicit none
  public :: checkstop
  private

contains

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

end module check
