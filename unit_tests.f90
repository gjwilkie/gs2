
!> A module containing functions for running unit tests
module unit_tests

  !> Is the relative error of the first argument with respect
  !! to the correct answer less than err? 
  public :: agrees_with

  !> Announce the start of the given test (to let the user know 
  !! which one failed!)
  public :: announce_test

  !> Take the logical result of the test (true for pass), and 
  !! either announce its success or announce its failure then stop.
  public :: process_test

contains

  function agrees_with(val, correct, err)
    real, intent(in) :: val, correct, err
    logical :: agrees_with
    write (*,*) val, ' should be ', correct 
    agrees_with = (abs((val-correct)/correct) .lt. err)
  end function agrees_with

  subroutine announce_test(test_name)
    character(*), intent(in) :: test_name
    write (*,*) '--> Testing ', test_name
  end subroutine announce_test
  subroutine process_test(rslt, test_name)
    logical, intent (in) :: rslt
    character(*), intent(in) :: test_name
    if (.not. rslt) then 
      write(*,*) '--> ', test_name, ' failed'
      stop 1
    end if

    write (*,*) '--> ', test_name, ' passed'
  end subroutine process_test

end module unit_tests
