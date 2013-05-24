
!> A module containing functions for running unit tests
module unit_tests

  !> Is the relative error of the first argument with respect
  !! to the correct answer less than err? 
  !! If the correct value is 0.0, then err is treated as an 
  !! absolute error
  public :: agrees_with

  !> Announce the start of the given test (to let the user know 
  !! which one failed!)
  public :: announce_test

  !> Take the logical result of the test (true for pass), and 
  !! either announce its success or announce its failure then stop.
  public :: process_test

  public :: announce_module_test

contains

  function agrees_with(val, correct, err)
    real, intent(in) :: val, correct, err
    logical :: agrees_with
    write (*,*) val, ' should be ', correct 
    if (correct .eq. 0.0) then 
      agrees_with = abs(val) .lt. err 
    else
      agrees_with = (abs((val-correct)/correct) .lt. err)
    end if
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

  subroutine announce_module_test(module_name)
    character(*), intent(in) :: module_name
    character(8) :: message = 'Testing '
    call print_with_stars(message, module_name)


  end subroutine announce_module_test

  subroutine close_module_test(module_name)
    character(*), intent(in) :: module_name
    character(17) :: message = 'Finished testing '
    call print_with_stars(message, module_name)
    write (*,*)


  end subroutine close_module_test

  subroutine print_with_stars(str1, str2)
    character(*), intent (in) :: str1, str2
    character, dimension(:), allocatable :: stars
    integer :: i

    allocate(stars(len(str1) + len(str2) ))
    do i = 1,len(str1)+len(str2)
      stars(i) = '*'
    end do

    write (*,*) stars
    write (*,*) str1, str2
    write (*,*) stars
    
  end subroutine print_with_stars

end module unit_tests
