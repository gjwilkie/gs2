
!> A module containing functions for running unit tests
module unit_tests

  implicit none 
  private

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


  !> Announce the start of the given partial test (to let the user know 
  !! which one failed!)
  public :: announce_check

  !> Take the logical result of the check (partial test) (true for pass), and 
  !! either announce its success or announce its failure.
  public :: process_check

  public :: announce_module_test
  public :: close_module_test

  interface agrees_with
    module procedure agrees_with_real
    module procedure agrees_with_real_1d_array
    module procedure agrees_with_integer
  end interface

  interface should_be
    module procedure should_be_int
    module procedure should_be_real
  end interface



contains

  function proc_message()
    use mp, only: iproc
    character(16) :: proc_message
    write(proc_message, '(A10,I2)') ' on iproc ', iproc
  end function proc_message

  function agrees_with_integer(val, correct)
    integer, intent(in) :: val, correct
    logical :: agrees_with_integer
    call should_be(val, correct)
    agrees_with_integer = (val .eq. correct)
  end function agrees_with_integer

  function agrees_with_real_1d_array(val, correct, err)
    real, dimension(:), intent(in) :: val, correct
    real, intent(in) :: err
    logical :: agrees_with_real_1d_array
    integer :: n, i
    n = size(val)
    agrees_with_real_1d_array = .true.
    do i = 1,n

      call should_be(val(i), correct(i))
      if (correct(i) .eq. 0.0) then 
        agrees_with_real_1d_array = agrees_with_real_1d_array .and. (abs(val(i)) .lt. err) 
      else
        agrees_with_real_1d_array = agrees_with_real_1d_array .and. &
                                        (abs((val(i)-correct(i))/correct(i)) .lt. err)
      end if
      if (.not. agrees_with_real_1d_array) exit
    end do 
  end function agrees_with_real_1d_array
  function agrees_with_real(val, correct, err)
    real, intent(in) :: val, correct, err
    logical :: agrees_with_real
    call should_be(val, correct)
    if (correct .eq. 0.0) then 
      agrees_with_real = abs(val) .lt. err 
    else
      agrees_with_real = (abs((val-correct)/correct) .lt. err)
    end if
  end function agrees_with_real

  subroutine should_be_int(val, rslt)
    integer, intent(in) :: val, rslt
    write (*,*) '      Value: ', val, ' should be ', rslt, proc_message()
  end subroutine should_be_int

  subroutine should_be_real(val, rslt)
    real, intent(in) :: val, rslt
    write (*,*) '      Value: ', val, ' should be ', rslt, proc_message()
  end subroutine should_be_real

  subroutine announce_test(test_name)
    character(*), intent(in) :: test_name
    write (*,*) '--> Testing ', test_name, proc_message()
  end subroutine announce_test

  subroutine process_test(rslt, test_name)
    logical, intent (in) :: rslt
    character(*), intent(in) :: test_name
    if (.not. rslt) then 
      write(*,*) '--> ', test_name, ' failed', proc_message()
      stop 1
    end if

    write (*,*) '--> ', test_name, ' passed', proc_message()
  end subroutine process_test

  subroutine announce_check(test_name)
    character(*), intent(in) :: test_name
    write (*,*) '   --> Checking ', test_name, proc_message()
  end subroutine announce_check

  subroutine process_check(rslt, test_name)
    logical, intent (in) :: rslt
    character(*), intent(in) :: test_name
    if (.not. rslt) then 
      write(*,*) '   --> ', test_name, ' failed', proc_message()
    else
      write (*,*) '   --> ', test_name, ' passed', proc_message()
    end if
  end subroutine process_check

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
    use mp, only: iproc
    character(*), intent (in) :: str1, str2
    character, dimension(:), allocatable :: stars
    integer :: i


    allocate(stars(len(str1) + len(str2) + len(proc_message()) ))
    do i = 1,len(str1)+len(str2)+len(proc_message())
      stars(i) = '*'
    end do

    write (*,*) stars
    write (*,*) str1, str2, proc_message()
    write (*,*) stars
    
  end subroutine print_with_stars

end module unit_tests
