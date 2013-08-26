
!> A module which contains a series of high level tests used in the
!! linear and nonlinear test cases, as well as a driver function 
!! which runs GS2.
module functional_tests
  use unit_tests

  !> Check that the relative error of the growth rates as a 
  !! function of ky with respect to the given rslt is less than
  !! err. rslt should be an array of length naky which contains
  !! the growth rates as a function of ky.
  !! Returns .true. for pass and .false. for fail
  public :: check_growth_rate
 
  !> Run gs2 and then call the test_function to check the results 
  !! corresponding to the input file provided. test_name should 
  !! be a character(*) containing the title of the tests, and 
  !! test_function should be a logical function which returns
  !! .true. for pass and false for .fail.
  public :: test_gs2


contains

  subroutine announce_functional_test(test_name)
    character(*), intent(in) :: test_name
    if (verbosity() .gt. 0) call print_with_stars('Starting functional test: ', test_name)
  end subroutine announce_functional_test
  subroutine close_functional_test(test_name)
    character(*), intent(in) :: test_name
    if (verbosity() .gt. 0) call print_with_stars('Completed functional test: ', test_name)
  end subroutine close_functional_test
  function check_growth_rate(rslt, err)
    use gs2_main, only: ilast_step
    use gs2_diagnostics, only: get_omegaavg
    use kt_grids, only: ntheta0, naky
    use mp, only: proc0
    real, dimension(naky), intent(in) :: rslt
    real, intent(in) :: err
    logical :: check_growth_rate

    logical :: dummy
    logical :: check_result
    complex, dimension(ntheta0, naky) :: omegaavg

    check_growth_rate = .true.
  

    if (proc0) then
      call announce_check('growth rate')
      call get_omegaavg(ilast_step-1, dummy, omegaavg)
      check_result =  agrees_with(aimag(omegaavg(1,:)), rslt, err)
      call process_check(check_growth_rate, check_result, 'growth rate')
    end if


  end function check_growth_rate

  subroutine test_gs2(test_name, test_function)
    use gs2_main, only: run_gs2, functional_test_flag, ilast_step, finish_gs2
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
    use gs2_diagnostics, only: finish_gs2_diagnostics
    implicit none
    character(*), intent(in) :: test_name
    logical, external :: test_function


    call init_mp

    test_driver_flag = .true.
    functional_test_flag = .true.

    if (proc0) call announce_functional_test(test_name)

    call run_gs2(mp_comm)


    call announce_test('results')
    call process_test(test_function(), 'results')

    call finish_gs2_diagnostics(ilast_step)
    call finish_gs2

    if (proc0) call close_functional_test(test_name)

    call finish_mp
  end subroutine test_gs2

end module functional_tests

