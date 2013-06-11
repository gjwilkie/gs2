
!> A module which contains a series of high level tests used in the
!! linear and nonlinear test cases.
module functional_tests
  use unit_tests

  public :: check_growth_rate
  public :: announce_functional_test
  public :: close_functional_test
  public :: test_gs2


contains

  subroutine announce_functional_test(test_name)
    character(*), intent(in) :: test_name
    call print_with_stars('Starting functional test: ', test_name)
  end subroutine announce_functional_test
  subroutine close_functional_test(test_name)
    character(*), intent(in) :: test_name
    call print_with_stars('Completed functional test: ', test_name)
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
    complex, dimension(ntheta0, naky) :: omegaavg

    check_growth_rate = .true.

    if (proc0) then
      call announce_check('growth rate')
      call get_omegaavg(ilast_step-1, dummy, omegaavg)
      check_growth_rate = check_growth_rate .and. agrees_with(aimag(omegaavg(1,:)), rslt, err)
      call process_check(check_growth_rate, 'growth rate')
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

