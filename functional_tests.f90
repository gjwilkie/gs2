
!> A module which contains a series of high level tests used in the
!! linear and nonlinear test cases.
module functional_tests
  use unit_tests

  public :: test_growth_rate
  public :: announce_functional_test
  public :: close_functional_test


contains

  subroutine announce_functional_test(test_name)
    character(*), intent(in) :: test_name
    call print_with_stars('Starting functional test: ', test_name)
  end subroutine announce_functional_test
  subroutine close_functional_test(test_name)
    character(*), intent(in) :: test_name
    call print_with_stars('Completed functional test: ', test_name)
  end subroutine close_functional_test
  subroutine test_growth_rate(istep, rslt, err)
    use gs2_diagnostics, only: get_omegaavg
    use kt_grids, only: ntheta0, naky
    use mp, only: proc0
    integer, intent(in) :: istep
    real, dimension(naky), intent(in) :: rslt
    real, intent(in) :: err

    logical :: dummy
    complex, dimension(ntheta0, naky) :: omegaavg

    if (proc0) then
      call announce_test('growth rate')
      call get_omegaavg(istep, dummy, omegaavg)
      call process_test(agrees_with(aimag(omegaavg(1,:)), rslt, err), 'growth rate')
    end if


  end subroutine test_growth_rate

end module functional_tests

