!> A module which contains a series of high level tests used in the
!! linear and nonlinear test cases, as well as a driver function 
!! which runs GS2.
module functional_tests
  use unit_tests, only: announce_check, agrees_with, process_check
  use unit_tests, only: process_test, announce_test, print_with_stars
  use runtime_tests, only: verbosity

  implicit none

  private

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

  public :: check_growth_rates_equal_in_list
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
    use unit_tests, only: ilast_step
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

  function check_growth_rates_equal_in_list(err)
    use unit_tests, only: ilast_step
    use gs2_diagnostics, only: get_omegaavg
    use kt_grids, only: ntheta0, naky
    use mp, only: proc0, grp0, scope, subprocs, allprocs
    use mp, only: iproc, nproc, sum_allreduce
    use kt_grids, only: aky
    use run_parameters, only: wstar_units
    use job_manage, only: njobs
    implicit none
    real, intent(in) :: err
    real, dimension(naky, njobs) :: omegas
    logical :: check_growth_rates_equal_in_list

    logical :: dummy
    logical :: check_result
    integer :: i, ijob
    complex, dimension(ntheta0, naky) :: omegaavg

    check_growth_rates_equal_in_list = .true.

    omegas=0.
    omegaavg=0.

    if (proc0) then
       !call announce_check('growth rate')
       call get_omegaavg(ilast_step-1, dummy, omegaavg)
       !check_result =  agrees_with(aimag(omegaavg(1,:)), rslt, err)
       !call process_check(check_growth_rate, check_result, 'growth rate')
       call scope(allprocs)
       ! work out which job in the list we are
       do i = 1,njobs
          if (iproc==grp0(i-1)) ijob = i 
       end do
       ! Get the omegas from this job
       omegas(:,ijob) = aimag(omegaavg(1,:))
       if (wstar_units) then
          omegas(:, ijob) = omegas(:,ijob) * aky(:) / 2.0
       end if
    else
       call scope(allprocs)
    end if

    call sum_allreduce(omegas)
    call scope(subprocs)

    do i = 2,njobs
       call announce_check('growth rate')
       check_result = agrees_with(omegas(:,1), omegas(:,i), err)
       call process_check( &
            check_growth_rates_equal_in_list, check_result, 'growth rate')
    end do

  end function check_growth_rates_equal_in_list

  subroutine test_gs2(test_name, test_function)
    use gs2_main, only: run_gs2, finish_gs2
    use unit_tests, only: functional_test_flag, ilast_step
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
    use gs2_diagnostics, only: finish_gs2_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
#endif
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
#ifdef NEW_DIAG
    call finish_gs2_diagnostics_new
#endif
    call finish_gs2

    if (proc0) call close_functional_test(test_name)

    call finish_mp
  end subroutine test_gs2
end module functional_tests

