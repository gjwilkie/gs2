#define CONCAT //

!> A program that tests the new diagnostics module. It  runs 
!! a  linear cyclone test case and then checks that the old and
!! new diagnostics give the same results
!!
!! This is free software released under GPLv3
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)

module checks_mod
  use unit_tests
  public checks
  contains
    function checks()
      logical :: checks
      checks = .true.
    end function checks
end module checks_mod

program test_gs2_diagnostics_new
  !use functional_tests
  !use checks_mod
  !call test_gs2('Linear CBC (unit test) to test new diagnostics', checks)
    use gs2_main, only: run_gs2, finish_gs2
    use unit_tests
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
    use gs2_diagnostics, only: finish_gs2_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
#endif
    implicit none
    integer :: n_vars 
    integer :: i

    character(len=40), dimension(200) :: variables, new_variables, n_lines
      !variables = (/'lambda', 'phi'/), &
      !n_lines = (/'3', '30'/)

    variables(1) = 'lambda'
    n_lines(1) = '4'
    variables(2) = 'phi'
    n_lines(2) = '78'
    variables(3) = 'kx'
    n_lines(3) = '3'
    variables(4) = 'phi2_by_mode'
    n_lines(4) = '10'
    variables(5) = 't'
    n_lines(5) = '5'
    variables(6) = 'phi2_by_ky'
    n_lines(6) = '12'
    variables(7) = 'phi2_by_kx'
    n_lines(7) = '12'
    variables(8) = 'phi2_by_kx'
    n_lines(8) = '8'
    variables(9) = 'gds21'
    n_lines(9) = '4'
    variables(10) = 'bmag'
    n_lines(10) = '4'
    new_variables(1:10) = variables(1:10)

    variables(11) = 'es_heat_by_k'
    new_variables(11) = 'es_heat_flux_by_mode'
    n_lines(11) = '20'
    variables(12) = 'hflux_tot'
    new_variables(12) = 'heat_flux_tot'
    n_lines(12) = '4'
    variables(13) = 'es_mom_flux'
    new_variables(13) = 'es_mom_flux'
    n_lines(13) = '4'

    variables(14) = 'omega'
    new_variables(14) = 'omega'
    n_lines(14) = '4'
    variables(15) = 'omegaavg'
    new_variables(15) = 'omega_average'
    n_lines(15) = '4'

    n_vars = 15



    call init_mp

    test_driver_flag = .true.
    functional_test_flag = .true.

   call announce_module_test("gs2_diagnostics_new")

    call run_gs2(mp_comm)


    !call announce_test('results')
    !call process_test(test_function(), 'results')

    call finish_gs2_diagnostics(ilast_step)
#ifdef NEW_DIAG
    call finish_gs2_diagnostics_new
#endif
    call finish_gs2

    if (proc0) then
      do i = 1,n_vars
        call announce_test("value of "//trim(new_variables(i)))
        call process_test(test_variable(trim(variables(i)), trim(new_variables(i)), &
          trim(n_lines(i))), &
          "value of "//trim(new_variables(i)))
      end do
    end if


   call close_module_test("gs2_diagnostics_new")

    call finish_mp

contains
  
  function test_variable(var_name, new_var_name, n_lines)
    use unit_tests, only: should_print
    character(*), intent(in) :: var_name, new_var_name, n_lines
    logical :: test_variable
    character(200) ::  command 
    
    test_variable=.true.
#ifdef NEW_DIAG
    command = "if [ ""`ncdump -v "//var_name//" test_gs2_diagnostics_new.out.nc  | tail -n "//n_lines//"`"" = &
     &  ""`ncdump -v "//new_var_name//" test_gs2_diagnostics_new.cdf | tail -n "//n_lines//" `"" ]; &
     & then echo ""T"" > test_tmp.txt; fi"
    
    if (should_print(3)) write(*,*) command
    call system(" echo ""F"" > test_tmp.txt")
    call system(command)
    test_variable = .true.
    open(120349, file='test_tmp.txt')
    read(120349, '(L)') test_variable
    close(120349)
#endif
  end function test_variable


end program test_gs2_diagnostics_new
