#define CONCAT //

!> A program that tests the new diagnostics module. It  runs 
!! a  linear cyclone test case and then checks that the old and
!! new diagnostics give the same results
!!
!! This is free software released under the MIT license
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
    use unit_tests, only: should_print
#ifdef NEW_DIAG
    use diagnostics_config, only: override_screen_printout_options
#endif
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp
    use mp, only: broadcast
    use gs2_diagnostics, only: finish_gs2_diagnostics
    use gs2_diagnostics, only: pflux_avg, qflux_avg, heat_avg, vflux_avg
    use gs2_diagnostics, only: diffusivity
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
    use gs2_diagnostics_new, only: gnostics
#endif
    use run_parameters, only: use_old_diagnostics, nstep
    use species, only: nspec
    implicit none
    integer :: n_vars, n_file_names
    integer :: i
    real :: eps
    real, dimension(:), allocatable :: pfluxav, qfluxav, heatav, vfluxav
    real :: diff
    !real :: vfluxav

    character(len=40), dimension(200) :: variables, new_variables, n_lines, file_names, n_lines_files
      !variables = (/'lambda', 'phi'/), &
      !n_lines = (/'3', '30'/)
  ! General config


  eps = 1.0e-7

  if (precision(eps).lt. 11) eps = eps * 100.0

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
    variables(6) = 'phi2_by_ky -p 7,11'
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
    variables(13) = 'vflux_tot'
    new_variables(13) = 'mom_flux_tot'
    n_lines(13) = '4'

    variables(14) = 'omega'
    new_variables(14) = 'omega'
    n_lines(14) = '4'
    variables(15) = 'omegaavg'
    new_variables(15) = 'omega_average'
    n_lines(15) = '4'
    
    variables(16) = 'ntot00  -p 4,8'
    new_variables(16) = 'ntot_flxsurf_avg  -p 4,8'
    n_lines(16) = '20'

    variables(17) = 'phi0'
    new_variables(17) = 'phi_igomega_by_mode'
    n_lines(17) = '15'

    variables(18) = 'apar_heat_by_k'
    new_variables(18) = 'apar_heat_flux_by_mode'
    n_lines(18) = '20'

    n_vars = 18


      ! Note that heat and heat2 no longer pass because the comparison
      ! is no longer between files output from the same run but between
      ! files output from two successive runs. This causes changes greater
      ! than the precision of the heat and heat2 files. EGH
      ! Thus they are skipped.
    file_names(1) = 'heat'
    n_lines_files(1) = '16'
    file_names(2) = 'heat2'
    n_lines_files(2) = '50'
    file_names(3) = 'lpc'
    n_lines_files(3) = '16'
    file_names(4) = 'vres'
    n_lines_files(4) = '16'
    file_names(5) = 'phase'
    n_lines_files(5) = '16'
    file_names(6) = 'jext'
    n_lines_files(6) = '0'
    file_names(7) = 'parity'
    n_lines_files(7) = '34'

    n_file_names = 7

    call init_mp

  ! Here we switch on print_line and print_flux_line if we have 
  ! high verbosity... tests if they are working.

#ifdef NEW_DIAG  
  if (proc0) override_screen_printout_options = should_print(3)
  call broadcast(override_screen_printout_options)
#endif

  test_driver_flag = .true.
  functional_test_flag = .true.

  call announce_module_test("gs2_diagnostics_new")

  call run_gs2(mp_comm)


  !call announce_test('results')
  !call process_test(test_function(), 'results')

#ifdef NEW_DIAG
    if (proc0) then
      if (use_old_diagnostics) then 
        open(120349, file='averages.dat', status="replace", action="write")
        write(120349, *) qflux_avg
        write(120349, *) pflux_avg
        write(120349, *) vflux_avg
        write(120349, *) heat_avg
        write(120349, *) diffusivity()
        close(120349)
      else
        if (nstep.eq.200) then
          allocate(qfluxav(nspec), pfluxav(nspec), heatav(nspec), vfluxav(nspec))
          open(120349, file='averages.dat')
          read(120349, *) qfluxav
          read(120349, *) pfluxav
          read(120349, *) vfluxav
          read(120349, *) heatav
          read(120349, *) diff
          close(120349)
          call announce_test("average heat flux")
          call process_test( &
            agrees_with(gnostics%current_results%species_heat_flux_avg, qfluxav, eps), &
            "average heat flux")
          call announce_test("average momentum flux")
          call process_test( &
            agrees_with(gnostics%current_results%species_momentum_flux_avg, vfluxav, eps), &
            "average momentum flux")
          call announce_test("average particle flux")
          call process_test( &
            agrees_with(gnostics%current_results%species_particle_flux_avg, pfluxav, eps), &
            "average particle flux")
          call announce_test("diffusivity")
          call process_test( &
            agrees_with(gnostics%current_results%diffusivity, diff, eps), &
            "diffusivity")
        else
          call announce_test("Size of t array")
        end if
      end if
    end if
#endif




#ifdef NEW_DIAG
    if (proc0 .and. .not. use_old_diagnostics) then
      if (nstep==200) then 
        do i = 1,n_vars
          call announce_test("value of "//trim(new_variables(i)))
          call process_test(test_variable(trim(variables(i)), trim(new_variables(i)), &
            trim(n_lines(i))), &
            "value of "//trim(new_variables(i)))
        end do
      else if (gnostics%appending) then 
        do i = 1,n_vars
          ! omega and omega_average won't work because the 
          ! history is not stored in the restart file
          if (i.eq.14.or.i.eq.15) cycle
          call announce_test("value of "//trim(new_variables(i)))
          call process_test(&
            test_variable(trim(new_variables(i)), trim(new_variables(i)), &
            trim(n_lines(i))), &
            "value of "//trim(new_variables(i)))
        end do
      end if
      if (nstep==200) then 
        ! Note that heat and heat2 no longer pass because the comparison
        ! is no longer between files output from the same run but between
        ! files output from two successive runs. This causes changes greater
        ! than the precision of the heat and heat2 files. EGH
        do i = 3,n_file_names
          call announce_test("content of "//trim(file_names(i)))
          call process_test(test_file(trim(file_names(i)), trim(n_lines_files(i))), &
            "content of "//trim(file_names(i)))
        end do
      end if
    end if
#endif

    if (use_old_diagnostics) then
      call finish_gs2_diagnostics(ilast_step)

#ifdef NEW_DIAG
    else
      call finish_gs2_diagnostics_new
#endif
    end if
    call finish_gs2

   call close_module_test("gs2_diagnostics_new")

    call finish_mp

contains
  
  function test_variable(var_name, new_var_name, n_lines)
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: gnostics
#endif
    use unit_tests, only: should_print
    character(*), intent(in) :: var_name, new_var_name, n_lines
    logical :: test_variable
    character(1000) ::  command 
    
    test_variable=.true.
#ifdef NEW_DIAG
!    command = "if [ ""`ncdump -v "//var_name//" test_gs2_diagnostics_new.out.nc  | tail -n "//n_lines//"`"" = &
!     &  ""`ncdump -v "//new_var_name//" test_gs2_diagnostics_new.cdf | tail -n "//n_lines//" `"" ]; &
!     & then echo ""T"" > tmpdata.dat; fi"
!    write(command,'("if [ ",A,"$(./getncdat ",A," ",A,")",A," &
!      & = ",A,"$(./getncdat ",A," ",A,")",A," ] ; &
!      & then echo ",A," > tmpdata.dat ; fi")') &
!         '"',"test_gs2_diagnostics_new.out.nc",var_name,'"',&
!         '"',"test_gs2_diagnostics_new.cdf",new_var_name,'"',"'T'"
    
    command = ''
    if (gnostics%append_old) then 
      ! Here we are comparing the new diagnostics 200 step run with 
      ! a 100-step run with an appended 100 step run.
      write(command,'("./compare ",A,A," ",A,A," ",A,A," ",A,A," ",A,A)') &
         '"',"test_gs2_diagnostics_new.out.nc",var_name,'"',&
         '"',"test_gs2_diagnostics_new_append.out.nc",new_var_name,'"'
    else 
      write(command,'("./compare ",A,A," ",A,A," ",A,A," ",A,A," ",A,A)') &
         '"',"old_diagnostics/test_gs2_diagnostics_new.out.nc",var_name,'"',&
         '"',"test_gs2_diagnostics_new.out.nc",new_var_name,'"'
     end if
    
    if (should_print(3)) write(*,*) trim(command)
    call system(" echo ""F"" > tmpdata.dat")
    call system(command)
    test_variable = .true.
    open(120349, file='tmpdata.dat')
    read(120349, '(L)') test_variable
    close(120349)
#endif
  end function test_variable
  function test_file(file_name, n_lines)
    use unit_tests, only: should_print
    character(*), intent(in) :: file_name, n_lines
    logical :: test_file
    character(1000) ::  command 
    
    test_file=.true.
#ifdef NEW_DIAG
    command = "if [ ""`cat old_diagnostics/test_gs2_diagnostics_new."//file_name//"  &
     & | tail -n "//n_lines//"`"" = &
     &  ""`cat test_gs2_diagnostics_new."//file_name//"   | tail -n "//n_lines//"`"" ]; &
     & then echo ""T"" > tmpdata.dat; fi"
    
    if (should_print(3)) write(*,*) trim(command)
    call system(" echo ""F"" > tmpdata.dat")
    call system(command)
    !test_= .true.
    open(120349, file='tmpdata.dat')
    read(120349, '(L)') test_file
    close(120349)
#endif
  end function test_file


end program test_gs2_diagnostics_new
