!> A test program to run a linear cyclone itg benchmark
!! at low resolution and check the growth rate

module checks_mod
  use functional_tests
  public checks
  contains
    function checks()
      logical :: checks
      ! We only expect this low res test to get the growth rate to 3%
      checks =  check_growth_rate((/0.1687/), 0.03)
    end function checks
end module checks_mod

program cyclone_itg_low_res
  use functional_tests
  use checks_mod
  use command_line
  use fields_local, only: fields_local_functional
  use mp, only: init_mp, proc0, finish_mp

  character(len=8) :: field_type
  integer :: length, ierr

  ! For this linear test there is an additional command line argument (see
  ! Makefile) which tells the test program which input file it has been passed.
  ! If it has been passed the input file which specifies local fields, it sees
  ! if fields_local is functional; if not, it exits.
  ! This is because we don't want the unit test to cause an error when the tests
  ! have been built with a compiler that doesn't support fields_local.
  ! 
  ! If you want to see the minimal linear test look at slab_itg_low_res
  call init_mp
  call cl_getarg(2, field_type, length, ierr)
  if (trim(field_type) == "local" .and. .not. fields_local_functional()) then
    if (proc0) then
      write (*,*)  
      write (*,*) "WARNING: fields_local is non-functional... skipping test&
       & cyclone_itg_low_res with fields local. "
      write (*,*)
    end if
    call finish_mp
    stop
  end if


  call test_gs2('Linear CBC ITG in single mode with field algorithm '//trim(field_type), checks)


end program cyclone_itg_low_res
