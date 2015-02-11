!> A test program to run a linear cyclone itg benchmark
!! at low resolution and check the growth rate

module checks_mod
  use functional_tests
  public checks
  contains
    function checks()
      use gs2_time, only: code_dt_cfl
      use unit_tests, only: agrees_with
      implicit none
      logical :: checks
      !Check if the cfl condition matches expected value with 3%
      !Checking against a predefined "correct" value, this may not
      !be a fatal error if agreement not found.
      checks =  agrees_with(code_dt_cfl,0.11400483899841848,0.03)
    end function checks
end module checks_mod

program cyclone_itg_low_res
  use functional_tests, only: test_gs2
  use checks_mod, only: checks
  use fields_local, only: fields_local_functional
  use mp, only: init_mp, proc0, finish_mp
  use gs2_time, only: code_dt_cfl
  implicit none
  character(len=8) :: field_type
  integer :: length, ierr

  ! If you want to see the minimal linear test look at slab_itg_low_res
  call init_mp
  call test_gs2('Linear CBC ITG in single mode with field algorithm '//trim(field_type), checks)
  if(proc0) print*,code_dt_cfl," ",0.11400483899841848
  call finish_mp
end program cyclone_itg_low_res
