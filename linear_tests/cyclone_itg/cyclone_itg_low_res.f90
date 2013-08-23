!> A test program to run a linear cyclone itg benchmark
!! at low resolution and check the growth rate

program cyclone_itg_low_res
  use functional_tests


  call test_gs2('Linear cyclone base case ITG in single mode', checks)

contains
  function checks()
    logical :: checks
    ! We only expect this low res test to get the growth rate to 1%
    checks =  check_growth_rate((/0.1687/), 0.01)
  end function checks

end program cyclone_itg_low_res
