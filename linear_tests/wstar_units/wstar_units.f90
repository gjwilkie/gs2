
module checks_mod
  use functional_tests
  public checks
  contains
    function checks()
      logical :: checks
      ! We only expect this low res test to get the growth rates equal to 10%
      checks =  check_growth_rates_equal_in_list(0.1)
    end function checks
end module checks_mod

!> A test program to run a linear test case in list mode 
!! with and without wstar_units and compare the growth
!! rates in the two cases

program wstar_units
  use functional_tests
  use checks_mod


  call test_gs2('List mode and wstar units', checks)


end program wstar_units
