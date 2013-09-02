!> A test program to run a linear slab itg benchmark
!! at low resolution and check the growth rate

module checks_mod
  use functional_tests
  public checks
  contains
    function checks()
      logical :: checks
      ! We only expect this low res test to get the growth rate to 3%
      checks =  check_growth_rate((/0.164/), 0.3)
    end function checks
end module checks_mod

program slab_itg_s_hat_0_low_res
  use functional_tests
  use checks_mod


  call test_gs2('Linear slab ITG with zero magnetic shear', checks)


end program slab_itg_s_hat_0_low_res
