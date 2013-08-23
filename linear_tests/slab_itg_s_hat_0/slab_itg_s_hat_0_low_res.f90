!> A test program to run a linear slab itg benchmark
!! at low resolution and check the growth rate

program slab_itg_s_hat_0_low_res
  use functional_tests


  call test_gs2('Linear slab ITG with zero magnetic shear', checks)

contains
  function checks()
    logical :: checks
    ! We only expect this low res test to get the growth rate to 10%
    checks =  check_growth_rate((/0.164/), 0.1)
  end function checks

end program slab_itg_s_hat_0_low_res
