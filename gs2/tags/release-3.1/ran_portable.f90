module ran
! D. Ernst: change Cray/SGI ranf to standard f90 calls to random_number
contains

  function ranf()
    real :: random

    call random_number(random)
    ranf = random

  end function ranf

end module ran

