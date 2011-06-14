# include "define.inc"

module ran

  ! <doc>
  !  A wrapper module for random number generator.
  !  Thie module provides real function ranf
  !  using intrinsic random_number/random_seed or
  !  Mersenne Twister 19937 (see mt19937.f90)
  ! </doc>

  implicit none

  private

  public :: ranf

contains

  function ranf (seed)
    
    ! <doc>
    !  returns a uniform deviate in [0., 1.)
    !  The generator is initialized with the given seed if exists,
    !  otherwise uses the default seed.
    ! </doc>
    
# if RANDOM == _RANMT_
    use mt19937, only: sgrnd, grnd
# endif
    integer, intent(in), optional :: seed
    real :: ranf
    integer :: l
    integer, allocatable :: seed_in(:)
# if RANDOM == _RANMT_    
    if (present(seed)) call sgrnd(seed)
    ranf = grnd()
# else
    if (present(seed)) then
       call random_seed(size=l)
       allocate(seed_in(l))
       seed_in(:)=seed
       call random_seed(put=seed_in)
    endif
    call random_number(ranf)
# endif

  end function ranf

end module ran
