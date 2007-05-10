module constants

  implicit none


! Symbolic names for kind type of single and double-precision reals:
! (with at least 6 and 12 digits of accuracy)

  integer, parameter :: sp = selected_real_kind(6)
  integer, parameter :: dp = selected_real_kind(12)

! Symbolic names for kind type of single and double-precision complex:

  integer, parameter :: spc = kind((1.0_sp,1.0_sp))
  integer, parameter :: dpc = kind((1.0_dp,1.0_dp))

!  complex(dp), parameter :: ii = (0._dp, 1._dp)
!  real(dp), parameter :: pi=3.141592653589793238_dp

  complex, parameter :: zi = ( 0.0 , 1.0 )
  real, parameter :: pi = 3.1415926535897931

! Note: we will use dp="double precision" for almost everything.
!
! The fortran-90 "kind" types is kind of awkward.  But the old trick of
! using a "-r8" compiler switch to promote all real variables to 64 bits 
! doesn't work on some fortran 90 compilers, and so the above use of 
! the standard fortran-90 routine selected_real_kind is more portable.
!
! It may not be a good idea to mimic "-r8" by making sp to be identical
! to dp, or to write single and double-precision versions of 
! generic subroutines, since on the Cray computers both single and
! "double" precision are 64 bits, and the compiler will complain that
! it can't distinguish the two specific subroutines.  In some cases,
! the cray compiler may be able to distinguish between two real "kinds"
! for the purposes of distinguishing overloaded procedure names,
! even though the two real kinds map to the same precision (64 bits).
!
! If this ever does become a problem, then you may be able to get around it by
! commenting out the double precision function names from the list of 
! overloaded procedures (i.e., the "module procedure" statements).
!
end module constants


