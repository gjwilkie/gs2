# include "define.inc"

!
! special function wrapper routine written by Tomo Tatsuno (5/7/08)
! only has Bessel function at the moment
!
! Unfortunately we do not support elemental feature for ifort...
! XL-fortran does not seem to have intrinsic Bessel function
!
! To do: avoid explicit specification of kind and use kind_rs or kind_rd
!        support of other compilers such as absoft, lahay etc...
!

module spfunc

  implicit none

  public :: j0, j1

  private

# if (FCOMPILER == _G95_ || FCOMPILER == _GFORTRAN_ || FCOMPILER == _INTEL_ \
  || FCOMPILER == _PATHSCALE_)

  interface j0
     module procedure sj0, dj0
  end interface
  interface j1
     module procedure sj1, dj1
  end interface

# elif FCOMPILER == _PGI_

  interface j0
     elemental function besj0(x)
       real (4), intent (in) :: x
       real (4) :: besj0
     end function besj0
     elemental function dbesj0(x)
       real (8), intent (in) :: x
       real (8) :: dbesj0
     end function dbesj0
  end interface
  ! j1 is below

# endif

contains

# if (FCOMPILER == _G95_ || FCOMPILER == _GFORTRAN_ || FCOMPILER == _INTEL_ \
  || FCOMPILER == _PATHSCALE_)

# if FCOMPILER == _INTEL_
  function sj0 (x)
    use ifport, only: besj0
# else
  elemental function sj0 (x)
# endif
    real (4), intent (in) :: x
    real (4) :: sj0
    sj0 = besj0(x)
  end function sj0

# if FCOMPILER == _INTEL_
  function dj0 (x)
    use ifport, only: dbesj0
# else
  elemental function dj0 (x)
# endif
    real (8), intent (in) :: x
    real (8) :: dj0
    dj0 = dbesj0(x)
  end function dj0

# if FCOMPILER == _INTEL_
  function sj1 (x)
    use ifport, only: besj1
# else
  elemental function sj1 (x)
# endif
    real (4), intent (in) :: x
    real (4) :: sj1
    if (x == 0.0) then
       sj1 = 0.5
    else
       sj1 = besj1(x) / x
    end if
  end function sj1

# if FCOMPILER == _INTEL_
  function dj1 (x)
    use ifport, only: dbesj1
# else
  elemental function dj1 (x)
# endif
    real (8), intent (in) :: x
    real (8) :: dj1
    if (x == 0.0) then
       dj1 = 0.5
    else
       dj1 = dbesj1(x) / x
    end if
  end function dj1

# elif FCOMPILER == _PGI_

  elemental function j1 (x)
    real, intent (in) :: x
    real :: j1
    interface besj1
       elemental function besj1(x)
         real (4), intent (in) :: x
         real (4) :: besj1
       end function besj1
       elemental function dbesj1(x)
         real (8), intent (in) :: x
         real (8) :: dbesj1
       end function dbesj1
    end interface
    if (x == 0.0) then
       j1 = 0.5
    else
       j1 = besj1(x) / x
    end if
  end function j1

# else

  elemental function j0 (x)
! A&S, p. 369, 9.4
    implicit none
    real, intent (in) :: x
    real :: j0
    real, parameter, dimension (7) :: a = &
         (/ 1.0000000, -2.2499997, 1.2656208, -0.3163866, &
            0.0444479, -0.0039444, 0.0002100 /)
    real, parameter, dimension (7) :: b = &
         (/  0.79788456, -0.00000770, -0.00552740, -0.00009512, &
             0.00137237, -0.00072805,  0.00014476 /)
    real, parameter, dimension (7) :: c = &
         (/ -0.78539816, -0.04166397, -0.00003954,  0.00262573, &
            -0.00054125, -0.00029333,  0.00013558 /)
    real :: y

    if (x <= 3.0) then
       y = (x/3.0)**2
       j0 = a(1)+y*(a(2)+y*(a(3)+y*(a(4)+y*(a(5)+y*(a(6)+y*a(7))))))
    else
       y = 3.0/x
       j0 = (b(1)+y*(b(2)+y*(b(3)+y*(b(4)+y*(b(5)+y*(b(6)+y*b(7))))))) &
            *cos(x+c(1)+y*(c(2)+y*(c(3)+y*(c(4)+y*(c(5)+y*(c(6)+y*c(7))))))) &
            /sqrt(x)
    end if
  end function j0

  elemental function j1 (x)
! A&S, p. 370, 9.4 j1 = 1/x J_1(x)
    implicit none
    real, intent (in) :: x
    real :: j1
    real, parameter, dimension (7) :: a = &
         (/  0.50000000, -0.56249985,  0.21093573, -0.03954289, &
             0.00443319, -0.00031761,  0.00001109 /)
    real, parameter, dimension (7) :: b = &
         (/  0.79788456,  0.00000156,  0.01659667,  0.00017105, &
            -0.00249511,  0.00113653,  0.00020033 /)
    real, parameter, dimension (7) :: c = &
         (/ -2.35619449,  0.12499612,  0.00005650,  -0.00637879, &
             0.00074348,  0.00079824, -0.00029166 /)
    real :: y

    if (x <= 3.0) then
       y = (x/3.0)**2
       j1 = a(1)+y*(a(2)+y*(a(3)+y*(a(4)+y*(a(5)+y*(a(6)+y*a(7))))))
    else
       y = 3.0/x
       j1 = (b(1)+y*(b(2)+y*(b(3)+y*(b(4)+y*(b(5)+y*(b(6)+y*b(7))))))) &
            *cos(x+c(1)+y*(c(2)+y*(c(3)+y*(c(4)+y*(c(5)+y*(c(6)+y*c(7))))))) &
            /x**1.5
    end if
  end function j1

# endif

# if FCOMPILER == -1
# this is always false at the momont
  elemental function erf(x)
! A&S, p.299 7.1.28 |epsilon|<=3.e-7
    implicit none
    real, intent(in) :: x
    real :: xerf
    real, parameter, dimension(6) :: a = (/ &
         0.0705230784, 0.0422820123, 0.0092705272, &
         0.0001520143, 0.0002765672, 0.0000430638 /)

    erf = 1.0 - 1.0/(1.0 + &
         x*(a(1) + x*(a(2) + x*(a(3) + x*(a(4) + x*(a(5) + x*(a(6))))))))**16

  end function erf
# endif

end module spfunc
