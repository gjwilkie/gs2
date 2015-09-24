module regression

  implicit none
  type :: reg_type
    integer :: n = 0
    real :: xy = 0.
    real :: xsum = 0.
    real :: ysum = 0.
    real :: x2sum = 0.
    real :: a, b
  end type reg_type
  
  interface regress
     module procedure regress_0, regress_1, regress_2
  end interface

contains

  subroutine regress_0 (reg, x, y, a, b)
    type (reg_type), target :: reg
    real, intent (in) :: x, y
    real, intent (out), optional :: a, b

    integer, pointer :: n
    real, pointer :: xy, xsum, ysum, x2sum
    
    n     => reg % n 
    xy    => reg % xy
    xsum  => reg % xsum
    ysum  => reg % ysum
    x2sum => reg % x2sum

    n = n + 1
    xy = xy + x*y
    xsum = xsum + x
    ysum = ysum + y
    x2sum = x2sum + x*x
    
    if (n > 1) then
       reg%b = (xy-xsum*ysum/real(n))/(x2sum-xsum*xsum/real(n))
       reg%a = ysum/real(n) - reg%b*xsum/real(n)
    else
       reg%b = 0.
       reg%a = 0.
    end if
    
    if (present(a)) then
       b=reg%b
       a=reg%a
    end if

  end subroutine regress_0

  subroutine regress_1 (reg, x, y, a, b)
    type (reg_type), target :: reg
    real, dimension(:), intent (in) :: x, y
    real, dimension(:), intent (out), optional :: a, b

    integer, pointer :: n
    real, pointer :: xy, xsum, ysum, x2sum
    
    n     => reg % n 
    xy    => reg % xy
    xsum  => reg % xsum
    ysum  => reg % ysum
    x2sum => reg % x2sum

    n = n + 1
    xy = xy + sum(x*y)
    xsum = xsum + sum(x)
    ysum = ysum + sum(y)
    x2sum = x2sum + sum(x*x)

    if (n>1) then
       reg%b = (xy-xsum*ysum/real(n))/(x2sum-xsum*xsum/real(n))
       reg%a = ysum/real(n) - reg%b*xsum/real(n)
    else
       reg%b = 0.
       reg%a = 0.
    end if
    
    if (present(a)) then
       b=reg%b
       a=reg%a
    end if

  end subroutine regress_1

  subroutine regress_2 (reg, x, y, a, b)
    type (reg_type), dimension(:), target :: reg
    real, dimension(:), intent (in) :: x, y
    real, dimension(:), intent (out), optional :: a, b

    integer, dimension(:), pointer :: n
    real, dimension(:), pointer :: xy, xsum, ysum, x2sum
    
    n     => reg % n 
    xy    => reg % xy
    xsum  => reg % xsum
    ysum  => reg % ysum
    x2sum => reg % x2sum

    n = n + 1
    xy = xy + x*y
    xsum = xsum + x
    ysum = ysum + y
    x2sum = x2sum + x*x

    if (n(1)>1) then
       reg%b = (xy-xsum*ysum/real(n))/(x2sum-xsum*xsum/real(n))
       reg%a = ysum/real(n) - reg%b*xsum/real(n)
    else
       reg%b = 0.
       reg%a = 0.
    end if
    
    if (present(a)) then
       b=reg%b
       a=reg%a
    end if

  end subroutine regress_2
  
  elemental function yp (reg, x)
    real :: yp
    type (reg_type), intent (in) :: reg
    real, intent (in) :: x

    yp = reg%a + reg%b*x

  end function yp

  elemental function xp (reg, y)
    real :: xp
    type (reg_type), intent (in) :: reg
    real, intent (in) :: y

    if (reg%b == 0.) then
       xp = 100.
    else
       xp = (y-reg%a)/reg%b
    end if

  end function xp

end module regression
