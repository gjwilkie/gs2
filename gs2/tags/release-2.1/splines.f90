module splines
  implicit none

  type :: spline
     integer :: n
     real, dimension (:), pointer :: x, y, y2
  end type spline

  type :: periodic_spline
     integer :: n
     real :: period
     real, dimension (:), pointer :: x, y, y2
  end type periodic_spline

contains

  subroutine new_spline (n, x, y, spl)
    implicit none
    integer, intent (in) :: n
    real, dimension (n), intent (in) :: x, y
    type (spline), intent (out) :: spl
    real, dimension (n) :: temp
    integer :: ierr
    interface
       subroutine fitp_curv1 (n,x,y,slp1,slpn,islpsw,y2,temp,sigma,ierr)
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y
         real, intent (in) :: slp1, slpn
         integer, intent (in) :: islpsw
         real, dimension (n), intent (out) :: y2
         real, dimension (n), intent (out) :: temp
         real, intent (in) :: sigma
         integer, intent (out) :: ierr
       end subroutine fitp_curv1
    end interface

    spl%n = n
    allocate (spl%x(n),spl%y(n))
    spl%x = x
    spl%y = y
    allocate (spl%y2(n))
    call fitp_curv1 (n, x, y, 0.0, 0.0, 3, spl%y2, temp, 1.0, ierr)
  end subroutine new_spline

  subroutine new_periodic_spline (n, x, y, period, spl)
    implicit none
    integer, intent (in) :: n
    real, dimension (n), intent (in) :: x, y
    real, intent (in) :: period
    type (periodic_spline), intent (out) :: spl
    real, dimension (2*n) :: temp
    integer :: ierr
    interface
       subroutine fitp_curvp1 (n,x,y,p,yp,temp,sigma,ierr)
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y
         real, intent (in) :: p
         real, dimension (n), intent (out) :: yp
         real, dimension (2*n), intent (out) :: temp
         real, intent (in) :: sigma
         integer, intent (out) :: ierr
       end subroutine fitp_curvp1
    end interface

    spl%n = n
    spl%period = period
    allocate (spl%x(n),spl%y(n))
    spl%x = x
    spl%y = y
    allocate (spl%y2(n))
    call fitp_curvp1 (n,x,y,period,spl%y2,temp,1.0,ierr)
  end subroutine new_periodic_spline

  subroutine delete_spline (spl)
    implicit none
    type (spline), intent (in out) :: spl
    spl%n = 0
    deallocate (spl%x,spl%y)
    nullify (spl%x)
    nullify (spl%y)
    deallocate (spl%y2)
    nullify (spl%y2)
  end subroutine delete_spline

  subroutine delete_periodic_spline (spl)
    implicit none
    type (periodic_spline), intent (in out) :: spl
    spl%n = 0
    spl%period = 0.0
    deallocate (spl%x,spl%y)
    nullify (spl%x)
    nullify (spl%y)
    deallocate (spl%y2)
    nullify (spl%y2)
  end subroutine delete_periodic_spline

  function splint (x, spl)
    implicit none
    real, intent (in) :: x
    type (spline), intent (in) :: spl
    integer i
    real :: splint
    interface
       function fitp_curv2 (x0, n, x, y, y2, sigma)
         real, intent (in) :: x0
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y, y2
         real, intent (in) :: sigma
         real :: fitp_curv2
       end function fitp_curv2
    end interface
!    do i=1,spl%n
!       write(*,*) spl%x(i), spl%y(i), spl%y2(i)
!    enddo
    splint = fitp_curv2 (x, spl%n, spl%x, spl%y, spl%y2, 1.0)
  end function splint

  function periodic_splint (x, spl)
    implicit none
    real, intent (in) :: x
    type (periodic_spline), intent (in) :: spl
    real :: periodic_splint
    interface
       function fitp_curvp2 (x0, n, x, y, p, y2, sigma)
         real, intent (in) :: x0
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y, y2
         real, intent (in) :: p, sigma
         real :: fitp_curvp2
       end function fitp_curvp2
    end interface
    periodic_splint = fitp_curvp2 &
         (x, spl%n, spl%x, spl%y, spl%period, spl%y2, 1.0)
  end function periodic_splint

  function dsplint (x, spl)
    implicit none
    real, intent (in) :: x
    type (spline), intent (in) :: spl
    real :: dsplint
    interface
       function fitp_curvd (x0, n, x, y, y2, sigma)
         real, intent (in) :: x0
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y, y2
         real, intent (in) :: sigma
         real :: fitp_curvd
       end function fitp_curvd
    end interface
    dsplint = fitp_curvd (x, spl%n, spl%x, spl%y, spl%y2, 1.0)
  end function dsplint

  function splintint (x0, x1, spl)
    implicit none
    real, intent (in) :: x0, x1
    type (spline), intent (in) :: spl
    real :: splintint
    interface
       function fitp_curvi (x0, x1, n, x, y, y2, sigma)
         real, intent (in) :: x0, x1
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y, y2
         real, intent (in) :: sigma
         real :: fitp_curvi
       end function fitp_curvi
    end interface
    splintint = fitp_curvi (x0,x1,spl%n,spl%x,spl%y,spl%y2,1.0)
  end function splintint

  function periodic_splintint (x0, x1, spl)
    implicit none
    real, intent (in) :: x0, x1
    type (periodic_spline), intent (in) :: spl
    real :: periodic_splintint
    interface
       function fitp_curvpi (x0, x1, n, x, y, p, y2, sigma)
         real, intent (in) :: x0, x1
         integer, intent (in) :: n
         real, dimension (n), intent (in) :: x, y, y2
         real, intent (in) :: p, sigma
         real :: fitp_curvpi
       end function fitp_curvpi
    end interface
    periodic_splintint = fitp_curvpi &
         (x0,x1,spl%n,spl%x,spl%y,spl%period,spl%y2, 1.0)
  end function periodic_splintint

end module splines
