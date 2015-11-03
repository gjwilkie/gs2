module leq
  implicit none

  integer :: nr, nt
  private
  
  real, allocatable, dimension (:)     :: eqpsi, fp, beta, pressure
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi
  real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm
  real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbtcart
  real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbtbish

  real :: beta_0
  
  type :: flux_surface
     real :: R_center, R_geo, k, kp, d, dp, r, dr, delp, q, shat, pp, a, ap
     integer :: nt
  end type flux_surface

  type (flux_surface) :: surf

  public :: leq_init, leqin, gradient, eqitem, bgradient

  public :: invR, Rpos, Zpos, diameter, btori, dbtori,  qfun, pfun, &
       dpfun, betafun, psi, rcenter, dpdrhofun

contains

  subroutine leqin(R0, Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap, nt_used)
        
    real :: R0, Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap
    integer :: nt_used

    surf%R_center = R0
    surf%R_geo = Ra
    surf%delp = s
    surf%k = k
    surf%kp = kp
    surf%q = qq
    surf%shat = qs
    surf%d = d
    surf%dp = dp
    surf%a = a
    surf%ap = ap
    surf%r = r
    surf%dr = dr
    surf%pp = 0.
    
    beta_0 = 1.

    nr = 3
    nt = nt_used
    if(.not.allocated(beta)) call alloc_arrays(3, nt)
    surf%nt = nt
    call leq_init

  end subroutine leqin

  subroutine alloc_arrays(nr, nt)

    integer :: nr, nt

    allocate(eqpsi(nr), fp(nr), beta(nr), pressure(nr))
    allocate(R_psi(nr, nt), Z_psi(nr, nt))
    allocate(drm(nr, nt, 2), dzm(nr, nt, 2), dbtm(nr, nt, 2), &
         dpm(nr, nt, 2), dtm(nr, nt, 2))
    allocate(dpcart(nr, nt, 2), dtcart(nr, nt, 2), dbtcart(nr, nt, 2))
    allocate(dpbish(nr, nt, 2), dtbish(nr, nt, 2), dbtbish(nr, nt, 2))

  end subroutine alloc_arrays

  subroutine leq_init

    implicit none
    real, dimension(nr, nt) :: eqpsi1, eqth, eqbtor

    real dr(3)
    real pi, t, r
    integer i, j
   
    pi=2*acos(0.)
    dr(1) = -surf%dr
    dr(2) = 0.
    dr(3) = surf%dr
    
    do j=1,nt
       do i=1,nr
          r = surf%r + dr(i)
          t = (j-1)*pi/real(nt-1)
          R_psi(i,j) = Rpos(r, t) 
          Z_psi(i,j) = Zpos(r, t)
          eqth(i,j) = t
          eqpsi1(i,j) = 1 + dr(i)
          eqbtor(i,j) = surf%r_geo/R_psi(i,j)
       enddo
    enddo

    do i=1,nr
       pressure(i) = -dr(i)
    enddo

    eqpsi(:) = eqpsi1(:,1)

    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
    call derm(eqbtor, dbtm, 'E')
    call derm(eqpsi1, dpm,  'E')
    
! diagnostics
!      do j=1,nt
!         do i=1,nr
!            write(*,*) i,j
!            write(*,100) drm(i,j,1),drm(i,j,2),R_psi(i,j)
!            write(*,101) dzm(i,j,1),dzm(i,j,2),Z_psi(i,j)
!            write(*,102) dtm(i,j,1),dtm(i,j,2),eqth(i,j)
!         enddo
!      enddo
! 100  format('(gr R)1 ',g10.4,' (gr R)2 ',g10.4,' R ',g10.4)
! 101  format('(gr Z)1 ',g10.4,' (gr Z)2 ',g10.4,' Z ',g10.4)
! 102  format('(gr t)1 ',g10.4,' (gr t)2 ',g10.4,' t ',g10.4)
!      write(*,*) nr, nt
!      stop

! below is actually grad(rho) instead of grad(psi),
! and 'cartesian' refers to (R,Z) coordinates -- MAB
! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

! grad(BT) in cartesian form
    call eqdcart(dbtm, dbtcart)
! grad(BT) in Bishop form
    call eqdbish(dbtcart, dbtbish)

! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

  end subroutine leq_init

  subroutine derm(f, dfm, char)

    implicit none
    integer i, j
    character(1) :: char
    real f(:,:), dfm(:,:,:), pi

    pi = 2*acos(0.)
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))
   
! assume up-down symmetry for now:
 
    select case (char)
    case ('E') 
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))
       
       j=nt      
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
    case ('O')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))
       
       j=nt      
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j=1
       dfm(:,j,2) = f(:,j+1)
       
       j=nt      
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo
    
    do j=2,nt-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo
    
  end subroutine derm

  subroutine gradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character(1) char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real tmp(2), aa(1), daa(1), rp, rpt(1)
    real, dimension(nr,nt,2) :: dcart
    integer i
    
    select case(char)
    case('D')  ! diagnostic 
       dcart = dbtcart
    case('P') 
       dcart = dpcart
    case('R') 
       dcart = dpcart  ! dpcart is correct for 'R'
    case('T')
       dcart = dtcart
    end select
    
    do i=-nth_used,-1
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       if(char == 'T') then
          grad(i,1)=-tmp(1)
          grad(i,2)=-tmp(2)
       else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       endif
    enddo

    do i=0,nth_used
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1)*0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1)*0.5*beta_0
       enddo
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character(1) char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real :: aa(1), daa(1), rp, rpt(1)
    real, dimension(nr,nt,2) ::  dbish
    integer i

    select case(char)
    case('D')  ! diagnostic
       dbish = dbtbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       call eqitem(rgrid(i), theta(i), dbish(:,:,1), grad(i,1), 'R')
       call eqitem(rgrid(i), theta(i), dbish(:,:,2), grad(i,2), 'Z')
    enddo

    if (char == 'T') then
       where (theta(-nth_used:nth_used) < 0.0)
          grad(-nth_used:nth_used,1) = -grad(-nth_used:nth_used,1)
          grad(-nth_used:nth_used,2) = -grad(-nth_used:nth_used,2)
       end where
    end if

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif

  end subroutine bgradient

  subroutine eqitem(r, theta_in, f, fstar, char)
      
    integer :: j, istar, jstar
    character(1) :: char
    real :: r, thet, fstar, tp, tps, theta_in
    real :: st, dt, pi
    real, dimension(:,:) :: f
    real, dimension(size(f,2)) :: mtheta
    
    pi = 2.*acos(0.)
! find r on psi mesh
    
    istar = 2

! Now do theta direction

    thet = mod2pi(theta_in)

! assume up-down symmetry

    tp=abs(thet)
    tps=1.
    if(char == 'Z' .and. thet /= 0.) tps=thet/abs(thet)
        
! get thet on theta mesh

    mtheta = (/ ( real(j-1)*pi/real(nt-1), j=1,nt) /)
  
! note that theta(1)=0 for local_eq theta 

    jstar=-1    
    do j=1,nt
       if(tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar=j-1
          exit
       endif
       if(jstar /= -1) write(*,*) 'exit error j'
    enddo
      
! treat theta = pi separately
  
    if(jstar == -1) then
       jstar=nt-1
       dt=mtheta(jstar+1)-mtheta(jstar)
       st=0.
    endif

! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    )  * st &
         +f(istar    , jstar + 1)  * dt 
    fstar=fstar*tps &
         /(mtheta(jstar+1)-mtheta(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) eqpsi(istar),eqpsi(istar+1)
!     write(*,*) mtheta(jstar),mtheta(jstar+1)
      
  end subroutine eqitem

  subroutine eqdcart(dfm, dfcart)
      
    implicit none

    real, dimension (:,:,:), intent(in)  :: dfm
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (size(dfm,1),size(dfm,2)) :: denom
    integer :: i, j
      
    denom(:,:) = drm(:,:,1)*dzm(:,:,2) - drm(:,:,2)*dzm(:,:,1)

    dfcart = 0.
    
    dfcart(:,:,1) =   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)
    
    do j=1,nt
       do i=2,nr
          dfcart(i,j,:)=dfcart(i,j,:)/denom(i,j)
       enddo
    enddo
    
  end subroutine eqdcart

  subroutine eqdbish(dcart, dbish)
    
    implicit none
    real, dimension(:, :, :), intent (in) :: dcart
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j

    denom(:,:) = sqrt(dpcart(:,:,1)**2 + dpcart(:,:,2)**2)

    dbish(:,:,1) = dcart(:,:,1)*dpcart(:,:,1) + dcart(:,:,2)*dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1)*dpcart(:,:,2) + dcart(:,:,2)*dpcart(:,:,1)
    
    do j=1,nt
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
      enddo
    enddo

  end subroutine eqdbish

  function invR (r, theta)
   
    real, intent (in) :: r, theta
    real :: invR
    
    invR=1./Rpos(r, theta)
    
  end function invR

  function rcenter(rp)

    real, intent(in) :: rp
    real :: rcenter

    rcenter = surf%R_center
    
  end function rcenter

  function Rpos (r, theta)
   
    use constants, only: pi

    real, intent (in) :: r, theta
    real :: Rpos
    real :: g, gp, dr
    
    dr = r - surf%r

! For Y Xiao: 
!    g = surf%delp/surf%r + surf%d * sin(theta)**2
!    Rpos = surf%R_center*(1.+r*(cos(theta)-g)-g*dr)

    g = cos(theta + surf%d * sin(theta-0.5*pi*surf%a))
    gp = -sin(theta + surf%d * sin(theta-0.5*pi*surf%a))*(surf%dp*sin(theta-0.5*surf%a)-surf%d*0.5*pi*surf%ap*cos(theta-0.5*pi*surf%a))
    
    Rpos = surf%R_center + surf%delp*dr + g*surf%r + (g+surf%r*gp)*dr
    
  end function Rpos

  function Zpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: Zpos, dr
    
    dr = r - surf%r
    Zpos = surf%k*sin(theta)*surf%r + (surf%r*surf%kp + surf%k)*sin(theta)*dr
    
  end function Zpos

  function psi (r, theta)
   
    real, intent (in) :: r, theta
    real :: psi

    psi = r - surf%r
    
  end function psi

  function mod2pi (theta)
    
    real, intent(in) :: theta
    real :: pi, th, mod2pi
    real, parameter :: theta_tol = 1.e-6
    logical :: out
    
    pi=2.*acos(0.)
    
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
       return
    endif
    
    if(theta - theta_tol <= pi .and. theta >= -pi) then
       mod2pi = pi
       return
    endif

    if(theta <= pi .and. theta + theta_tol >= -pi) then
       mod2pi = -pi
       return
    endif

    th=theta
    out=.true.
    do while(out)
       if(th > pi) th = th - 2.*pi
       if(th <-pi) th = th + 2.*pi
       if(th <= pi .and. th >= -pi) out=.false.
    enddo
    mod2pi=th
    
  end function mod2pi
   
  function diameter (rp)
  
    real :: rp, diameter

    diameter = 2.*rp

  end function diameter

  function dbtori (pbar)
    real :: pbar, dbtori
    dbtori = 1.
  end function dbtori

  function btori (pbar)
    real :: pbar, btori
    btori = surf%r_geo
  end function btori

  function qfun (pbar)
    real :: pbar, qfun
    qfun = surf%q
  end function qfun

  function pfun (pbar)
    real :: pbar, pfun
    pfun = 0.5*beta_0
  end function pfun
  
  function dpfun (pbar)  
    real :: pbar, dpfun    

       dpfun = -1.

  end function dpfun

  function dpdrhofun(rho)

    real :: rho, dpdrhofun

    dpdrhofun = surf%pp

  end function dpdrhofun
  
  function betafun (pbar)  
    real :: pbar, betafun
    betafun = beta_0
  end function betafun

end module leq
