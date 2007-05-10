module veq
  implicit none

  public :: psi_bar, eqpsi, fp, qsf, pressure, beta, B_psi
  public :: vmomin 
  public :: rho_d, vc, R_0, Z_0
  public :: nr, md

  public :: bgradient, gradient, eqitem, veq_init

!  integer, parameter :: nrx = 41
!  integer, parameter :: ntx = 41
!  integer, parameter :: vmnmx = 50

  public :: invR
  public :: Rpos     
  public :: Zpos     
  public :: diameter, initialize_diameter
  public :: btori,    initialize_btori
  public :: dbtori,   initialize_dbtori
  public :: qfun,     initialize_q
  public :: rcenter,  initialize_rcenter
  public :: pfun,     initialize_pressure
  public :: dpfun,    initialize_dpressure
  public :: betafun,  initialize_beta
  public :: psi,      initialize_psi
  public :: psiofrho

  logical :: init_diameter = .true.
  logical :: init_btori = .true.
  logical :: init_dbtori = .true.
  logical :: init_q = .true.
  logical :: init_pressure = .true.
  logical :: init_dpressure = .true.
  logical :: init_beta = .true.
  logical :: init_psi = .true.
  logical :: init_rcenter = .true.

  private 
  
  integer :: nr, md
  real :: eqpsi_0, eqpsi_a, vr_min, beta_0
  real :: vr_max, vRmag, vz_max, vZmag, vavgrmid
  real :: vB_T, R_0, Z_0, vI_0, vp_0
  real, allocatable, dimension (:) :: rho_d, qsf, fp, psi_bar, beta, pressure, cur, eqpsi, &
	vdvdrho, varea
  real, allocatable, dimension (:,:,:) :: drm, dzm, dpm, dtm, dbm
  real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbcart
  real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbbish
  real, allocatable, dimension (:,:) :: R_psi, Z_psi, B_psi, vc

  logical :: alloced = .false.

!  real, dimension (vmnmx) :: vdvdrho, varea

!  real, dimension (nrx, 6) :: vc

contains

  subroutine vmomin(psi_0, psi_a, rmaj, B_T0, avgrmid, vmeqinit, in_nt, nthg)
    
!     This subroutine reads a VMOM output file containing 
!     the axisymmetric magnetic field geometry on a rectangular 
!     domain defined by the cylindrical coordinates (R,Z).
!
    use vdimstub, only: vdat, vdim
    implicit none
!    include 'vshare.inc'
    real :: psi_0, psi_a, rmaj, B_T0, avgrmid, psi_N, f_N
    integer :: vmeqinit, nthg
    logical :: in_nt
    integer :: i, j

! Need to generalize initialization condition if equilibrium changes

    if(vmeqinit.eq.0) return

! Read the data

    if(.not.in_nt) then
       open(unit=10,file='test.psi',status='old')
       read(10,*) eqpsi_0,eqpsi_a,vB_T,vI_0,vp_0
       read(10,*) vr_min,vr_max,vRmag,vz_max,vZmag,vavgrmid
       read(10,*) nr, md

       if(.not. alloced) then 
          call alloc_module_arrays(nr, md)
          alloced=.true.
       endif

       do i=1,nr
          read(10,*) rho_d(i),psi_bar(i),fp(i),qsf(i),beta(i), &
               pressure(i),cur(i),vdvdrho(i),varea(i),eqpsi(i)
       enddo
       read(10,*) nr, md, i
       read(10,1000) ((B_psi(i,j), i=1,nr), j=1,md)
1000 format(20(1x,g16.10))
       close(10)
       open(unit=10,file='test.vpsi')
       read(10,*) i  ! dummy
       do i=1,nr
          read(10,*) (vc(i,j),j=1,6),rho_d(i)
       enddo
    else
       call vdim(nr, md)
!
       if(.not. alloced) then
          call alloc_module_arrays(nr, md)
          alloced=.true.
       endif

       call vdat(eqpsi_0, eqpsi_a, vB_T, vI_0, vp_0, &
            vr_min, vr_max, vRmag, vz_max, vZmag, vavgrmid, &
            rho_d, psi_bar, fp, qsf, beta, &
            pressure, cur, vdvdrho, varea, eqpsi, B_psi, vc)
    endif

    beta_0 = beta(1)

    R_0=vRmag 
    Z_0=vZmag 

    avgrmid=vavgrmid
    rmaj=vRmag
    B_T0=vB_T

    psi_N=5*avgrmid*B_T0/vI_0

    do i=1,nr
       eqpsi(i)=eqpsi(i)/psi_N
    enddo

    eqpsi_a=eqpsi_a/psi_N
    eqpsi_0=eqpsi_0/psi_N

    psi_a=eqpsi_a
    psi_0=eqpsi_0

    i = initialize_rcenter(1) 
    f_N=1./rcenter(psi_a)
    do i=1,nr
       fp(i)=fp(i)/f_N
    enddo
    nthg = md - 1

  end subroutine vmomin

  subroutine veq_init
    
    real :: pi, th
    real, dimension(nr, md) :: eq_psi, eqth
    integer :: i, j
    
    pi = 2*acos(0.)

    do j = 1, md
       do i = 1, nr
          eq_psi(i, j) = eqpsi(i)
          th = float(j-1)*pi/float(md-1)
          eqth(i, j) = pi - th
          R_psi(i, j) = vc(i,1) - rho_d(i)*cos(th) + vc(i,5)*cos(2*(th))
          Z_psi(i, j) = rho_d(i)*vc(i,3)*sin(th) + vc(i,3)*vc(i,5)*sin(2*(th))
       enddo
    enddo
      
    call derm(eqth,   dtm, 'T')
    call derm(R_psi,  drm, 'E')
    call derm(Z_psi,  dzm, 'O')
    call derm(B_psi,  dbm, 'E')
    call derm(eq_psi, dpm, 'E')

    call eqdcart(dpm, dpcart)
    call eqdcart(dbm, dbcart)
    call eqdcart(dtm, dtcart)

    call eqdbish(dpcart, dpbish)
    call eqdbish(dbcart, dbbish)
    call eqdbish(dtcart, dtbish)

  end subroutine veq_init

  subroutine derm(f, dfm, char)

    implicit none
    character*1 :: char
    integer i, j
    real f(:,:), dfm(:,:,:), pi
    
    pi = 2*acos(0.)

    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))
   
! assume up-down symmetry for now:
 
    select case(char)
    case('E')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))

       j=md  
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))

    case('O')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))
       
       j=md  
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))

    case('T')
       j=1
       dfm(:,j,2) = f(:,j+1)
       
       j=md  
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo
    
    do j=2,md-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo
    
  end subroutine derm

  subroutine gradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
!    include 'eik.par'
    
    integer nth_used, ntm
    character*1 char, eqf
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real tmp(2), aa(1), daa(1), rp, rpt(1)
    real dcart(nr, md, 2)      
    integer i
    
    select case (char)
    case ('B')
       dcart = dbcart	 
    case ('P')
       dcart = dpcart
    case ('R')
       dcart = dpcart
    case ('T')
       dcart = dtcart
    end select
    
    do i=-nth_used,-1
       eqf='R'
!       f(:,:) = dcart(:,:,1)
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), eqf)
       eqf='Z'
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), eqf)
       if(char == 'T') then
          grad(i,1)=-tmp(1)
          grad(i,2)=-tmp(2)
       else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       endif
    enddo

    do i=0,nth_used
       eqf='R'
       call eqitem(rgrid(i),theta(i),dcart(:,:,1),tmp(1),eqf)
       eqf='Z'
       call eqitem(rgrid(i),theta(i),dcart(:,:,2),tmp(2),eqf)
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5 * beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5 * beta_0
       enddo
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character*1 char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real tmp(2), aa(1), daa(1), rp, rpt(1)
    real, dimension(nr, md, 2) ::  dbish
    real, dimension(nr, md) :: f
    integer i

    select case(char)
    case('B') 
       dbish = dbbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       f(:,:) = dbish(:,:,1)
       call eqitem(rgrid(i), theta(i), f, tmp(1), 'R')
       f(:,:) = dbish(:,:,2)
       call eqitem(rgrid(i), theta(i), f, tmp(2), 'Z')
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
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
      
    integer :: i, j, istar, jstar
    character*1 :: char
    real :: r, thet, fstar, sign, tp, tps, theta_in
    real :: st, dt, sr, dr, pi
    real, dimension(:,:) :: f
    real, dimension(size(f,2)) :: vtheta
!    logical out
    
    pi = 2.*acos(0.)

! check for axis evaluation
      
    if(r == eqpsi(1)) then
       write(*,*) 'no evaluation at axis allowed in eqitem'
       write(*,*) r, theta_in, eqpsi(1)
       stop
    endif
    
! allow psi(r) to be a decreasing function

    sign=1.
    if(eqpsi(2) < eqpsi(1)) sign=-1.
    
    if(r < sign*eqpsi(1)) then
       write(*,*) 'r < Psi_0 in eqitem'
       write(*,*) r,sign,eqpsi(1)
       stop
    endif
      
! find r on psi mesh

! disallow evaluations on or outside the surface for now

    if(r >= eqpsi(nr)) then
       write(*,*) 'No evaluation of eqitem allowed on or outside surface'
       write(*,*) '(Could this be relaxed a bit?)'
       write(*,*) r, theta_in, eqpsi(nr), sign
       stop      
    endif
    
    istar=0
    do i=2,nr
       if(istar /= 0) cycle
       if(r < sign*eqpsi(i)) then
          dr = r - sign*eqpsi(i-1)
          sr = sign*eqpsi(i) - r
          istar=i-1
       endif
    enddo
    
! no gradients available at axis, so do not get too close

    if(istar == 1) then
       write(*,*) 'Too close to axis in eqitem'
       write(*,*) r, theta_in, eqpsi(1), eqpsi(2)
       stop
    endif
  
! Now do theta direction

    thet = mod2pi(theta_in)

! assume up-down symmetry

    tp=abs(thet)
    tps=1.
    if(char == 'Z' .and. thet /= 0.) tps=thet/abs(thet)
    
! get thet on theta mesh

    do j=1,md
       vtheta(j)=(j-1)*pi/float(md-1)
    enddo
  
! note that theta(1)=0 for vmom_eq theta   !needs to be rechecked!

    jstar=-1
    do j=1,md
       if(jstar /= -1) cycle
       if(tp < vtheta(j)) then
          dt = tp - vtheta(j-1)
          st = vtheta(j) - tp
          jstar=j-1
!        write(*,*) tp,dt,vtheta(j-1),j,vtheta(j),st
!        write(*,*) 
       endif
    enddo
      
! treat theta = pi separately
  
    if(tp == pi) then
       jstar=md-1
       dt=vtheta(jstar+1)-vtheta(jstar)
       st=0.
    endif

! use opposite area stencil to interpolate

!    if(char == 'R') i=1
!    if(char == 'Z') i=2
    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar=fstar*tps &
         /abs(eqpsi(istar+1)-eqpsi(istar)) &
         /(vtheta(jstar+1)-vtheta(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) eqpsi(istar),eqpsi(istar+1)
!     write(*,*) vtheta(jstar),vtheta(jstar+1)
      
  end subroutine eqitem

  subroutine eqdcart(dfm, dfcart)
      
    implicit none
    integer :: i,j

    real, dimension (:,:,:) :: dfm, dfcart
    real, dimension (size(dfm,1),size(dfm,2)) :: denom

    denom(:,:)=drm(:,:,1)*dzm(:,:,2)-drm(:,:,2)*dzm(:,:,1)

    dfcart(:,:,1)=   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2)= - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)

    do j=1,md
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
    
    do j=1,md
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
       enddo
    enddo

  end subroutine eqdbish

  function initialize_rcenter (init) 

    integer :: init, initialize_rcenter

    init_rcenter = .false.
    if(init == 1) init_rcenter = .true.
    initialize_rcenter = 1

  end function initialize_rcenter

  function rcenter(rp)

    use splines

    real, intent (in) :: rp
    real :: rcenter
!    integer :: init_rc = 1
    type (spline) :: s1, s5

!    if(init_rc == 1) then
    if(init_rcenter) then
       call new_spline(nr, eqpsi, vc(:,1), s1)
       call new_spline(nr, eqpsi, vc(:,5), s5)
    endif

    rcenter = splint(rp, s1) + splint(rp, s5)
    
  end function rcenter

  function invR( r, theta) 

    real :: invR, r, theta, th, f
    
    th = mod2pi(theta) 

    call eqitem(r, th, R_psi, f, 'R')
    invR = 1./f

  end function invR

  function Rpos( r, theta) 

    real :: Rpos, r, theta, th, f
    
    th = mod2pi(theta) 

    call eqitem(r, th, R_psi, f, 'R')
    Rpos = f

  end function Rpos

  function Zpos( r, theta) 

    real :: Zpos, r, theta, th, f
    
    th = mod2pi(theta) 

    call eqitem(r, th, Z_psi, f, 'Z')
    Zpos = f

  end function Zpos

  function initialize_psi (init) 

    integer :: init, initialize_psi
    
    init_psi = .false.
    if(init == 1) init_psi = .true.
    initialize_psi = 1

  end function initialize_psi

  function psi (r, theta)
   
    real, intent (in) :: r, theta
    real :: psi

    psi = r
    
  end function psi
   
  function initialize_diameter (init) 

    integer :: init, initialize_diameter

    init_diameter = .false.
    if(init == 1) init_diameter = .true.
    initialize_diameter = 1

  end function initialize_diameter

  function diameter (rp)
    
    use splines
    real :: rp, diameter
    type (spline), save :: spl

    if(init_diameter) then
       call new_spline(nr, eqpsi, rho_d, spl)
       init_diameter = .false.
    endif

    diameter = splint(rp, spl)

  end function diameter
      
  function initialize_btori (init) 

    integer :: init, initialize_btori
    
    init_btori = .false.
    if(init == 1) init_btori = .true.
    initialize_btori = 1

  end function initialize_btori

  function btori (pbar)
  
    use splines
    real :: pbar, btori
    type (spline), save :: spl

    if(init_btori) call new_spline(nr, psi_bar, fp, spl)

    btori = splint(pbar, spl)

  end function btori

  function initialize_dbtori (init) 

    integer :: init, initialize_dbtori
    
    init_dbtori = .false.
    if(init == 1) init_dbtori = .true.
    initialize_dbtori = 1

  end function initialize_dbtori

  function dbtori (pbar)
  
    use splines
    real :: pbar, dbtori
    type (spline), save :: spl

    if(init_dbtori) call new_spline(nr, psi_bar, fp, spl)

    dbtori = dsplint(pbar, spl)/(eqpsi_a-eqpsi_0)

  end function dbtori

  function initialize_q (init) 

    integer :: init, initialize_q
    
    init_q = .false.
    if(init == 1) init_q = .true.
    initialize_q = 1

  end function initialize_q

  function qfun (pbar)
  
    use splines
    real :: pbar, qfun
    type (spline), save :: spl

    if(init_q) call new_spline(nr, psi_bar, qsf, spl)

    qfun = splint(pbar, spl)

  end function qfun

  function initialize_pressure (init) 

    integer :: init, initialize_pressure
    
    init_pressure = .false.
    if(init == 1) init_pressure = .true.
    initialize_pressure = 1

  end function initialize_pressure

  function pfun (pbar)
  
    use splines
    real :: pbar, pfun
    type (spline), save :: spl

    if(init_pressure) call new_spline(nr, psi_bar, pressure, spl)

    pfun = splint(pbar, spl) * beta_0/2.

  end function pfun

  function initialize_dpressure (init) 

    integer :: init, initialize_dpressure
    
    init_dpressure = .false.
    if(init == 1) init_dpressure = .true.
    initialize_dpressure = 1

  end function initialize_dpressure

  function dpfun (pbar)
  
    use splines
    real :: pbar, dpfun
    type (spline), save :: spl

    if(init_dpressure) call new_spline(nr, psi_bar, pressure, spl)

    dpfun = dsplint(pbar, spl)/(eqpsi_a-eqpsi_0) * beta_0/2.

  end function dpfun

  function initialize_beta (init) 

    integer :: init, initialize_beta
    
    init_beta = .false.
    if(init == 1) init_beta = .true.
    initialize_beta = 1

  end function initialize_beta

  function betafun (pbar)
  
    use splines
    real :: pbar, betafun
    type (spline), save :: spl

    if(init_beta) call new_spline(nr, psi_bar, beta, spl)

    betafun = splint(pbar, spl)

  end function betafun

  function rhoofpsi( psi) 

    use splines
    real :: rhoofpsi, psi
    type (spline) :: spl
    integer :: initrho = 1

    if(initrho == 1) then
       initrho = 0
       call new_spline(nr, eqpsi, rho_d, spl)
    endif
    
    rhoofpsi = splint (psi, spl)

  end function rhoofpsi

  subroutine psiofrho(rho, psi, drhodpsi)

    use splines
    real :: rho, psi, drhodpsi, df
    type (spline) :: spl
    integer :: initpsi = 1
    
    if(initpsi == 1) then
       initpsi = 0
       call new_spline(nr, rho_d, eqpsi, spl)
    endif

! error inside of NT because not reinitialized
    
    psi = splint(rho, spl)
    df = dsplint(rho, spl)
    
    if(abs(df) <= 1.e-5) then
       write(*,*) 'drhodpsi is very large.'
       write(*,*) 'dpsidrho= ',df
       write(*,*) 'Stopping in psiofrho.'
       stop
    else
       drhodpsi=1/df
    endif
        
  end subroutine psiofrho

  function modpi (theta) 

    real, intent(in) :: theta
    real :: pi, th, modpi

    logical out 

    pi = 2*acos(0.)
    
    if(abs(theta) <= pi) then
       modpi = abs(theta)
       return
    endif

    th = abs(theta)
    out = .true.
    do while(out) 
       if(th > pi) th = th - pi
       if(th <= pi) out = .false.
    enddo
    modpi = th
    
  end function modpi

  function mod2pi (theta)
    
    real, intent(in) :: theta
    real :: pi, th, mod2pi
    logical :: out
    
    pi=2.*acos(0.)
    
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
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

  subroutine alloc_module_arrays(nr, md)
     integer :: nr, md

     allocate(rho_d(nr), qsf(nr), fp(nr), psi_bar(nr), beta(nr), pressure(nr), cur(nr))
     allocate(eqpsi(nr), vdvdrho(nr), varea(nr))
     allocate(drm(nr,md,2), dzm(nr,md,2), dpm(nr,md,2), dtm(nr,md,2), dbm(nr,md,2))
     allocate(dpcart(nr,md,2), dtcart(nr,md,2), dbcart(nr,md,2))
     allocate(R_psi(nr,md), Z_psi(nr,md), B_psi(nr,md), vc(nr, 6))

  end subroutine alloc_module_arrays

end module veq
