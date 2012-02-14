module rk_schemes

  !  Define generic Runge-Kutta scheme types of different order
  !  using the 'Butcher tableau' notation
  !  The number of stages is nearly, but not quite, the same as the scheme's order
  !
  !  John Butcher (2007)
  !  <a href='http://www.scholarpedia.org/article/Runge-Kutta_methods'>
  !  Runge-Kutta methods</a>.
  !  <a href='http://www.scholarpedia.org'>Scholarpedia</a>, 2(9):3147
  !
  !  T&M/PKNIGHT/LOGBOOK22, pp.59,60,62-64

  type :: rk1_scheme
     integer :: stages = 1
     real, dimension(1,1) :: aij = 0.0
     real, dimension(1) :: bi = 0.0
     real, dimension(1) :: ci = 0.0
  end type rk1_scheme

  type :: rk2_scheme
     integer :: stages = 2
     real, dimension(2,2) :: aij = 0.0
     real, dimension(2) :: bi = 0.0
     real, dimension(2) :: ci = 0.0
  end type rk2_scheme

  type :: rk3_scheme
     integer :: stages = 3
     real, dimension(3,3) :: aij = 0.0
     real, dimension(3) :: bi = 0.0
     real, dimension(3) :: ci = 0.0
  end type rk3_scheme

  type :: rk4_scheme
     integer :: stages = 4
     real, dimension(4,4) :: aij = 0.0
     real, dimension(4) :: bi = 0.0
     real, dimension(4) :: ci = 0.0
  end type rk4_scheme

  !  Runge-Kutta pair schemes
  !  Suitable for use in adaptive methods to control truncation errors

  type :: rk12_scheme
     !  First and second order schemes embedded
     integer :: stages = 2
     real, dimension(2,2) :: aij = 0.0
     real, dimension(2,2) :: bi = 0.0
     real, dimension(2) :: ci = 0.0
  end type rk12_scheme

  type :: rk23_scheme
     !  Second and third order schemes embedded
     integer :: stages = 4
     real, dimension(4,4) :: aij = 0.0
     real, dimension(2,4) :: bi = 0.0
     real, dimension(4) :: ci = 0.0
  end type rk23_scheme

  type :: rk45_scheme
     !  Fourth and fifth order schemes embedded
     integer :: stages = 6
     real, dimension(6,6) :: aij = 0.0
     real, dimension(2,6) :: bi = 0.0
     real, dimension(6) :: ci = 0.0
  end type rk45_scheme

  type :: rk45a_scheme
     !  Fourth and fifth order schemes embedded, but with 7 stages
     integer :: stages = 7
     real, dimension(7,7) :: aij = 0.0
     real, dimension(2,7) :: bi = 0.0
     real, dimension(7) :: ci = 0.0
  end type rk45a_scheme

  !  Define specific schemes

  real, parameter, private :: sixth = 1.0/6.0
  real, parameter, private :: third = 1.0/3.0
  real, parameter, private :: twothirds = 2.0/3.0

  type(rk1_scheme), parameter :: euler = rk1_scheme( &
       1, (/0.0/), (/1.0/), (/0.0/) )

  type(rk2_scheme), parameter :: midpoint = rk2_scheme( &
       2, &
       (/ 0.0, 0.5, 0.0, 0.0 /), &
       (/ 0.0, 1.0 /), &
       (/ 0.0, 0.5 /) )
  type(rk2_scheme), parameter :: heun = rk2_scheme( &
       2, &
       (/ 0.0, 1.0, 0.0, 0.0 /), &
       (/ 0.5, 0.5 /), &
       (/ 0.0, 1.0 /) )
  type(rk2_scheme), parameter :: atkinson = rk2_scheme( &
       2, &
       (/ 0.0, twothirds, 0.0, 0.0 /), &
       (/ 0.25, 0.75 /), &
       (/ 0.0, twothirds /) )

  type(rk4_scheme), parameter :: rk4 = rk4_scheme( & !  classical Runge-Kutta method
       4, &
       (/ 0.0, 0.5, 0.0, 0.0, &
          0.0, 0.0, 0.5, 0.0, &
          0.0, 0.0, 0.0, 1.0, &
          0.0, 0.0, 0.0, 0.0 /), &
       (/ sixth, third, third, sixth /), &
       (/ 0.0, 0.5, 0.5, 1.0/) )

  type(rk12_scheme), parameter :: simple_adaptive = rk12_scheme( &
       2, &
       (/ 0.0, 1.0, &
          0.0, 0.0 /), &
       (/ 1.0, 0.5, &
          0.0, 0.5 /), &
       (/ 0.0, 1.0 /) )

  type(rk23_scheme), parameter :: bogacki_shampine = rk23_scheme( &
       4, &
       (/ 0.0, 0.5, 0.0, 2.0/9.0, &
          0.0, 0.0, 0.75, third, &
          0.0, 0.0, 0.0, 4.0/9.0, &
          0.0, 0.0, 0.0, 0.0 /), &
       (/ 2.0/9.0, 7.0/24.0, &
          third, 0.25, &
          4.0/9.0, third, &
          0.0, 0.125 /), &
       (/ 0.0, 0.5, 0.75, 1.0 /) )

  type(rk23_scheme), parameter :: rkf3 = rk23_scheme( & !  Stoer & Bulirsch
       4, &
       (/ 0.0, 0.25, -189.0/800.0, 214.0/891.0, &
          0.0, 0.0, 729.0/800.0, 1.0/33.0, &
          0.0, 0.0, 0.0, 650.0/891.0, &
          0.0, 0.0, 0.0, 0.0 /), &
       (/ 214.0/891.0, 533.0/2106.0, &
          1.0/33.0, 0.0, &
          650.0/891.0, 800.0/1053.0, &
          0.0, -1.0/78.0 /), &
       (/ 0.0, 0.0, 0.0, 0.0 /) )

  type(rk45_scheme), parameter :: fehlberg = rk45_scheme( &
       6, &
       (/ 0.0, 0.25, 0.09375, 1932.0/2197.0, 439.0/216.0, -8.0/27.0, &
          0.0, 0.0, 0.28125, -7200.0/2197.0, -8.0, 2.0, &
          0.0, 0.0, 0.0, 7296.0/2197.0, 3680.0/513.0, -3544.0/2565.0, &
          0.0, 0.0, 0.0, 0.0, -845.0/4104.0, 1859.0/4104.0, &
          0.0, 0.0, 0.0, 0.0, 0.0, -0.275, &
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), &
       (/ 25.0/216.0, 16.0/135.0, &
          0.0, 0.0, &
          1408.0/2565.0, 6656.0/12825.0, &
          2197.0/4104.0, 28561.0/56430.0, &
          -0.2, -0.18, &
          0.0, 2.0/55.0 /), &
       (/ 0.0, 0.25, 0.375, 12.0/13.0, 1.0, 0.5 /) )

  type(rk45_scheme), parameter :: cash_karp = rk45_scheme( &
       6, &
       (/ 0.0, 0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, &
          0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875, &
          0.0, 0.0, 0.0, 1.2, -70.0/27.0, 575.0/13824.0, &
          0.0, 0.0, 0.0, 0.0, 35.0/27.0, 44275.0/110592.0, &
          0.0, 0.0, 0.0, 0.0, 0.0, 253.0/4096.0, &
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), &
       (/ 37.0/378.0, 2825.0/27648.0, &
          0.0, 0.0, &
          250.0/621.0, 18575.0/48384.0, &
          125.0/594.0, 13525.0/55296.0, &
          0.0, 277.0/14336.0, &
          512.0/1771.0, 0.25 /), &
       (/ 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /) )

  type(rk45a_scheme), parameter :: dormand_prince = rk45a_scheme( &
       7, &
       (/ 0.0, 0.2, 0.075, 44.0/45.0, 19372.0/6561.0, 9017.0/3168.0, 35.0/384.0, &
          0.0, 0.0, 0.225, -56.0/15.0, -25360.0/2187.0, -355.0/33.0, 0.0, &
          0.0, 0.0, 0.0, 32.0/9.0, 64448.0/6561.0, 46732.0/5247.0, 500.0/1113.0, &
          0.0, 0.0, 0.0, 0.0, -212.0/729.0, 49.0/176.0, 125.0/192.0, &
          0.0, 0.0, 0.0, 0.0, 0.0, -5103.0/18656.0, -2187.0/6784.0, &
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0/84.0, &
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), &
       (/ 5179.0/57600.0, 35.0/384.0, &
          0.0, 0.0, &
          7571.0/16695.0, 500.0/1113.0, &
          393.0/640.0, 125.0/192.0, &
          -92097.0/339200.0, -2187.0/6784.0, &
          187.0/2100.0, 11.0/84.0, &
          0.025, 0.0 /), &
       (/ 0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0 /) )

end module rk_schemes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fields_explicit

! NOT UP TO DATE ... DO NOT USE

  implicit none

  public :: init_fields_explicit
  public :: advance_explicit
  public :: init_phi_explicit
  public :: reset_init
!+PJK
  public :: adaptive_dt, adaptive_dt_reset, adaptive_dt_new
!-PJK
  private

  interface nodal2modal
     module procedure nodal2modal_complex3d
     module procedure nodal2modal_real3d
  end interface nodal2modal

  interface modal2nodal
     module procedure modal2nodal_complex3d
     module procedure modal2nodal_real3d
  end interface modal2nodal

  logical :: initialized = .false.
!+PJK
  logical :: adaptive_dt = .true.  !  switch for adaptive timestep algorithm
  logical :: adaptive_dt_reset = .false.
  real :: adaptive_dt_new
!-PJK

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_fields_explicit
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    
  end subroutine init_fields_explicit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_phi_explicit
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
!+PJK
    use fields_arrays, only: phiold, aparold, bparold
    use dist_fn_arrays, only: gold, g
!-PJK
    use dist_fn, only: getfieldexp

!+PJK
!    call getfieldexp (phinew, aparnew, bparnew)
!    phi = phinew; apar = aparnew; bpar = bparnew
!-PJK
    gold = g
    call getfieldexp (gold, phiold, aparold, bparold, 'g')
    phi = phiold; apar = aparold; bpar = bparold

  end subroutine init_phi_explicit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advance_explicit (istep)

    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use fields_arrays, only: phitmp1, apartmp1, bpartmp1
    use dist_fn, only: timeadv, getfieldexp, g_adjust_exp
    use dist_fn, only: gnl_1, gnl_2, gnl_3, def_parity, even
    use dist_fn_arrays, only: g, gnew, gold, vpar
!+PJK
    use collisions, only: solfp1
    use dist_fn_arrays, only: gwork, aj0
    use fields_arrays, only: phiold, aparold, bparold
    use gs2_layouts, only: g_lo, ik_idx, is_idx
    use gs2_time, only: save_dt, code_time, code_dt
    use hyper, only: hyper_diff
    use nonlinear_terms, only: add_nonlinear_terms
    use run_parameters, only: tunits, fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    use species, only: spec
    use dist_fn, only: wdrift, wcoriolis
    use mp, only: max_allreduce
!-PJK
    implicit none
    integer, intent (in) :: istep
!+PJK

    real :: time, dt, dtnew, alpha = 0.9, eps = 1.0e-5, rcushion = 0.3
    integer :: ig, ie, mz, p = 3, isgn, iglo, ik, ileft, ip, is
    integer :: modep
    complex, parameter :: z1 = (1.0, 0.0)
    complex, allocatable, dimension (:,:,:) :: gnew1, gnew2
    real, allocatable, dimension (:,:,:) :: v
    real, allocatable, dimension (:,:) :: wd

    !  Adaptive timestep variables

    real, save :: epsmach = 1.0E-15, epsr = 1.0E-6, epsa = 1.0E-6 ! both from 1.0E-6
    real, save :: ffac = 0.9 ! = (1.0 - eps1) in write-up
    real, save :: dt0,dtmin,dtmax,epsrmin,epsbig,g3dn,g23n, dtmax_damped = 10.0
    real, save :: reps,epsar,zfacu,zfacut,zfacd,zfacdt,zrsord
    real, save :: rs,zhf,zdt
    real :: gomax,gnmax
    integer, save :: nfacup = 5, nfacdn = 10, first_call = 1
    logical, save :: lfail = .false., ldummy

    !  Initialise adaptive timestep stuff

    if (adaptive_dt) then
       if (first_call == 1) then
          epsrmin = 1.0E-12 + 2.0*epsmach
          epsbig = 33.0*epsmach
          epsr = max(epsrmin,epsr)
          reps = 2.0/epsr
          epsar = 2.0*epsa/epsr
          zfacu = real(nfacup)
          zfacut = (ffac/zfacu)**p
          zfacd = 1.0/nfacdn
          zfacdt = (ffac/zfacd)**p
          zrsord = 1.0/p

          dt0 = code_dt
          dtmin = epsbig*dt0
       end if
    end if
!-PJK

    gold = g! ; phitmp1 = phi ; apartmp1 = apar ; bpartmp1 = bpar

!+PJK Following two lines are present in advance_implicit...
!    call antenna_amplitudes (apar_ext)
!    if (allocated(kx_shift)) call exb_shear (gnew, phinew, aparnew, bparnew) 
!-PJK

    !  Original code... 'Simple Explicit'
    !g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew
    !call getfieldexp (gnew, phinew, aparnew, bparnew, 'g')
    !call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, -1)
    !g = gnew
    !gold = gnew ; phiold = phinew ; aparold = aparnew ; bparold = bparnew
    !goto 100

    !  Discontinuous Galerkin + Runge-Kutta, adaptive timestep method

    !  mz is the number of finite elements,
    !  each representing p finite difference grid points

    mz = (2*ntgrid+1)/p

    allocate( &
         gnew1(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         gnew2(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         v(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         wd(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc) )

    modep = 0
    !if (present(mode)) modep = mode

    !  Apply normalisations to vpar, wdrift

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid,ntgrid
          v(ig,:,iglo) = spec(is)%tz * vpar(ig,:,iglo)
          wd(ig,iglo)  = spec(is)%tz * (wdrift(ig,iglo) + wcoriolis(ig,iglo))
       end do
    end do

    !  Evaluate fields at old timestep

    call getfieldexp (g, phi, apar, bpar, 'g')
    call add_nonlinear_terms (gnl_1, gnl_2, gnl_3, &
         phi, apar, bpar, istep, 0.0, z1)

    !  Original call being replaced...
    !call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)

    time = code_time  !  scaling by tunits(ik) occurs inside rk_advance2
    dt = code_dt      !  ditto

    call g_adjust_exp (g, apar, fapar)  !  convert g to i

    !  On first call only, approximate the 'best' choice for the
    !  initial time step

    if (adaptive_dt) then
       adaptive_dt_reset = .false.
       if (first_call == 1) then
          call adaptive_dt0(g,phi,apar,bpar,v,wd,istep,p,mz,time, &
               dt0,dt,dtmin,dtmax,epsr,epsa)
          if (dt /= code_dt) then
             adaptive_dt_reset = .true.
             adaptive_dt_new = dt
          end if
          first_call = 0
       end if
    end if

    !  Runge-Kutta loop

!    call pjk_advance(g,gnew1,gnew2,phi,apar,bpar,v,wd,istep,p,mz,time,dt)
    call rk_advance2(g,gnew1,gnew2,phi,apar,bpar,v,wd,istep,p,mz,time,dt)

    !  Adaptive timestep algorithm

    if (adaptive_dt) then

       !  Actual difference between 2nd and 3rd order estimates
       g3dn = maxval(abs(gnew2-gnew1)) ; call max_allreduce(g3dn)

       gomax = maxval(abs(g)) ; call max_allreduce(gomax)
       gnmax = maxval(abs(gnew2)) ; call max_allreduce(gnmax)

       !  (reps *) Maximum allowed difference between 2nd and 3rd order estimates
       g23n = gomax + gnmax + epsar
       rs = reps*g3dn/g23n

       !write(*,*) '    |g2|, |g3| = ',maxval(abs(gnew1)),maxval(abs(gnew2))
       !write(*,*) 'diff, max diff = ',g3dn,g23n/reps,rs

       if (rs <= 1.0) then
          !  Step succeeded, but if possible increase timestep for next step

          zhf = min(ffac/rs**zrsord, zfacu)
          if (lfail) then
             zhf = 1.0
             dtmax_damped = min(dtmax_damped,dt)
          end if
          zdt = max(zhf*dt,dtmin)
          !zdt = min(zdt,dtmax_damped)  !  Uncomment for 'damped' dt
          dt = zdt
          gnew = gnew1
          lfail = .false.
          !write(*,*) 'Success: old dt, new dt = ',code_dt, dt

       else
          !  Step failed, reduce timestep

          zhf = max(ffac/rs**zrsord, zfacd)
          zdt = zhf*dt
          dt = zdt

          !  Try again unless timestep too small
          if (dt <= dtmin) then
             write(*,*) 'Adaptive timestep too small... stopping!'
             stop
          end if
          gnew = gold
          lfail = .true.
          !write(*,*) 'Failure: old dt, new dt = ',code_dt, dt

       end if

       !  If necessary, force time-step change in main program

       if (dt /= code_dt) then
          adaptive_dt_reset = .true.
          adaptive_dt_new = dt
       end if

    else
       gnew = gnew1
    end if

    !  Evaluate fields at new timestep

    call getfieldexp (gnew, phinew, aparnew, bparnew, 'i')
    call g_adjust_exp (gnew, aparnew, -fapar)  !  convert from i back to g
    g = gnew

    !  From timeadv (implicit solver) - replaced g0 with gwork
    call hyper_diff (gnew, gwork, phinew, bparnew)
    call kill (gnew, gwork, phinew, bparnew)
    call solfp1 (gnew, g, gwork, phi, apar, bpar, phinew, aparnew, bparnew, modep)

    if (def_parity) then
       if (even) then
          gnew(-ntgrid:-1, 1,:) = gnew( ntgrid: 1:-1,2,:)
          gnew( 1: ntgrid, 1,:) = gnew(-1:-ntgrid:-1,2,:)
       else
          gnew( 1: ntgrid, 1,:) = -gnew(-1:-ntgrid:-1,2,:)
          gnew(-ntgrid:-1, 1,:) = -gnew( ntgrid: 1:-1,2,:)
       end if
    end if

100 continue

    if (allocated(gnew1)) then
       deallocate(gnew1,gnew2,v,wd)
    end if

  end subroutine advance_explicit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine reset_init
    initialized = .false.
  end subroutine reset_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nodal2modal_complex3d(y,mz,p,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: mz, p  !  mz finite elements, each with p points
    integer, intent(in) :: lb1, ub1, lb2, ub2, lb3, ub3
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    integer :: i, j, k, kk, element
    complex, dimension(mz*p) :: ytemp
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3) :: ycopy

    ycopy = y

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = ycopy(:,i,j)

          do element = 1,mz
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 2) then
                y(kk,i,j) = 0.5*(ytemp(k)+ytemp(k+1))
                y(kk+1,i,j) = ytemp(k+1) - ytemp(k)
             else if (p == 3) then
                y(kk,i,j) = 0.125*(3.0*ytemp(k) + 2.0*ytemp(k+1) + 3.0*ytemp(k+2))
                y(kk+1,i,j) = 0.75*(ytemp(k+2) - ytemp(k))
                y(kk+2,i,j) = 0.75*(ytemp(k) - 2.0*ytemp(k+1) + ytemp(k+2))
             else
                write(*,*) 'Bad p value...'
                stop
             end if
          end do

       end do
    end do

  end subroutine nodal2modal_complex3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nodal2modal_real3d(y,mz,p,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: mz, p  !  mz finite elements, each with p points
    integer, intent(in) :: lb1, ub1, lb2, ub2, lb3, ub3
    real, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    integer :: i, j, k, kk, element
    real, dimension(mz*p) :: ytemp
    real, dimension(lb1:ub1,lb2:ub2,lb3:ub3) :: ycopy

    ycopy = y

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = ycopy(:,i,j)

          do element = 1,mz
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 2) then
                y(kk,i,j) = 0.5*(ytemp(k)+ytemp(k+1))
                y(kk+1,i,j) = ytemp(k+1) - ytemp(k)
             else if (p == 3) then
                y(kk,i,j) = 0.125*(3.0*ytemp(k) + 2.0*ytemp(k+1) + 3.0*ytemp(k+2))
                y(kk+1,i,j) = 0.75*(ytemp(k+2) - ytemp(k))
                y(kk+2,i,j) = 0.75*(ytemp(k) - 2.0*ytemp(k+1) + ytemp(k+2))
             else
                write(*,*) 'Bad p value...'
                stop
             end if
          end do

       end do
    end do

  end subroutine nodal2modal_real3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine modal2nodal_complex3d(y,mz,p,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: mz, p  !  mz finite elements, each with p points
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    real, parameter :: twothirds = 2.0/3.0
    real, parameter :: sixth = 1.0/6.0
    integer :: i, j, k, kk, element
    complex, dimension(mz*p) :: ytemp

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = y(:,i,j)

          do element = 1,mz
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 2) then
                y(kk,i,j) = ytemp(k) - 0.5*ytemp(k+1)
                y(kk+1,i,j) = ytemp(k) + 0.5*ytemp(k+1)
             else if (p == 3) then
                y(kk,i,j) = ytemp(k) - twothirds*ytemp(k+1) + sixth*ytemp(k+2)
                y(kk+1,i,j) = ytemp(k) - 0.5*ytemp(k+2)
                y(kk+2,i,j) = ytemp(k) + twothirds*ytemp(k+1) + sixth*ytemp(k+2)
             else
                write(*,*) 'Bad p value...'
                stop
             end if
          end do

       end do
    end do

  end subroutine modal2nodal_complex3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine modal2nodal_real3d(y,mz,p,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: mz, p  !  mz finite elements, each with p points
    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    real, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    real, parameter :: twothirds = 2.0/3.0
    real, parameter :: sixth = 1.0/6.0
    integer :: i, j, k, kk, element
    real, dimension(mz*p) :: ytemp

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = y(:,i,j)

          do element = 1,mz
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 2) then
                y(kk,i,j) = ytemp(k) - 0.5*ytemp(k+1)
                y(kk+1,i,j) = ytemp(k) + 0.5*ytemp(k+1)
             else if (p == 3) then
                y(kk,i,j) = ytemp(k) - twothirds*ytemp(k+1) + sixth*ytemp(k+2)
                y(kk+1,i,j) = ytemp(k) - 0.5*ytemp(k+2)
                y(kk+2,i,j) = ytemp(k) + twothirds*ytemp(k+1) + sixth*ytemp(k+2)
             else
                write(*,*) 'Bad p value...'
                stop
             end if
          end do

       end do
    end do

  end subroutine modal2nodal_real3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine adaptive_dt0(y,phi,apar,bpar,v,wd,istep,p,mz,t,dt0,dt,dtmin,dtmax,epsr,epsa)

    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use kt_grids, only: naky, ntheta0
    use mp, only: max_allreduce

    implicit none

    !  Arguments

    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: y  !  g_old (actually i_old)
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), &
         intent(in) :: phi, apar, bpar  !  fields
    real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: v
    real, dimension(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: wd
    integer, intent(in) :: istep
    integer, intent(in) :: p   !  order of the spatial Legendre fit
    integer, intent(in) :: mz  !  number of finite elements
    real, intent(in) :: t  !  code_time
    real, intent(in) :: dt0  !  code_dt
    real, intent(out) :: dt  !  'best' choice for initial timestep
    real, intent(in) :: dtmin !  minimum allowed timestep
    real, intent(out) :: dtmax !  estimate of maximum allowed timestep
    real, intent(in) :: epsr  !  relative tolerance level
    real, intent(in) :: epsa  !  absolute tolerance level

    !  Local variables

    integer :: ix, i, j, ie, ik, iglo
    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    real :: gdotn,gn,epst
    complex, allocatable, dimension(:,:,:) :: fluxfn
    complex, allocatable, dimension(:,:,:) :: src
    complex, allocatable, dimension(:,:,:) :: dy

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lb1 = -ntgrid ; ub1 = ntgrid
    lb2 = 1 ; ub2 = 2
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_alloc

    allocate(fluxfn(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(src(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy(lb1:ub1,lb2:ub2,lb3:ub3))

    !  Evaluate the source terms

    call fluxfn_dg(dt0,phi,apar,bpar,istep,fluxfn,src)  !  All nodal

    !  Calculate dg/dt using the Discontinuous Galerkin scheme

    call dydt_dg(t,dt0,y,dy,p,mz,v,wd,fluxfn,src,lb1,ub1,lb3,ub3)

    call modal2nodal(dy,mz,p,lb1,ub1,1,2,lb3,ub3)

    gdotn = maxval(abs(dy)) ; call max_allreduce(gdotn)
    gn = maxval(abs(y)) ; call max_allreduce(gn)

    epst = epsr*gn + epsa

    dtmax = (epst/gdotn)**(1.0/p)

    dt = min(dt0,dtmax)
    dt = max(dt,dtmin)

  end subroutine adaptive_dt0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rk_advance2(y,ynew1,ynew2,phi,apar,bpar,v,wd,istep,p,mz,t,h)

    use rk_schemes, rk => rkf3
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use dist_fn, only: getfieldexp
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx
    use run_parameters, only: tunits
    use kt_grids, only: naky, ntheta0

    implicit none

    !  Arguments

    !  Two output values for ynew: one for each of the two orders of the
    !  chosen adaptive RK scheme

    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: y  !  g_old (actually i_old)
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(out) :: ynew1
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(out) :: ynew2
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), &
         intent(in) :: phi, apar, bpar  !  fields
    real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: v
    real, dimension(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: wd
    integer, intent(in) :: istep
    integer, intent(in) :: p   !  order of the spatial Legendre fit
    integer, intent(in) :: mz  !  number of finite elements
    real, intent(in) :: t  !  code_time
    real, intent(in) :: h  !  code_dt

    !  Local variables

    integer :: ix, i, j, ie, ik, iglo
    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    complex, allocatable, dimension(:,:,:) :: fluxfn
    complex, allocatable, dimension(:,:,:) :: src
    complex, allocatable, dimension(:,:,:) :: dy, dy_modal
    complex, allocatable, dimension(:,:,:) :: aksum
    complex, allocatable, dimension(:,:,:,:) :: k
    complex, allocatable, dimension(:,:,:) :: ytmp, ymodal
    complex, allocatable, dimension(:,:,:) :: ynew1_modal, ynew2_modal

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lb1 = -ntgrid ; ub1 = ntgrid
    lb2 = 1 ; ub2 = 2
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_alloc

    allocate(ytmp(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ymodal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew1_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew2_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(fluxfn(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(src(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(aksum(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(k(rk%stages,lb1:ub1,lb2:ub2,lb3:ub3))

    !  Evaluate the source terms

    call fluxfn_dg(h,phi,apar,bpar,istep,fluxfn,src)  !  All nodal

    !  Calculate dg/dt using the Discontinuous Galerkin scheme

    call dydt_dg(t,h,y,dy_modal,p,mz,v,wd,fluxfn,src,lb1,ub1,lb3,ub3)
    k(1,:,:,:) = dy_modal

!    !  Euler advance (modal space)
!    ymodal = y
!    call nodal2modal(ymodal,mz,p,lb1,ub1,1,2,lb3,ub3)
!    ynew1_modal = ymodal + h * dy_modal
!    ynew1 = ynew1_modal
!    call modal2nodal(ynew1,mz,p,lb1,ub1,1,2,lb3,ub3)
!    ynew2 = ynew1
!    return

!    !  Euler advance (nodal space)
!    dy = dy_modal
!    call modal2nodal(dy,mz,p,lb1,ub1,1,2,lb3,ub3)
!    ynew1 = y + h * dy
!    ynew2 = ynew1
!    return

    !  Runge-Kutta loop

    do i = 2,rk%stages
       aksum = 0.0
       do j = 1,i-1
          aksum(:,:,:) = aksum(:,:,:) + rk%aij(i,j)*k(j,:,:,:)
       end do

       dy = h*aksum
       call modal2nodal(dy,mz,p,lb1,ub1,1,2,lb3,ub3)

       ytmp = y + dy

!       !  Multiply aksum by tunits(ik) so that y+h*aksum are okay
!       !  although it appears that tunits(ik) == 1.0 at present
!       do iglo = g_lo%llim_proc, g_lo%ulim_proc
!          ik = ik_idx(g_lo,iglo)
!          aksum(:,:,iglo) = aksum(:,:,iglo)*tunits(ik)
!       end do

       !  Calculate here the field equations for g=y+h*aksum
       !  and pass these into fluxfn_dg for use in evaluating the source terms

       call getfieldexp(ytmp,phitmp,apartmp,bpartmp,'i')
       call fluxfn_dg(h,phitmp,apartmp,bpartmp,istep,fluxfn,src)
       !  rk%ci(*) = 0.0 for rkf3, so no need to mult by tunits - but beware!

       call dydt_dg(t+h*rk%ci(i),h,ytmp,dy_modal,p,mz,v,wd,fluxfn,src, &
            lb1,ub1,lb3,ub3)
       k(i,:,:,:) = dy_modal
    end do

    ymodal = y
    call nodal2modal(ymodal,mz,p,lb1,ub1,1,2,lb3,ub3)
    ynew1_modal = ymodal
    ynew2_modal = ymodal

    do iglo = g_lo%llim_proc,g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       do i = 1,rk%stages
          ynew1_modal(:,:,iglo) = ynew1_modal(:,:,iglo) + &
               h*tunits(ik)*rk%bi(1,i)*k(i,:,:,iglo)
          ynew2_modal(:,:,iglo) = ynew2_modal(:,:,iglo) + &
               h*tunits(ik)*rk%bi(2,i)*k(i,:,:,iglo)
!  Non-adaptive RK types only...
!          ynew1_modal(:,:,iglo) = ynew1_modal(:,:,iglo) + &
!               h*tunits(ik)*rk%bi(i)*k(i,:,:,iglo)
!          ynew2_modal = ynew1_modal
       end do
    end do

    ynew1 = ynew1_modal
    call modal2nodal(ynew1,mz,p,lb1,ub1,1,2,lb3,ub3)
    ynew2 = ynew2_modal
    call modal2nodal(ynew2,mz,p,lb1,ub1,1,2,lb3,ub3)

    !  Boundary conditions...

!    ynew1(-ntgrid,1,:) = y(-ntgrid,1,:)
!    ynew2(-ntgrid,1,:) = y(-ntgrid,1,:)
!    ynew1( ntgrid,2,:) = y( ntgrid,2,:)
!    ynew2( ntgrid,2,:) = y( ntgrid,2,:)

  end subroutine rk_advance2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dydt_dg(t,dt,g,dgdt,p,mz,v,wd,f,source,lb1,ub1,lb3,ub3)

!    use dist_fn, only: wdrift, wcoriolis
    use gs2_layouts, only: g_lo, ik_idx!, is_idx
!    use run_parameters, only: tunits
!    use species, only: spec

    implicit none

    integer, intent(in) :: p, mz
    integer, intent(in) :: lb1,ub1,lb3,ub3
    real, intent(in) :: t, dt
    real, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: v
    real, dimension(lb1:ub1,lb3:ub3), intent(in) :: wd
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: g
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: f
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: source
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(out) :: dgdt  !  modal

    integer :: ip, ip2
!    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    integer :: i, j, kk, element, ig, ik, is, isgn, iglo

    complex, parameter :: zi = (0.0,1.0)
    real :: ddt

    complex, dimension(lb1:ub1,2,lb3:ub3) :: iwg, iwg_modal
    complex, dimension(lb1:ub1,2,lb3:ub3) :: gplusf, gplusf_modal
    complex, dimension(lb1:ub1,2,lb3:ub3) :: mgplusf, mgplusf_modal
    complex, dimension(lb1:ub1,2,lb3:ub3) :: vterm, vterm_modal
    complex, dimension(lb1:ub1,2,lb3:ub3) :: smodal

    real, dimension(2,4) :: mat_p2p, mat_p2m
    real, dimension(3,6) :: mat_p3p, mat_p3m

    complex, dimension(2*p) :: ftemp

    !  Matrix, p=2, positive v
    mat_p2p = real(reshape( &
         source = (/ 1, -3, 1, -3, -1, 3, -1, -3 /), &
         shape = (/ 2,4 /) ))
    !  Matrix, p=2, negative v
    mat_p2m = real(reshape( &
         source = (/ 1, 3, -1, 3, -1, -3, 1, 3 /), &
         shape = (/ 2,4 /) ))
    !  Matrix, p=3, positive v
    mat_p3p = real(reshape( &
         source = (/ 1, -3, 5, 1, -3, 5, 1, -3, 5, -1, 3, -5, -1, -3, 5, -1, -3, -5 /), &
         shape = (/ 3,6 /) ))
    !  Matrix, p=3, negative v
    mat_p3m = real(reshape( &
         source = (/ 1, 3, 5, -1, 3, 5, 1, -3, 5, -1, -3, -5, 1, 3, 5, -1, -3, -5 /), &
         shape = (/ 3,6 /) ))

    !  Evaluate -i*wd*g term

    do iglo = lb3,ub3
       do isgn = 1,2
          do ig = lb1,ub1
             iwg(ig,isgn,iglo) = -zi * wd(ig,iglo) * g(ig,isgn,iglo) / dt
          end do
       end do
    end do
    iwg_modal = iwg
    call nodal2modal(iwg_modal,mz,p,lb1,ub1,1,2,lb3,ub3)

    !  Add g to F

    gplusf = g + f
    gplusf_modal = gplusf
    call nodal2modal(gplusf_modal,mz,p,lb1,ub1,1,2,lb3,ub3)

    !  Perform matrix multiplication, using the correct matrix above
    !  For speed, ought to store and use transpose...

    mgplusf_modal = 0.0

    if (p == 2) then

       ! do iglo = lb3,ub3
       !    do i = 1,2

       !       if (i == 1) then  !  positive parallel velocity

       !          do element = 1,mz
       !             !  Fill ftemp with correct elements from the flux function
       !             if (element /= 1) then
       !                kk = lb1 + p*(element-2)  !  take values from element to the left
       !                ftemp(1:p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             else
       !                kk = lb1 + p*(mz-1)  !  assumes periodic BCs for now
       !                ftemp(1:p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             end if
       !             kk = lb1 + p*(element-1)  !  take values from this element
       !             ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1,i,iglo)

       !             do ip = 1,p
       !                kk = lb1 + p*(element-1) + ip-1
       !                do ip2 = 1,2*p
       !                   mgplusf_modal(kk,i,iglo) = mgplusf_modal(kk,i,iglo) + mat_p2p(ip,ip2)*ftemp(ip2)
       !                end do
       !             end do

       !          end do

       !       else  !  v < 0.0

       !          do element = 1,mz
       !             !  Fill ftemp with correct elements from the flux function
       !             kk = lb1 + p*(element-1)  !  take values from this element
       !             ftemp(1:p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             if (element /= mz) then
       !                kk = lb1 + p*element  !  take values from element to the right
       !                ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             else
       !                kk = lb1  !  assumes periodic BCs for now
       !                ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             end if

       !             do ip = 1,p
       !                kk = lb1 + p*(element-1) + ip-1
       !                do ip2 = 1,2*p
       !                   mgplusf_modal(kk,i,iglo) = mgplusf_modal(kk,i,iglo) + mat_p2m(ip,ip2)*ftemp(ip2)
       !                end do
       !             end do

       !          end do

       !       end if

       !    end do
       ! end do

    else if (p == 3) then

       do iglo = lb3,ub3
          do isgn = 1,2

             if (isgn == 1) then  !  positive parallel velocity
                                  !  v(*,isgn==1,*) >= 0.0, <= 0.0 for isgn==2
                                  !  Checked to be true.

                do element = 1,mz
                   !  Fill ftemp with correct elements from the flux function
                   if (element /= 1) then
                      kk = lb1 + p*(element-2)  !  take values from element to left
                      ftemp(1:p) = gplusf_modal(kk:kk+p-1,isgn,iglo)
                   else
                      ftemp(1:p) = 0.0  !  'ghost' element to left contains zeroes
                   end if
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1,isgn,iglo)

                   do ip = 1,p
                      kk = lb1 + p*(element-1) + ip-1
                      do ip2 = 1,2*p
                         mgplusf_modal(kk,isgn,iglo) = mgplusf_modal(kk,isgn,iglo) &
                              + mat_p3p(ip,ip2)*ftemp(ip2)
                      end do
                   end do

                end do

             else  !  v < 0.0

                do element = 1,mz
                   !  Fill ftemp with correct elements from the flux function
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(1:p) = gplusf_modal(kk:kk+p-1,isgn,iglo)
                   if (element /= mz) then
                      kk = lb1 + p*element  !  take values from element to right
                      ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1,isgn,iglo)
                   else
                      ftemp(p+1:2*p) = 0.0  !  'ghost' element to right contains zeros
                   end if

                   do ip = 1,p
                      kk = lb1 + p*(element-1) + ip-1
                      do ip2 = 1,2*p
                         mgplusf_modal(kk,isgn,iglo) = mgplusf_modal(kk,isgn,iglo) &
                              + mat_p3m(ip,ip2)*ftemp(ip2)
                      end do
                   end do

                end do

             end if

          end do
       end do

    end if

    !  Multiply this term by v/dz; must be done in nodal form
    !  Denominator dz is element size, not finite difference mesh separation,
    !  so need to divide by p.
    !  v is non-dimensional, and contains a factor dt/(F.D. dz), so we also have
    !  to divide by dt

    mgplusf = mgplusf_modal
    call modal2nodal(mgplusf,mz,p,lb1,ub1,1,2,lb3,ub3)

    do iglo = lb3,ub3 
       do isgn = 1,2
          do ig = lb1,ub1
             vterm(ig,isgn,iglo) = 1.0/(p*dt) * v(ig,isgn,iglo) &
                  * mgplusf(ig,isgn,iglo)
          end do
       end do
    end do
    vterm_modal = vterm
    call nodal2modal(vterm_modal,mz,p,lb1,ub1,1,2,lb3,ub3)

    !  Source term

    smodal = source
    call nodal2modal(smodal,mz,p,lb1,ub1,1,2,lb3,ub3)

    !  Finally sum all the modal terms to get dg/dt (modal)

    do iglo = lb3,ub3
       do isgn = 1,2
          do i = lb1,ub1 ! modes
             dgdt(i,isgn,iglo) = &
                  iwg_modal(i,isgn,iglo) &
                  + vterm_modal(i,isgn,iglo) &
                  + smodal(i,isgn,iglo)
!             dgdt(i,isgn,iglo) = &
!                  iwg(i,isgn,iglo) &
!                  + vterm(i,isgn,iglo) &
!                  + source(i,isgn,iglo)
          end do
       end do
    end do

!    isgn=1 ; iglo = lb3+4
!    do i = lb1,ub1
!       write(*,10) i, &
!            real(iwg(i,isgn,iglo)), &
!            real(vterm(i,isgn,iglo)), &
!            real(source(i,isgn,iglo))
!    end do
!10 format(i4, 1pe13.4e3,1pe13.4e3,1pe13.4e3)
!stop

  end subroutine dydt_dg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fluxfn_dg(dt,phi,apar,bpar,istep,f,src)

    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx
    use dist_fn, only: get_source_term_exp, wdrift, wstar
    use run_parameters, only: tunits
    use kt_grids, only: naky, ntheta0

    implicit none

    !  Arguments

    real, intent(in) :: dt  !  code_dt
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), intent(in) :: phi,apar,bpar
    integer, intent (in) :: istep
    complex, dimension(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(out) :: f  !  flux function
    complex, dimension(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(out) :: src  !  source term

    !  Local variables

    integer :: iglo, isgn, ik

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)

       do isgn = 1,2
          call get_source_term_exp(phi,apar,bpar,istep,isgn,iglo, &
               f(:,isgn,iglo), src(:,isgn,iglo))
          !  source is actually (2.dt.S), therefore need to divide by dt
          src(:,isgn,iglo) = src(:,isgn,iglo)/(2.0*dt*tunits(ik))
       end do
    end do

  end subroutine fluxfn_dg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pjk_advance(y,ynew1,ynew2,phi,apar,bpar,v,wd,istep,p,mz,t,h)

    !  Simple attempted replication of the original method for advancing the
    !  distribution function

    use rk_schemes, rk => rkf3
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use dist_fn, only: getfieldexp, get_source_term_exp, get_source_term
    use dist_fn, only: bkdiff, t0, source0, gamma0, omega0
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, is_idx
    use run_parameters, only: tunits
    use kt_grids, only: naky, ntheta0
    use constants
    use species, only: spec

    implicit none

    !  Arguments

    !  Two output values for ynew: currently set to be the same...

    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: y  !  g_old (actually i_old)
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(out) :: ynew1
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(out) :: ynew2
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), &
         intent(in) :: phi, apar, bpar  !  fields
    real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: v
    real, dimension(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc), &
         intent(in) :: wd
    integer, intent(in) :: istep
    integer, intent(in) :: p   !  order of the spatial Legendre fit
    integer, intent(in) :: mz  !  number of finite elements
    real, intent(in) :: t  !  code_time
    real, intent(in) :: h  !  code_dt

    !  Local variables

    integer :: ix, i, j, ie, ik, iglo, isgn, ig, is
    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    complex, allocatable, dimension(:,:,:) :: fluxfn
    complex, allocatable, dimension(:,:,:) :: src
    real :: sss
    complex :: sourcefac

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lb1 = -ntgrid ; ub1 = ntgrid
    lb2 = 1 ; ub2 = 2
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc

    allocate(fluxfn(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(src(lb1:ub1,lb2:ub2,lb3:ub3))

    fluxfn = (0.0, 0.0)

    !  wd = omega.dt
    !  v = v*dt/dz

    !  Evaluate the source terms
    !  src = 2.S0.dt as far as I can make out...

    !  value of sourcefac makes no apparent difference to results
    !  sourcefac = (1.0, 0.0)
    if (t > t0) then
       sourcefac = source0*exp(-zi*omega0*t+gamma0*t)
    else
       sourcefac = (0.5 - 0.5*cos(pi*t/t0))*exp(-zi*omega0*t+gamma0*t)
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1,2
          call get_source_term_exp(phi,apar,bpar,istep,isgn,iglo, &
               fluxfn(:,isgn,iglo), src(:,isgn,iglo))

          !  N.B. Need to comment out matrix multiplications section in dist_fn.f90
          !  if the original get_source_term call is used here

!          call get_source_term(phi,apar,bpar,phi,apar,bpar,istep,isgn,iglo, &
!               sourcefac, src(:,isgn,iglo))
       end do
    end do

    !  Do the time-advance; left to right for v > 0, right to left for v < 0

    sss = bkdiff(1)  !  = bakdif

    do iglo = g_lo%llim_proc, g_lo%ulim_proc

       !  v > 0
       ynew1(-ntgrid,1,iglo) = y(-ntgrid,1,iglo)
       do ig = -ntgrid,ntgrid-1
          ynew1(ig+1,1,iglo) = 1.0/(1.0+sss) * ( &
               src(ig,1,iglo) - (1.0-sss)*ynew1(ig,1,iglo) + &
               (1.0 - zi*wd(ig,iglo))*((1.0+sss)*y(ig+1,1,iglo) + (1.0-sss)*y(ig,1,iglo)) &
               - 2.0*v(ig,1,iglo)*(y(ig+1,1,iglo)+fluxfn(ig+1,1,iglo)-y(ig,1,iglo)-fluxfn(ig,1,iglo)) &
               )
       end do

       !  v < 0
       ynew1(ntgrid,2,iglo) = y(ntgrid,2,iglo)
       do ig = ntgrid-1,-ntgrid,-1
          ynew1(ig,2,iglo) = 1.0/(1.0+sss) * ( &
               src(ig,2,iglo) - (1.0-sss)*ynew1(ig+1,2,iglo) + &
               (1.0 - zi*wd(ig,iglo)) * ((1.0-sss)*y(ig+1,2,iglo) + (1.0+sss)*y(ig,2,iglo)) &
               - 2.0*v(ig,2,iglo)*(y(ig+1,2,iglo)+fluxfn(ig+1,2,iglo)-y(ig,2,iglo)-fluxfn(ig,2,iglo)) &
               )
       end do

    end do

    ynew2 = ynew1

  end subroutine pjk_advance

end module fields_explicit
