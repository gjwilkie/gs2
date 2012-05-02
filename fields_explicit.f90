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

module dg_scheme

  implicit none

  !  Variables and routines required for use with the Discontinuous Galerkin scheme

  public :: adaptive_dt, adaptive_dt_reset, adaptive_dt_new

  interface nodal2modal
     module procedure nodal2modal_complex3d
     module procedure nodal2modal_real3d
  end interface nodal2modal

  interface modal2nodal
     module procedure modal2nodal_complex3d
     module procedure modal2nodal_real3d
  end interface modal2nodal

  !  Primary DG variables

  integer :: p = 3  !  scheme order, = number of FD grid points per finite element
  integer :: ne     !  number of finite elements

  !  Normalised parallel velocity, omega drift

  real, allocatable, dimension (:,:,:) :: v
  real, allocatable, dimension (:,:) :: wd

  !  Other stuff needed to be passed to the low-level routines

  integer :: istep_dg

  !  Adaptive timestep variables

  logical :: adaptive_dt = .false.  !  master switch for adaptive timestep algorithm
  logical :: adaptive_dt_reset = .false.
  real :: adaptive_dt_new

  real, save :: epsmach = 1.0E-15, epsr = 1.0E-6, epsa = 1.0E-6 ! both from 1.0E-6
  real, save :: ffac = 0.9 ! = (1.0 - eps1) in write-up
  real, save :: dt0,dtmin,dtmax,epsrmin,epsbig,g3dn,g23n, dtmax_damped = 10.0
  real, save :: reps,epsar,zfacu,zfacut,zfacd,zfacdt,zrsord
  real, save :: rs,zhf,zdt
  real :: gomax,gnmax
  integer, save :: nfacup = 5, nfacdn = 10, first_call = 1
  logical, save :: lfail = .false., ldummy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nodal2modal_complex3d(y,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: lb1, ub1, lb2, ub2, lb3, ub3
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    integer :: i, j, k, kk, element
    complex, dimension(ne*p) :: ytemp
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3) :: ycopy

    !  There are ne finite elements, each with p points

    ycopy = y

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = ycopy(:,i,j)

          do element = 1,ne
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 3) then
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

  subroutine nodal2modal_real3d(y,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: lb1, ub1, lb2, ub2, lb3, ub3
    real, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    integer :: i, j, k, kk, element
    real, dimension(ne*p) :: ytemp
    real, dimension(lb1:ub1,lb2:ub2,lb3:ub3) :: ycopy

    !  There are ne finite elements, each with p points

    ycopy = y

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = ycopy(:,i,j)

          do element = 1,ne
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 3) then
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

  subroutine modal2nodal_complex3d(y,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    real, parameter :: twothirds = 2.0/3.0
    real, parameter :: sixth = 1.0/6.0
    integer :: i, j, k, kk, element
    complex, dimension(ne*p) :: ytemp

    !  There are ne finite elements, each with p points

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = y(:,i,j)

          do element = 1,ne
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 3) then
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

  subroutine modal2nodal_real3d(y,lb1,ub1,lb2,ub2,lb3,ub3)

    implicit none

    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    real, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y

    real, parameter :: twothirds = 2.0/3.0
    real, parameter :: sixth = 1.0/6.0
    integer :: i, j, k, kk, element
    real, dimension(ne*p) :: ytemp

    !  There are ne finite elements, each with p points

    do j = lb3,ub3
       do i = lb2,ub2
          ytemp = y(:,i,j)

          do element = 1,ne
             !  k is the index of first point of this finite element in ytemp;
             !  ytemp starts at index number 1.  kk is the equivalent index in y
             !  since this does not necessarily start at 1
             k = p*(element-1) + 1
             kk = lb1 + k-1

             if (p == 3) then
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

end module dg_scheme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fields_explicit

  use dg_scheme

  implicit none

  public :: init_fields_explicit
  public :: advance_explicit
  public :: init_phi_explicit
  public :: reset_init

  private

  logical :: initialized = .false.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_fields_explicit

    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids

    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    call init_dg_scheme
    
  end subroutine init_fields_explicit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_phi_explicit

    use fields_arrays, only: phi, apar, bpar
    use fields_arrays, only: phiold, aparold, bparold
    use dist_fn, only: getfieldexp
    use dist_fn_arrays, only: gold, g

    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    gold = g
    call getfieldexp (gold, phiold, aparold, bparold, 'g')
    phi = phiold; apar = aparold; bpar = bparold

  end subroutine init_phi_explicit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dg_scheme

    use gs2_time, only: code_dt
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo

    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Initialise adaptive timestep stuff

    if (adaptive_dt) then
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

    !  ne is the number of finite elements,
    !  each representing p finite difference grid points

    ne = (2*ntgrid+1)/p

    allocate( &
         v(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         wd(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc) )

  end subroutine init_dg_scheme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advance_explicit (istep)

    use collisions, only: solfp1
    use dist_fn, only: getfieldexp, g_adjust_exp
    use dist_fn, only: gnl_1, gnl_2, gnl_3, def_parity, even
    use dist_fn, only: wdrift, wcoriolis
    use dist_fn_arrays, only: g, gnew, gold, gwork, vpar
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_layouts, only: g_lo, ik_idx, is_idx, il_idx, ie_idx
    use gs2_time, only: code_time, code_dt
    use hyper, only: hyper_diff
    use le_grids, only: forbid
    use mp, only: max_allreduce
    use nonlinear_terms, only: add_nonlinear_terms
    use run_parameters, only: tunits, fphi, fapar, fbpar
    use species, only: spec
    use theta_grid, only: ntgrid

    implicit none

    !  Arguments

    integer, intent (in) :: istep

    !  Local variables

    real :: dt, dtnew
    integer :: nx, ig, ie, isgn, iglo, ik, is, il
    integer :: lb1, ub1, lb3, ub3
    integer :: modep
    complex, parameter :: z1 = (1.0, 0.0)
    complex, allocatable, dimension (:,:,:) :: gmodal, gnew1, gnew2

    integer, save :: first_call = 1

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    istep_dg = istep  !  to pass into dydt routine

    lb1 = -ntgrid ; ub1 = ntgrid
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc
    nx = (ub1-lb1+1)*2*(ub3-lb3+1)

    gold = g

!+PJK Following two lines are present in advance_implicit...
!    call antenna_amplitudes (apar_ext)
!    if (allocated(kx_shift)) call exb_shear (gnew, phinew, aparnew, bparnew) 
!-PJK

    write(65,*) 2*ntgrid+1

    allocate( gmodal(lb1:ub1,1:2,lb3:ub3), &
         gnew1(lb1:ub1,1:2,lb3:ub3), gnew2(lb1:ub1,1:2,lb3:ub3) )

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

    !  Evaluate fields at old timestep (still need to do this because of
    !  call to g_adjust_exp)

    call getfieldexp (g, phi, apar, bpar, 'g')
    call add_nonlinear_terms (gnl_1, gnl_2, gnl_3, &
         phi, apar, bpar, istep, 0.0, z1)

    !  Original call being replaced...
    !call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)

    call g_adjust_exp (g, apar, fapar)  !  convert g to i

!    !  On first call only, approximate the 'best' choice for the
!    !  initial time step
!
!    if (adaptive_dt) then
!       dt = code_dt
!       adaptive_dt_reset = .false.
!       if (first_call == 1) then
!          call adaptive_dt0(g,phi,apar,bpar,v,wd,istep,p,ne,code_time, &
!               dt0,dt,dtmin,dtmax,epsr,epsa)
!          if (dt /= code_dt) then
!             adaptive_dt_reset = .true.
!             adaptive_dt_new = dt
!          end if
!          first_call = 0
!       end if
!    end if

    gmodal = g
    call nodal2modal(gmodal,lb1,ub1,1,2,lb3,ub3)

    !  Time-advance using Runge-Kutta pair method

    call rk_advance3(gmodal,gnew1,gnew2,nx,code_time,code_dt)

    call modal2nodal(gnew1,lb1,ub1,1,2,lb3,ub3)
    call modal2nodal(gnew2,lb1,ub1,1,2,lb3,ub3)

    !  Boundary conditions...
!+PJK SHOULD THESE BE TURNED ON???
!+PJK or gnew(-ntgrid,1,iglo) = 0.0, gnew(ntgrid,2,iglo) = 0.0 ???
!    ynew1(-ntgrid,1,:) = y(-ntgrid,1,:)
!    ynew2(-ntgrid,1,:) = y(-ntgrid,1,:)
!    ynew1( ntgrid,2,:) = y( ntgrid,2,:)
!    ynew2( ntgrid,2,:) = y( ntgrid,2,:)

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

!do iglo = g_lo%llim_proc, g_lo%ulim_proc
!   il = il_idx(g_lo,iglo) ! il = ittp(ig) = 25 at ig = +-56, +-28, 0
!   ie = ie_idx(g_lo,iglo)
!   is = is_idx(g_lo,iglo)
!   if ((il == 14).and.(ie == 13).and.(is == 1)) then
!      do ig = -ntgrid,ntgrid
!         write(65,*) real(ig), real(g(ig,1,iglo))
!      end do
!   end if
!end do

    if (allocated(gnew1)) then
       deallocate(gmodal,gnew1,gnew2)
    end if

  end subroutine advance_explicit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine reset_init
    initialized = .false.
  end subroutine reset_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine adaptive_dt0(y,phi,apar,bpar,v,wd,istep,p,ne,t,dt0,dt,dtmin,dtmax,epsr,epsa)

    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use kt_grids, only: naky, ntheta0
    use mp, only: max_allreduce

    implicit none

    !  Arguments

    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: y  !  g_old (actually i_old)
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), &
         intent(in) :: phi, apar, bpar  !  fields
    real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: v
    real, dimension(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: wd
    integer, intent(in) :: istep
    integer, intent(in) :: p   !  order of the spatial Legendre fit
    integer, intent(in) :: ne  !  number of finite elements
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
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc

    allocate(fluxfn(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(src(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy(lb1:ub1,lb2:ub2,lb3:ub3))

    !  Evaluate the source terms

    call fluxfn_dg(dt0,phi,apar,bpar,istep,fluxfn,src)  !  All nodal

    !  Calculate dg/dt using the Discontinuous Galerkin scheme

    call dydt_dg(t,dt0,y,dy,p,ne,v,wd,fluxfn,src,lb1,ub1,lb3,ub3)

    call modal2nodal(dy,lb1,ub1,1,2,lb3,ub3)

    gdotn = maxval(abs(dy)) ; call max_allreduce(gdotn)
    gn = maxval(abs(y)) ; call max_allreduce(gn)

    epst = epsr*gn + epsa

    dtmax = (epst/gdotn)**(1.0/p)

    dt = min(dt0,dtmax)
    dt = max(dt,dtmin)

  end subroutine adaptive_dt0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rk_advance3(y,ynew1,ynew2,nx,t,dt)

    !  Rename to rk_adaptive_complex..

    use rk_schemes, rk => rkf3

    implicit none

    !  Arguments

    !  Two output values for ynew: one for each of the two orders of the
    !  chosen adaptive RK scheme

    integer, intent(in) :: nx
    complex, dimension(nx), intent(in) :: y
    complex, dimension(nx), intent(out) :: ynew1
    complex, dimension(nx), intent(out) :: ynew2
    real, intent(in) :: t  !  code_time
    real, intent(in) :: dt !  code_dt

    !  Local variables

    real :: ttest
    integer :: i, j
    complex, allocatable, dimension(:) :: aksum
    complex, allocatable, dimension(:,:) :: k
    complex, allocatable, dimension(:) :: ytest, ydot

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(ytest(1:nx))
    allocate(ydot(1:nx))
    allocate(aksum(1:nx))
    allocate(k(1:nx,rk%stages))

    !  Runge-Kutta loop

    do i = 1,rk%stages
       aksum = 0.0
       do j = 1,i-1
          aksum(:) = aksum(:) + rk%aij(i,j)*k(:,j)
       end do

       ytest = y + dt*aksum
       ttest = t + dt*rk%ci(i)

       !  Calculate dy/dt

       call dydt_dggs2(ytest,ttest,dt,ydot)
       k(:,i) = ydot
    end do

    ynew1 = y
    ynew2 = y

    do i = 1,rk%stages
       ynew1 = ynew1 + dt*rk%bi(1,i)*k(:,i)
       ynew2 = ynew2 + dt*rk%bi(2,i)*k(:,i)
    end do

  end subroutine rk_advance3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dydt_dggs2(gmodal,t,dt,dgdt_modal)

    use dist_fn, only: getfieldexp, get_source_term_exp
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use gs2_layouts, only: g_lo, ik_idx, il_idx
    use le_grids, only: nlambda, ng2, lmax, forbid
    use run_parameters, only: tunits
    use theta_grid, only: ntgrid

    implicit none

    !  Arguments

    real, intent(in) :: t, dt
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: gmodal
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: dgdt_modal

    !  Local variables

    integer :: ip, ip2
    integer :: i, j, kk, element, ig, ik, is, il, isgn, iglo
    integer :: lb1, ub1, lb3, ub3

    complex, parameter :: zi = (0.0,1.0)
    real :: ddt

    complex, allocatable, dimension(:,:,:) :: g, gnew, gnewh
    complex, allocatable, dimension(:,:,:) :: dgdt
    complex, allocatable, dimension(:,:,:) :: iwg, iwg_modal
    complex, allocatable, dimension(:,:,:) :: gplusf, gplusf_modal
    complex, allocatable, dimension(:,:,:) :: mgplusf, mgplusf_modal
    complex, allocatable, dimension(:,:,:) :: vterm, vterm_modal
    complex, allocatable, dimension(:,:,:) :: fluxfn
    complex, allocatable, dimension(:,:,:) :: src, smodal

    !real, dimension(2,4) :: mat_p2p, mat_p2m
    real, dimension(3,6) :: mat_p3p, mat_p3m

    complex, dimension(2*p) :: ftemp

    !  Constraints variables (to control dg/dt at certain points on an element)

    real, parameter :: sixth = 1.0/6.0
    real, parameter :: twothirds = 2.0/3.0
  
    complex, dimension(0:2) :: amodal, bmodal, cnodal
    complex, dimension(0:2,0:2) :: n2m, m2n

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  local matrices for nodal to modal, and vice versa

    n2m(0,0) = 0.375D0
    n2m(0,1) = 0.25D0
    n2m(0,2) = 0.375D0
    n2m(1,0) = -0.75D0
    n2m(1,1) = 0.0D0
    n2m(1,2) = 0.75D0
    n2m(2,0) = 0.75D0
    n2m(2,1) = -1.5D0
    n2m(2,2) = 0.75D0

    m2n(0,0) = 1.0D0
    m2n(0,1) = -twothirds
    m2n(0,2) = sixth
    m2n(1,0) = 1.0D0
    m2n(1,1) = 0.0D0
    m2n(1,2) = -0.5D0
    m2n(2,0) = 1.0D0
    m2n(2,1) = twothirds
    m2n(2,2) = sixth

    lb1 = -ntgrid ; ub1 = ntgrid
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc

    allocate( &
         g(lb1:ub1,2,lb3:ub3), &
         dgdt(lb1:ub1,2,lb3:ub3), &
         gnew(lb1:ub1,2,lb3:ub3), &
         iwg(lb1:ub1,2,lb3:ub3), &
         iwg_modal(lb1:ub1,2,lb3:ub3), &
         gplusf(lb1:ub1,2,lb3:ub3), &
         gplusf_modal(lb1:ub1,2,lb3:ub3), &
         mgplusf(lb1:ub1,2,lb3:ub3), &
         mgplusf_modal(lb1:ub1,2,lb3:ub3), &
         vterm(lb1:ub1,2,lb3:ub3), &
         vterm_modal(lb1:ub1,2,lb3:ub3), &
         fluxfn(lb1:ub1,2,lb3:ub3), &
         src(lb1:ub1,2,lb3:ub3), &
         smodal(lb1:ub1,2,lb3:ub3) )

    !  Obtain g in nodal form

    g = gmodal
    call modal2nodal(g,lb1,ub1,1,2,lb3,ub3)

    !  Evaluate field equations

    call getfieldexp(g,phitmp,apartmp,bpartmp,'i')

    !  Calculate source terms

    do iglo = lb3, ub3
       ik = ik_idx(g_lo,iglo)

       do isgn = 1,2
          call get_source_term_exp(phitmp,apartmp,bpartmp,istep_dg,isgn,iglo, &
               fluxfn(:,isgn,iglo), src(:,isgn,iglo))
          !  src is actually (2.dt.S), therefore need to divide by dt
          src(:,isgn,iglo) = src(:,isgn,iglo)/(2.0*dt*tunits(ik))
       end do
    end do

    !  Evaluate -i*wd*g term

    do iglo = lb3,ub3
       il = il_idx(g_lo,iglo)
       do isgn = 1,2
          do ig = lb1,ub1
             if (forbid(ig,il)) then
                iwg(ig,isgn,iglo) = 0.0
             else
                iwg(ig,isgn,iglo) = -zi * wd(ig,iglo) * g(ig,isgn,iglo) / dt
             end if
          end do
       end do
    end do
    iwg_modal = iwg
    call nodal2modal(iwg_modal,lb1,ub1,1,2,lb3,ub3)

    !  Add g to F

    gplusf = g + fluxfn
    gplusf_modal = gplusf
    call nodal2modal(gplusf_modal,lb1,ub1,1,2,lb3,ub3)

    !  Perform matrix multiplication to calculate d(g+F)/dz, using the correct matrix
    !  For speed, ought to store and use transpose...

    mgplusf_modal = 0.0

!    !  Matrix, p=2, positive v
!    mat_p2p = real(reshape( &
!         source = (/ 1,-3, 1,-3, -1,3, -1,-3 /), &
!         shape = (/ 2,4 /) ))
!    !  Matrix, p=2, negative v
!    mat_p2m = real(reshape( &
!         source = (/ 1,3, -1,3, -1,-3, 1,3 /), &
!         shape = (/ 2,4 /) ))

    !  Matrix, p=3, positive v
    mat_p3p = real(reshape( &
         source = (/ 1,-3,5, 1,-3,5, 1,-3,5, -1,3,-5, -1,-3,5, -1,-3,-5 /), &
         shape = (/ 3,6 /) ))
    !  Matrix, p=3, negative v
    mat_p3m = real(reshape( &
         source = (/ 1,3,5, -1,3,5, 1,-3,5, -1,-3,-5, 1,3,5, -1,-3,-5 /), &
         shape = (/ 3,6 /) ))

    do iglo = lb3,ub3
       do isgn = 1,2

          if (isgn == 1) then  !  positive parallel velocity
                               !  v(*,isgn==1,*) >= 0.0, <= 0.0 for isgn==2
                               !  Checked to be true.

             do element = 1,ne
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

             do element = 1,ne
                !  Fill ftemp with correct elements from the flux function
                kk = lb1 + p*(element-1)  !  take values from this element
                ftemp(1:p) = gplusf_modal(kk:kk+p-1,isgn,iglo)
                if (element /= ne) then
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

    !  Multiply this term by v/dz; must be done in nodal form
    !  Denominator dz is element size, not finite difference mesh separation,
    !  so need to divide by p.
    !  v is non-dimensional, and contains a factor dt/(F.D. dz), so we also have
    !  to divide by dt

    mgplusf = mgplusf_modal
    call modal2nodal(mgplusf,lb1,ub1,1,2,lb3,ub3)

    do iglo = lb3,ub3 
       il = il_idx(g_lo,iglo)
       do isgn = 1,2
          do ig = lb1,ub1
             if (forbid(ig,il)) then
                vterm(ig,isgn,iglo) = 0.0
             else
                vterm(ig,isgn,iglo) = 1.0/(p*dt) * v(ig,isgn,iglo) &
                     * mgplusf(ig,isgn,iglo)
             end if
          end do
       end do
    end do
    vterm_modal = vterm
    call nodal2modal(vterm_modal,lb1,ub1,1,2,lb3,ub3)

    !  Source term

    smodal = src
    call nodal2modal(smodal,lb1,ub1,1,2,lb3,ub3)

    !  Finally sum all the modal terms to get dg/dt (modal)

    do iglo = lb3,ub3
       do isgn = 1,2
          do i = lb1,ub1 ! modes
             dgdt_modal(i,isgn,iglo) = &
                  iwg_modal(i,isgn,iglo) &
                  + vterm_modal(i,isgn,iglo) &
                  + smodal(i,isgn,iglo)
          end do
       end do
    end do

    !  Now deal with boundary conditions...
    !  Passing particles only

    !  v >= 0

    do iglo = lb3,ub3
       amodal(0:2) = dgdt_modal(lb1:lb1+2,1,iglo)
       cnodal = matmul(m2n,amodal)
       cnodal(0) = 0.0 - cnodal(0)  !  Want left-hand nodal point to be zero
       cnodal(1) = 0.0
       cnodal(2) = 0.0
       bmodal = matmul(n2m,cnodal)
       dgdt_modal(lb1:lb1+2,1,iglo) = amodal + bmodal
    end do

    !  v < 0

    do iglo = lb3,ub3
       amodal(0:2) = dgdt_modal(ub1-2:ub1,2,iglo)
       cnodal = matmul(m2n,amodal)
       cnodal(0) = 0.0
       cnodal(1) = 0.0
       cnodal(2) = 0.0 - cnodal(2)  !  Want right-hand nodal point to be zero
       bmodal = matmul(n2m,cnodal)
       dgdt_modal(ub1-2:ub1,2,iglo) = amodal + bmodal
    end do

!    dgdt = dgdt_modal
!    call modal2nodal(dgdt,lb1,ub1,1,2,lb3,ub3)

    !  No trapped particle case
    !  (nlambda > ng2) == .true. means 'There ARE trapped particles'

    if (nlambda <= ng2) then

       !  Enforce no change in g at each end of theta domain
!       dgdt(lb1,:,:) = 0.0
!       dgdt(ub1,:,:) = 0.0
!
!       dgdt_modal = dgdt
!       call nodal2modal(dgdt_modal,lb1,ub1,1,2,lb3,ub3)

    else
       stop
    end if
    !  We now have homogeneous and inhomogeneous solutions for both
    !  2nd order and 3rd order methods.
    !  Now, for trapped particles, we need to sum these contributions with a
    !  factor beta to ensure that ynew(v>0) = ynew(v<0) at each bounce point

    ! do iglo = g_lo%llim_proc,g_lo%ulim_proc
    !    il = il_idx(g_lo,iglo)

    !    !  (nlambda > ng2) == .true. means 'There ARE trapped particles'
    !    !  (il >= ng2+2 .and. il <= lmax) == .true.
    !    !     means 'This is a trapped particle orbit'
    !    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
    !       ynew1(:,:,iglo) = 0.0
    !       lbp = -9999 ; ubp = -9999
    !       do ig = -ntgrid, ntgrid-1
    !          if (forbid(ig,il).and..not.forbid(ig+1,il)) then
    !             lbp = ig+1  !  lower bounce point found
    !          end if
    !          if (.not.forbid(ig,il).and.forbid(ig+1,il)) then
    !             ubp = ig  !  upper bounce point found

    !             if (lbp == -9999) then
    !                write(*,*) 'oh dear... UBP found but no LBP!'
    !                write(*,*) lbp,ubp,il
    !                stop
    !             end if

    !             !  Add correct amount of homogeneous solution in the region
    !             !  between LBP and UBP

    !             beta1 = (ynew1i(lbp,2,iglo) + ynew1i(ubp,2,iglo) - &
    !                  ynew1i(lbp,1,iglo) - ynew1i(ubp,1,iglo)) / &
    !                  (ynew1h(lbp,1,iglo) + ynew1h(ubp,1,iglo) - &
    !                  ynew1h(lbp,2,iglo) - ynew1h(ubp,2,iglo))

    !             do igg = lbp,ubp
    !                ynew1(igg,:,iglo) = ynew1i(igg,:,iglo) + beta1*ynew1h(igg,:,iglo)
    !             end do
    !             write(*,*) 'LBP: ',lbp,ynew1(lbp,1,iglo)-ynew1(lbp,2,iglo)
    !             write(*,*) 'UBP: ',ubp,ynew1(ubp,1,iglo)-ynew1(ubp,2,iglo)
    !             lbp = -9999 ; ubp = -9999
    !          end if
    !       end do
    !    else  !  nothing trapped at this il
    !       ynew1(:,:,iglo) = ynew1i(:,:,iglo)
    !    end if
    ! end do

  end subroutine dydt_dggs2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rk_advance2(y,ynew1,ynew2,phi,apar,bpar,v,wd,istep,p,ne,t,h)

    use rk_schemes, rk => rkf3
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use dist_fn, only: getfieldexp
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, il_idx
    use run_parameters, only: tunits
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, ng2, lmax, forbid

    implicit none

    !  Arguments

    !  Two output values for ynew: one for each of the two orders of the
    !  chosen adaptive RK scheme

    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: y  !  g_old (actually i_old)
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: ynew1
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: ynew2
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), &
         intent(in) :: phi, apar, bpar  !  fields
    real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: v
    real, dimension(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: wd
    integer, intent(in) :: istep
    integer, intent(in) :: p   !  order of the spatial Legendre fit
    integer, intent(in) :: ne  !  number of finite elements
    real, intent(in) :: t  !  code_time
    real, intent(in) :: h  !  code_dt

    !  Local variables

    integer :: ix, i, j, ie, ik, il, iglo, ig, igg
    integer :: lb1, ub1, lb2, ub2, lb3, ub3, lbp, ubp
    complex :: beta1, beta2
    complex, allocatable, dimension(:,:,:) :: fluxfn, f0
    complex, allocatable, dimension(:,:,:) :: src, src0
    complex, allocatable, dimension(:,:,:) :: dy, dy_modal
    complex, allocatable, dimension(:,:,:) :: aksum, aksumh
    complex, allocatable, dimension(:,:,:,:) :: k, kh
    complex, allocatable, dimension(:,:,:) :: ytmp, ymodal, yh, yh0
    complex, allocatable, dimension(:,:,:) :: ynew1_modal, ynew2_modal
    complex, allocatable, dimension(:,:,:) :: ynew1h_modal, ynew2h_modal
    complex, allocatable, dimension(:,:,:) :: ynew1i, ynew2i
    complex, allocatable, dimension(:,:,:) :: ynew1h, ynew2h

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lb1 = -ntgrid ; ub1 = ntgrid
    lb2 = 1 ; ub2 = 2
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc

    allocate(ytmp(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(yh(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(yh0(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ymodal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew1_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew1h_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew2_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew2h_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(fluxfn(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(f0(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(src(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(src0(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy_modal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(aksum(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(aksumh(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(k(rk%stages,lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(kh(rk%stages,lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew1i(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew2i(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew1h(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(ynew2h(lb1:ub1,lb2:ub2,lb3:ub3))

    !  Storage for homogeneous solution

    yh0 = 1.0
    f0 = 0.0
    src0 = 0.0

    !  Evaluate the source terms

    call fluxfn_dg(h,phi,apar,bpar,istep,fluxfn,src)  !  All nodal

    !  Calculate dg/dt using the Discontinuous Galerkin scheme

    call dydt_dg(t,h,y,dy_modal,p,ne,v,wd,fluxfn,src,lb1,ub1,lb3,ub3)
    k(1,:,:,:) = dy_modal

    call dydt_dg(t,h,y,dy_modal,p,ne,v,wd,f0,src0,lb1,ub1,lb3,ub3)
    kh(1,:,:,:) = dy_modal

    !  Runge-Kutta loop

    do i = 2,rk%stages
       aksum = 0.0
       aksumh = 0.0
       do j = 1,i-1
          aksum(:,:,:) = aksum(:,:,:) + rk%aij(i,j)*k(j,:,:,:)
          aksumh(:,:,:) = aksumh(:,:,:) + rk%aij(i,j)*kh(j,:,:,:)
       end do

       dy = h*aksum
       call modal2nodal(dy,lb1,ub1,1,2,lb3,ub3)
       ytmp = y + dy

       dy = h*aksumh
       call modal2nodal(dy,lb1,ub1,1,2,lb3,ub3)
       yh = yh0 + dy

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

       call dydt_dg(t+h*rk%ci(i),h,ytmp,dy_modal,p,ne,v,wd,fluxfn,src, &
            lb1,ub1,lb3,ub3)
       k(i,:,:,:) = dy_modal

       call dydt_dg(t+h*rk%ci(i),h,yh,dy_modal,p,ne,v,wd,f0,src0, &
            lb1,ub1,lb3,ub3)
       kh(i,:,:,:) = dy_modal
    end do

    ymodal = y
    call nodal2modal(ymodal,lb1,ub1,1,2,lb3,ub3)
    ynew1_modal = ymodal
    ynew2_modal = ymodal

    ymodal = yh0
    call nodal2modal(ymodal,lb1,ub1,1,2,lb3,ub3)
    ynew1h_modal = ymodal
    ynew2h_modal = ymodal

    do iglo = g_lo%llim_proc,g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       do i = 1,rk%stages
          ynew1_modal(:,:,iglo) = ynew1_modal(:,:,iglo) + &
               h*tunits(ik)*rk%bi(1,i)*k(i,:,:,iglo)
          ynew2_modal(:,:,iglo) = ynew2_modal(:,:,iglo) + &
               h*tunits(ik)*rk%bi(2,i)*k(i,:,:,iglo)

          ynew1h_modal(:,:,iglo) = ynew1h_modal(:,:,iglo) + &
               h*tunits(ik)*rk%bi(1,i)*kh(i,:,:,iglo)
          ynew2h_modal(:,:,iglo) = ynew2h_modal(:,:,iglo) + &
               h*tunits(ik)*rk%bi(2,i)*kh(i,:,:,iglo)
!  Non-adaptive RK types only...
!          ynew1_modal(:,:,iglo) = ynew1_modal(:,:,iglo) + &
!               h*tunits(ik)*rk%bi(i)*k(i,:,:,iglo)
!          ynew2_modal = ynew1_modal
       end do
    end do

    ynew1i = ynew1_modal
    call modal2nodal(ynew1i,lb1,ub1,1,2,lb3,ub3)
    ynew2i = ynew2_modal
    call modal2nodal(ynew2i,lb1,ub1,1,2,lb3,ub3)

    ynew1h = ynew1h_modal
    call modal2nodal(ynew1h,lb1,ub1,1,2,lb3,ub3)
    ynew2h = ynew2h_modal
    call modal2nodal(ynew2h,lb1,ub1,1,2,lb3,ub3)

    !  We now have homogeneous and inhomogeneous solutions for both
    !  2nd order and 3rd order methods.
    !  Now, for trapped particles, we need to sum these contributions with a
    !  factor beta to ensure that ynew(v>0) = ynew(v<0) at each bounce point

    do iglo = g_lo%llim_proc,g_lo%ulim_proc
       il = il_idx(g_lo,iglo)

       !  (nlambda > ng2) == .true. means 'There ARE trapped particles'
       !  (il >= ng2+2 .and. il <= lmax) == .true.
       !     means 'This is a trapped particle orbit'
       if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
          ynew1(:,:,iglo) = 0.0
          ynew2(:,:,iglo) = 0.0
          lbp = -9999 ; ubp = -9999
          do ig = -ntgrid, ntgrid-1
             if (forbid(ig,il).and..not.forbid(ig+1,il)) then
                lbp = ig+1  !  lower bounce point found
             end if
             if (.not.forbid(ig,il).and.forbid(ig+1,il)) then
                ubp = ig  !  upper bounce point found

                if (lbp == -9999) then
                   write(*,*) 'oh dear... UBP found but no LBP!'
                   write(*,*) lbp,ubp,il
                   stop
                end if

                !  Add correct amount of homogeneous solution in the region
                !  between LBP and UBP

                beta1 = (ynew1i(lbp,2,iglo) + ynew1i(ubp,2,iglo) - &
                     ynew1i(lbp,1,iglo) - ynew1i(ubp,1,iglo)) / &
                     (ynew1h(lbp,1,iglo) + ynew1h(ubp,1,iglo) - &
                     ynew1h(lbp,2,iglo) - ynew1h(ubp,2,iglo))

                beta2 = (ynew2i(lbp,2,iglo) + ynew2i(ubp,2,iglo) - &
                     ynew2i(lbp,1,iglo) - ynew2i(ubp,1,iglo)) / &
                     (ynew2h(lbp,1,iglo) + ynew2h(ubp,1,iglo) - &
                     ynew2h(lbp,2,iglo) - ynew2h(ubp,2,iglo))

                do igg = lbp,ubp
                   ynew1(igg,:,iglo) = ynew1i(igg,:,iglo) + beta1*ynew1h(igg,:,iglo)
                   ynew2(igg,:,iglo) = ynew2i(igg,:,iglo) + beta2*ynew2h(igg,:,iglo)
                end do
                write(*,*) 'LBP: ',lbp,ynew1(lbp,1,iglo)-ynew1(lbp,2,iglo)
                write(*,*) 'UBP: ',ubp,ynew1(ubp,1,iglo)-ynew1(ubp,2,iglo)
                lbp = -9999 ; ubp = -9999
             end if
          end do
       else  !  nothing trapped at this il
          ynew1(:,:,iglo) = ynew1i(:,:,iglo)
          ynew2(:,:,iglo) = ynew2i(:,:,iglo)
       end if
    end do

    !  Boundary conditions...
!+PJK SHOULD THESE BE TURNED ON???
!+PJK or gnew(-ntgrid,1,iglo) = 0.0, gnew(ntgrid,2,iglo) = 0.0 ???
!    ynew1(-ntgrid,1,:) = y(-ntgrid,1,:)
!    ynew2(-ntgrid,1,:) = y(-ntgrid,1,:)
!    ynew1( ntgrid,2,:) = y( ntgrid,2,:)
!    ynew2( ntgrid,2,:) = y( ntgrid,2,:)

  end subroutine rk_advance2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dydt_dg(t,dt,g,dgdt,p,ne,v,wd,f,source,lb1,ub1,lb3,ub3)

    use gs2_layouts, only: g_lo, ik_idx, il_idx
    use le_grids, only: forbid

    implicit none

    integer, intent(in) :: p, ne
    integer, intent(in) :: lb1,ub1,lb3,ub3
    real, intent(in) :: t, dt
    real, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: v
    real, dimension(lb1:ub1,lb3:ub3), intent(in) :: wd
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: g
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: f
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(in) :: source
    complex, dimension(lb1:ub1,2,lb3:ub3), intent(out) :: dgdt  !  modal

    integer :: ip, ip2
    integer :: i, j, kk, element, ig, ik, is, il, isgn, iglo

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

! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!
! move trapped stuff into this routine !!!!!

    !  Evaluate -i*wd*g term

    do iglo = lb3,ub3
       il = il_idx(g_lo,iglo)
       do isgn = 1,2
          do ig = lb1,ub1
             if (forbid(ig,il)) then
                iwg(ig,isgn,iglo) = 0.0
             else
                iwg(ig,isgn,iglo) = -zi * wd(ig,iglo) * g(ig,isgn,iglo) / dt
             end if
          end do
       end do
    end do
    iwg_modal = iwg
    call nodal2modal(iwg_modal,lb1,ub1,1,2,lb3,ub3)

    !  Add g to F

    gplusf = g + f
    gplusf_modal = gplusf
    call nodal2modal(gplusf_modal,lb1,ub1,1,2,lb3,ub3)

    !  Perform matrix multiplication, using the correct matrix above
    !  For speed, ought to store and use transpose...

    mgplusf_modal = 0.0

    if (p == 2) then

       ! do iglo = lb3,ub3
       !    do i = 1,2

       !       if (i == 1) then  !  positive parallel velocity

       !          do element = 1,ne
       !             !  Fill ftemp with correct elements from the flux function
       !             if (element /= 1) then
       !                kk = lb1 + p*(element-2)  !  take values from element to the left
       !                ftemp(1:p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             else
       !                kk = lb1 + p*(ne-1)  !  assumes periodic BCs for now
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

       !          do element = 1,ne
       !             !  Fill ftemp with correct elements from the flux function
       !             kk = lb1 + p*(element-1)  !  take values from this element
       !             ftemp(1:p) = gplusf_modal(kk:kk+p-1,i,iglo)
       !             if (element /= ne) then
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

                do element = 1,ne
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

                do element = 1,ne
                   !  Fill ftemp with correct elements from the flux function
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(1:p) = gplusf_modal(kk:kk+p-1,isgn,iglo)
                   if (element /= ne) then
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
    call modal2nodal(mgplusf,lb1,ub1,1,2,lb3,ub3)

    do iglo = lb3,ub3 
       il = il_idx(g_lo,iglo)
       do isgn = 1,2
          do ig = lb1,ub1
             if (forbid(ig,il)) then
                vterm(ig,isgn,iglo) = 0.0
             else
                vterm(ig,isgn,iglo) = 1.0/(p*dt) * v(ig,isgn,iglo) &
                     * mgplusf(ig,isgn,iglo)
             end if
          end do
       end do
    end do
    vterm_modal = vterm
    call nodal2modal(vterm_modal,lb1,ub1,1,2,lb3,ub3)

    !  Source term

    smodal = source
    call nodal2modal(smodal,lb1,ub1,1,2,lb3,ub3)

    !  Finally sum all the modal terms to get dg/dt (modal)

    do iglo = lb3,ub3
       do isgn = 1,2
          do i = lb1,ub1 ! modes
             dgdt(i,isgn,iglo) = &
                  iwg_modal(i,isgn,iglo) &
                  + vterm_modal(i,isgn,iglo) &
                  + smodal(i,isgn,iglo)
          end do
       end do
    end do

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
    complex, dimension(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: f  !  flux function
    complex, dimension(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc), &
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

  subroutine pjk_advance(y,ynew1,ynew2,phi,apar,bpar,v,wd,istep,p,ne,t,h)

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

    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: y  !  g_old (actually i_old)
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: ynew1
    complex, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: ynew2
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), &
         intent(in) :: phi, apar, bpar  !  fields
    real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: v
    real, dimension(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: wd
    integer, intent(in) :: istep
    integer, intent(in) :: p   !  order of the spatial Legendre fit
    integer, intent(in) :: ne  !  number of finite elements
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
