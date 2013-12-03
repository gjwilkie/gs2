!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        1, reshape( (/0.0/),(/1,1/) ), (/1.0/), (/0.0/) )

   type(rk2_scheme), parameter :: midpoint = rk2_scheme( &
        2, &
        reshape( (/ 0.0, 0.5, 0.0, 0.0 /), (/2,2/) ), &
        (/ 0.0, 1.0 /), &
        (/ 0.0, 0.5 /) )
   type(rk2_scheme), parameter :: heun = rk2_scheme( &
        2, &
        reshape( (/ 0.0, 1.0, 0.0, 0.0 /), (/2,2/) ), &
        (/ 0.5, 0.5 /), &
        (/ 0.0, 1.0 /) )
   type(rk2_scheme), parameter :: atkinson = rk2_scheme( &
        2, &
        reshape( (/ 0.0, twothirds, 0.0, 0.0 /), (/2,2/) ), &
        (/ 0.25, 0.75 /), &
        (/ 0.0, twothirds /) )

   type(rk4_scheme), parameter :: rk4 = rk4_scheme( & !  classical Runge-Kutta method
        4, &
        reshape( (/ 0.0, 0.5, 0.0, 0.0, &
        0.0, 0.0, 0.5, 0.0, &
        0.0, 0.0, 0.0, 1.0, &
        0.0, 0.0, 0.0, 0.0 /), (/4,4/) ), &
        (/ sixth, third, third, sixth /), &
        (/ 0.0, 0.5, 0.5, 1.0/) )

   type(rk12_scheme), parameter :: simple_adaptive = rk12_scheme( &
        2, &
        reshape( (/ 0.0, 1.0, &
        0.0, 0.0 /), (/2,2/) ), &
        reshape( (/ 1.0, 0.5, &
        0.0, 0.5 /), (/2,2/) ), &
        (/ 0.0, 1.0 /) )

   type(rk23_scheme), parameter :: bogacki_shampine = rk23_scheme( &
        4, &
        reshape( (/ 0.0, 0.5, 0.0, 2.0/9.0, &
        0.0, 0.0, 0.75, third, &
        0.0, 0.0, 0.0, 4.0/9.0, &
        0.0, 0.0, 0.0, 0.0 /), (/4,4/) ), &
        reshape( (/ 2.0/9.0, 7.0/24.0, &
        third, 0.25, &
        4.0/9.0, third, &
        0.0, 0.125 /), (/2,4/) ), &
        (/ 0.0, 0.5, 0.75, 1.0 /) )

   type(rk23_scheme), parameter :: rkf3 = rk23_scheme( & !  Stoer & Bulirsch
        4, &
        reshape( (/ 0.0, 0.25, -189.0/800.0, 214.0/891.0, &
        0.0, 0.0, 729.0/800.0, 1.0/33.0, &
        0.0, 0.0, 0.0, 650.0/891.0, &
        0.0, 0.0, 0.0, 0.0 /), (/4,4/) ), &
        reshape( (/ 214.0/891.0, 533.0/2106.0, &
        1.0/33.0, 0.0, &
        650.0/891.0, 800.0/1053.0, &
        0.0, -1.0/78.0 /), (/2,4/) ), &
        (/ 0.0, 0.0, 0.0, 0.0 /) )

   type(rk45_scheme), parameter :: fehlberg = rk45_scheme( &
        6, &
        reshape( (/ 0.0, 0.25, 0.09375, 1932.0/2197.0, 439.0/216.0, -8.0/27.0, &
        0.0, 0.0, 0.28125, -7200.0/2197.0, -8.0, 2.0, &
        0.0, 0.0, 0.0, 7296.0/2197.0, 3680.0/513.0, -3544.0/2565.0, &
        0.0, 0.0, 0.0, 0.0, -845.0/4104.0, 1859.0/4104.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, -0.275, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), (/6,6/) ), &
        reshape( (/ 25.0/216.0, 16.0/135.0, &
        0.0, 0.0, &
        1408.0/2565.0, 6656.0/12825.0, &
        2197.0/4104.0, 28561.0/56430.0, &
        -0.2, -0.18, &
        0.0, 2.0/55.0 /), (/2,6/) ), &
        (/ 0.0, 0.25, 0.375, 12.0/13.0, 1.0, 0.5 /) )

   type(rk45_scheme), parameter :: cash_karp = rk45_scheme( &
        6, &
        reshape( (/ 0.0, 0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, &
        0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875, &
        0.0, 0.0, 0.0, 1.2, -70.0/27.0, 575.0/13824.0, &
        0.0, 0.0, 0.0, 0.0, 35.0/27.0, 44275.0/110592.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 253.0/4096.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), (/6,6/) ), &
        reshape( (/ 37.0/378.0, 2825.0/27648.0, &
        0.0, 0.0, &
        250.0/621.0, 18575.0/48384.0, &
        125.0/594.0, 13525.0/55296.0, &
        0.0, 277.0/14336.0, &
        512.0/1771.0, 0.25 /), (/2,6/) ), &
        (/ 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /) )

   type(rk45a_scheme), parameter :: dormand_prince = rk45a_scheme( &
        7, &
        reshape( (/ 0.0, 0.2, 0.075, 44.0/45.0, 19372.0/6561.0, 9017.0/3168.0, 35.0/384.0, &
        0.0, 0.0, 0.225, -56.0/15.0, -25360.0/2187.0, -355.0/33.0, 0.0, &
        0.0, 0.0, 0.0, 32.0/9.0, 64448.0/6561.0, 46732.0/5247.0, 500.0/1113.0, &
        0.0, 0.0, 0.0, 0.0, -212.0/729.0, 49.0/176.0, 125.0/192.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, -5103.0/18656.0, -2187.0/6784.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0/84.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), (/7,7/) ), &
        reshape( (/ 5179.0/57600.0, 35.0/384.0, &
        0.0, 0.0, &
        7571.0/16695.0, 500.0/1113.0, &
        393.0/640.0, 125.0/192.0, &
        -92097.0/339200.0, -2187.0/6784.0, &
        187.0/2100.0, 11.0/84.0, &
        0.025, 0.0 /), (/2,7/) ), &
        (/ 0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0 /) )

end module rk_schemes

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dg_scheme

  implicit none

  !  Variables and routines required for use with the Discontinuous Galerkin scheme

  public :: adaptive_dt, adaptive_dt_reset, adaptive_dt_new

  interface nodal2modal
     module procedure nodal2modal_full
     module procedure nodal2modal_1d
  end interface nodal2modal

  interface modal2nodal
     module procedure modal2nodal_full
     module procedure modal2nodal_1d
  end interface modal2nodal

  !  Primary DG variables

  integer, parameter :: p = 3  !  scheme order, = number of FD grid points
                               !  per finite element
  integer :: ne     !  number of finite elements
  integer :: ne2    !  ne, or ne+2 if flux-tube (to account for FEs 0 and ne+1)
  integer :: ntgriddgmin  !  array lower bound in theta direction, for *some* arrays
  integer :: ntgriddg  !  array upper bound in theta direction

  !  Normalised parallel velocity, omega drift

  real, allocatable, dimension (:,:,:) :: v
  !+PJK  21/08/2013  New isgn dimension
  real, allocatable, dimension (:,:,:) :: wd
  !-PJK

  !  Flux tube scenario flag

  logical :: flux_tube = .false.

  !  Other stuff needed to be passed to the low-level routines

  integer :: istep_dg

  !  Adaptive timestep variables
!=PJK TURN ON ADAPTIVE_DT WHEN gs2_main, gs2_reinit HAVE BEEN MODIFIED
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

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nodal2modal_full(y,lb1,ub1,lb2,ub2,lb3,ub3,ne)

    !  Version for y(ig,isgn,iglo)
    !  Scheme order p = 3 is assumed.

    implicit none

    !  Arguments

    integer, intent(in) :: lb1, ub1, lb2, ub2, lb3, ub3
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y
    integer, intent(in) :: ne

    !  Local variables

    integer :: i, j, k, kk, element
    complex, dimension(ne*p) :: ytemp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

             y(kk,i,j)   = 0.125*(3.0*ytemp(k) + 2.0*ytemp(k+1) + 3.0*ytemp(k+2))
             y(kk+1,i,j) = 0.75*(ytemp(k+2) - ytemp(k))
             y(kk+2,i,j) = 0.75*(ytemp(k) - 2.0*ytemp(k+1) + ytemp(k+2))

          end do

       end do
    end do

  end subroutine nodal2modal_full

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nodal2modal_1d(y,lbound,ubound,ne)

    !  Version for y(ig), i.e. no isgn or iglo dimensions
    !  Scheme order p = 3 is assumed.

    implicit none

    !  Arguments

    integer, intent(in) :: lbound, ubound
    complex, dimension(lbound:ubound), intent(inout) :: y
    integer, intent(in) :: ne

    !  Local variables

    integer :: k, kk, element
    complex, dimension(ne*p) :: ytemp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  There are ne finite elements, each with p points

    ytemp = y(:)

    do element = 1,ne

       !  k is the index of first point of this finite element in ytemp;
       !  ytemp starts at index number 1.  kk is the equivalent index in y
       !  since this does not necessarily start at 1
       k = p*(element-1) + 1
       kk = lbound + k-1

       y(kk)   = 0.125*(3.0*ytemp(k) + 2.0*ytemp(k+1) + 3.0*ytemp(k+2))
       y(kk+1) = 0.75*(ytemp(k+2) - ytemp(k))
       y(kk+2) = 0.75*(ytemp(k) - 2.0*ytemp(k+1) + ytemp(k+2))

    end do

  end subroutine nodal2modal_1d

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine modal2nodal_full(y,lb1,ub1,lb2,ub2,lb3,ub3,ne)

    !  Version for y(ig,isgn,iglo)
    !  Scheme order p = 3 is assumed.

    implicit none

    !  Arguments

    integer, intent(in) :: lb1,ub1,lb2,ub2,lb3,ub3
    complex, dimension(lb1:ub1,lb2:ub2,lb3:ub3), intent(inout) :: y
    integer, intent(in) :: ne

    !  Local variables

    real, parameter :: twothirds = 2.0/3.0
    real, parameter :: sixth = 1.0/6.0
    integer :: i, j, k, kk, element
    complex, dimension(ne*p) :: ytemp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

             y(kk,i,j)   = ytemp(k) - twothirds*ytemp(k+1) + sixth*ytemp(k+2)
             y(kk+1,i,j) = ytemp(k)                        -   0.5*ytemp(k+2)
             y(kk+2,i,j) = ytemp(k) + twothirds*ytemp(k+1) + sixth*ytemp(k+2)

          end do

       end do
    end do

  end subroutine modal2nodal_full

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine modal2nodal_1d(y,lbound,ubound,ne)

    !  Version for y(ig), i.e. no isgn or iglo dimensions
    !  Scheme order p = 3 is assumed.

    implicit none

    !  Arguments

    integer, intent(in) :: lbound,ubound
    complex, dimension(lbound:ubound), intent(inout) :: y
    integer, intent(in) :: ne

    !  Local variables

    real, parameter :: twothirds = 2.0/3.0
    real, parameter :: sixth = 1.0/6.0
    integer :: k, kk, element
    complex, dimension(ne*p) :: ytemp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ytemp = y(:)

    do element = 1,ne

       !  k is the index of first point of this finite element in ytemp;
       !  ytemp starts at index number 1.  kk is the equivalent index in y
       !  since this does not necessarily start at 1
       k = p*(element-1) + 1
       kk = lbound + k-1

       y(kk)   = ytemp(k) - twothirds*ytemp(k+1) + sixth*ytemp(k+2)
       y(kk+1) = ytemp(k)                        -   0.5*ytemp(k+2)
       y(kk+2) = ytemp(k) + twothirds*ytemp(k+1) + sixth*ytemp(k+2)

    end do

  end subroutine modal2nodal_1d

end module dg_scheme

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fields_explicit

  use dg_scheme

  implicit none

  public :: advance_explicit
  public :: init_fields_explicit
  public :: init_phi_explicit
  public :: reset_init

  private

  logical :: initialized = .false.

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_fields_explicit

    use antenna, only: init_antenna
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts

    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts
    call init_theta_grid
    call init_kt_grids
    call init_antenna

    call init_dg_scheme

  end subroutine init_fields_explicit

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_phi_explicit

    use fields_arrays, only: phi, apar, bpar
    use dist_fn, only: getfieldexp
    use dist_fn_arrays, only: gold, g

    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    gold = g
    call getfieldexp (gold, phi, apar, bpar, 'g')

  end subroutine init_phi_explicit

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dg_scheme

    use gs2_time, only: code_dt
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use dist_fn, only: boundary_option_switch, boundary_option_linked

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

    !  Local flag denoting whether this is a flux-tube run
    !  rather than a ballooning run

    if (boundary_option_switch == boundary_option_linked) then
       flux_tube = .true.
    else
       flux_tube = .false.
    end if

    !  Constrain theta grid size to certain values for convenience

    if (flux_tube) then
       if (mod(2*ntgrid,p) /= 0) then
          write(*,*) 'Number of theta points is currently incompatible with'
          write(*,*) 'explicit DG scheme in flux tube mode.'
          write(*,*) 'Adjust ntheta to be a multiple of ',2*p
          write(*,*) 'or 1 + (multiple of ',2*p,')'
          stop
       end if
    end if

    !  ne is the number of finite elements,
    !  each representing p finite difference grid points.
    !  ne*p is guaranteed to cover at least (-ntgrid:ntgrid) grid points

    if (.not.flux_tube) then
       ne = ceiling( real(2*ntgrid+1)/p )
       ne2 = ne
    else
       !  ne*p must exactly cover the range -ntgrid to ntgrid-1 inclusive
       !  (assuming the mod(...) stuff above applies)
       ne = int( real(2*ntgrid+1)/p )
       ne2 = ne+2
    end if

    !  ntgriddg is the new upper bound for arrays (local to this module)
    !  in the theta direction, guaranteeing that exactly p grid points
    !  lie on each finite element
    !
    !  ntgriddgmin is the lower bound for certain arrays in the theta direction
    !  to include a 'halo' finite element for flux tube message passing

    if (.not.flux_tube) then
       ntgriddgmin = -ntgrid
       ntgriddg = ne*p - ntgrid - 1
    else
       ntgriddgmin = -ntgrid-p
       ntgriddg = (ne+1)*p - ntgrid - 1  !  as ne is one less with FTs
    end if

    if (.not.allocated(v)) then
       allocate( &
            v(-ntgrid:ntgriddg,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
            !+PJK  21/08/2013  New isgn dimension
            wd(-ntgrid:ntgriddg,1:2,g_lo%llim_proc:g_lo%ulim_proc) )
            !-PJK
    end if

  end subroutine init_dg_scheme

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine reset_init
    initialized = .false.
  end subroutine reset_init

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advance_explicit (istep)

    use collisions, only: solfp1
    use dist_fn, only: getfieldexp, g_adjust_exp
    use dist_fn, only: gexp_1, gexp_2, gexp_3, def_parity, even
!+PJK  21/08/2013
    use dist_fn, only: g0, parity_redist, parity_redist_ik, enforce_parity
#ifdef LOWFLOW
    use dist_fn, only: wdfac
#else
    use dist_fn, only: wdrift
#endif
!-PJK
    !+PJK  21/08/2013  use dist_fn_arrays, only: g, gnew, gold, gwork, vpa
    use dist_fn_arrays, only: g, gnew, gold, vpa, g_fixpar
    !-PJK
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use gs2_layouts, only: g_lo, ik_idx, is_idx, il_idx, ie_idx
    use gs2_time, only: code_time, code_dt
    use hyper, only: hyper_diff
    !+PJK  21/08/2013
    use le_derivatives, only: vspace_derivatives
    !-PJK
    use le_grids, only: nlambda, ng2, lmax, forbid
    use mp, only: max_allreduce
    !+PJK  21/08/2013  use nonlinear_terms, only: add_nonlinear_terms
    use nonlinear_terms, only: add_explicit_terms
    !-PJK
    use run_parameters, only: tunits, fphi, fapar, fbpar
    !+PJK  21/08/2013
    use run_parameters, only: fixpar_secondary
    !-PJK
    use species, only: spec
    use theta_grid, only: ntgrid, delthet, gradpar

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
    complex, allocatable, dimension (:,:,:) :: g_dg  ! local copy of g
                                                     ! with extended theta range
    integer, save :: first_call = 1

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    istep_dg = istep  !  to pass into dydt routine

    lb1 = ntgriddgmin ; ub1 = ntgriddg
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc
    nx = (ub1-lb1+1)*2*(ub3-lb3+1)

    gold = g

    !+PJK Following two lines are present in advance_implicit...
    !    call antenna_amplitudes (apar_ext)
    !    if (allocated(kx_shift)) call exb_shear (gnew, phinew, aparnew, bparnew) 
    !-PJK

    allocate( gmodal(lb1:ub1,1:2,lb3:ub3), &
         gnew1(lb1:ub1,1:2,lb3:ub3), &
         gnew2(lb1:ub1,1:2,lb3:ub3), &
         g_dg(lb1:ub1,1:2,lb3:ub3))

    modep = 0
    !if (present(mode)) modep = mode

    !  Apply normalisations to vpar, wdrift

    !  v(ig,:,iglo) is spec(is)%tz * vpar(ig,:,iglo), but with vpar replaced by
    !     its equivalent evaluated *at* grid point ig, not grid-centred at (ig+0.5)
    !  wdrift(ig,isgn,iglo) is already evaluated at ig
    !  N.B. wcoriolis is grid-centred, but is (usually) zero everywhere anyway

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       do ig = -ntgrid,ntgrid
          v(ig,:,iglo) = spec(is)%tz * spec(is)%zstm * tunits(ik)*code_dt &
               * 1.0/delthet(ig) * abs(gradpar(ig))*vpa(ig,:,iglo)

          !  Shift delthet by half a grid point to make it symmetrical...
          !  Worse agreement with implicit case, though
          !if (ig == -ntgrid) then
          !   v(ig,:,iglo) = 0.0
          !else
          !   v(ig,:,iglo) = spec(is)%tz * spec(is)%zstm * tunits(ik)*code_dt &
          !        * 2.0/(delthet(ig)+delthet(ig-1)) * abs(gradpar(ig))*vpa(ig,:,iglo)
          !end if

!+PJK  21/08/2013
#ifdef LOWFLOW
          wd(ig,:,iglo)  = spec(is)%tz * wdfac(ig,:,iglo)
#else
          wd(ig,:,iglo)  = spec(is)%tz * wdrift(ig,:,iglo)
#endif
!-PJK
       end do
       v(ntgrid,:,iglo) = 0.0  !  since delthet(ntgrid) = 0.0...

       !  Zero values outside 'true' range
       v(ntgrid+1:ntgriddg,:,iglo) = 0.0
       wd(ntgrid+1:ntgriddg,:,iglo) = 0.0
    end do

    !  Evaluate fields at old timestep (still need to do this because of
    !  call to g_adjust_exp)

    call getfieldexp (g, phi, apar, bpar, 'g')
    call add_explicit_terms (gexp_1, gexp_2, gexp_3, &  !  used to be 'add_nonlinear_terms'
         phi, apar, bpar, istep, 0.0, z1)

    !  Original call being replaced...
    !  call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)

    call g_adjust_exp (g, apar, fapar)  !  convert g to i

    !  Initialise g_dg, = g but with expanded theta range for DG finite elements

    g_dg = 0.0
    g_dg(-ntgrid:ntgrid,:,:) = g

    !  On first call only, approximate the 'best' choice for the
    !  initial time step

    if (adaptive_dt) then
       dt = code_dt
       adaptive_dt_reset = .false.
       if (first_call == 1) then
          call adaptive_dt0(g_dg,p,code_time,dt0,dt,dtmin,dtmax,epsr,epsa)
          if (dt /= code_dt) then
             adaptive_dt_reset = .true.
             adaptive_dt_new = dt
          end if
          first_call = 0
       end if
    end if

    gmodal = g_dg

    call nodal2modal(gmodal,lb1,ub1,1,2,lb3,ub3,ne2)

    !  Time-advance using Runge-Kutta pair method

    call rk_adaptive_complex(dydt_dggs2,gmodal,gnew1,gnew2,nx,code_time,code_dt)

    call modal2nodal(gnew1,lb1,ub1,1,2,lb3,ub3,ne2)
    call modal2nodal(gnew2,lb1,ub1,1,2,lb3,ub3,ne2)

    !  Adaptive timestep algorithm

    if (adaptive_dt) then

       !  Actual difference between 2nd and 3rd order estimates
       g3dn = maxval(abs(gnew2-gnew1)) ; call max_allreduce(g3dn)

       gomax = maxval(abs(g_dg)) ; call max_allreduce(gomax)
       gnmax = maxval(abs(gnew2)) ; call max_allreduce(gnmax)

       !  (reps *) Maximum allowed difference between 2nd and 3rd order estimates
       g23n = gomax + gnmax + epsar
       rs = reps*g3dn/g23n

       if (rs <= 1.0) then
          !  Step succeeded, but if possible increase timestep for next step

          zhf = min(ffac/rs**zrsord, zfacu)
          if (lfail) then
             zhf = 1.0
             dtmax_damped = min(dtmax_damped,dt)
          end if
          zdt = max(zhf*dt,dtmin)
          zdt = min(zdt,dtmax_damped)  !  Uncomment for 'damped' dt
          dt = zdt
          gnew = gnew1(-ntgrid:ntgrid,:,:)
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
       gnew = gnew1(-ntgrid:ntgrid,:,:)
    end if

    !  Evaluate fields at new timestep

    call getfieldexp (gnew, phinew, aparnew, bparnew, 'i')
    call g_adjust_exp (gnew, aparnew, -fapar)  !  convert from i back to g
    g = gnew

    !  From timeadv (implicit solver)

    !+PJK  21/08/2013
    !Add hyper terms (damping)
    call hyper_diff (gnew, phinew, bparnew)
    !Add collisions
    call vspace_derivatives (gnew, g, g0, phi, apar, bpar, phinew, aparnew, bparnew, modep)
!    if (istep == nstep) call write_mpdist (gnew, '.gtmp', last=.true.)

    !Enforce parity if desired (also performed in invert_rhs, but this is required
    !as collisions etc. may break parity?)
    if (def_parity) call enforce_parity(parity_redist)
      
    !Ensure fixed ik has correct parity parts
    if (fixpar_secondary.gt.0)then
       call enforce_parity(parity_redist_ik,fixpar_secondary)
       gnew = gnew + g_fixpar
    end if
    ! Old stuff...
    !call hyper_diff (gnew, gwork, phinew, bparnew)
    !call kill (gnew, gwork, phinew, bparnew)
    !call solfp1 (gnew, g, gwork, phi, apar, bpar, phinew, aparnew, bparnew, modep)
    !
    !if (def_parity) then
    !   if (even) then
    !      gnew(-ntgrid:-1, 1,:) = gnew( ntgrid: 1:-1,2,:)
    !      gnew( 1: ntgrid, 1,:) = gnew(-1:-ntgrid:-1,2,:)
    !   else
    !      gnew( 1: ntgrid, 1,:) = -gnew(-1:-ntgrid:-1,2,:)
    !      gnew(-ntgrid:-1, 1,:) = -gnew( ntgrid: 1:-1,2,:)
    !   end if
    !end if
    !-PJK

    !+PJK diagnostic tests
    !do iglo = g_lo%llim_proc, g_lo%ulim_proc
    !   il = il_idx(g_lo,iglo) ! il = ittp(ig) = 25 at ig = +-56, +-28, 0
    !   ie = ie_idx(g_lo,iglo)
    !   is = is_idx(g_lo,iglo)
    !   if ((il == 20).and.(ie == 9).and.(is == 1)) then
    !      do ig = -ntgrid,ntgrid
    !         write(65,*) real(ig), abs(g(ig,1,iglo))
    !      end do
    !   end if
    !end do
    !-PJK

    if (allocated(gnew1)) then
       deallocate(gmodal,gnew1,gnew2,g_dg)
    end if

  end subroutine advance_explicit

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine adaptive_dt0(y,p,t,dt0,dt,dtmin,dtmax,epsr,epsa)

    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use mp, only: max_allreduce

    implicit none

    !  Arguments

    complex, dimension(ntgriddgmin:ntgriddg,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: y       !  g_old (actually i_old)

    integer, intent(in) :: p   !  order of the spatial Legendre fit
    real, intent(in) :: t      !  code_time
    real, intent(in) :: dt0    !  code_dt
    real, intent(out) :: dt    !  'best' choice for initial timestep
    real, intent(in) :: dtmin  !  minimum allowed timestep
    real, intent(out) :: dtmax !  estimate of maximum allowed timestep
    real, intent(in) :: epsr   !  relative tolerance level
    real, intent(in) :: epsa   !  absolute tolerance level

    !  Local variables

    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    real :: gdotn,gn,epst
    complex, allocatable, dimension(:,:,:) :: ymodal
    complex, allocatable, dimension(:,:,:) :: dy

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lb1 = ntgriddgmin ; ub1 = ntgriddg
    lb2 = 1 ; ub2 = 2
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc

    allocate(ymodal(lb1:ub1,lb2:ub2,lb3:ub3))
    allocate(dy(lb1:ub1,lb2:ub2,lb3:ub3))

    !  Calculate dg/dt using the Discontinuous Galerkin scheme

    ymodal = y
    call nodal2modal(ymodal,lb1,ub1,1,2,lb3,ub3,ne2)
    call dydt_dggs2(ymodal,t,dt0,dy)
    call modal2nodal(dy,lb1,ub1,1,2,lb3,ub3,ne2)

    gdotn = maxval(abs(dy)) ; call max_allreduce(gdotn)
    gn = maxval(abs(y)) ; call max_allreduce(gn)

    epst = epsr*gn + epsa

    dtmax = (epst/gdotn)**(1.0/p)

    dt = min(dt0,dtmax)
    dt = max(dt,dtmin)

    deallocate(ymodal,dy)

  end subroutine adaptive_dt0

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dydt_dggs2(gmodal,t,dt,dgdt_modal)

    use dist_fn, only: getfieldexp, get_source_term_exp
    use dist_fn_arrays, only: ittp
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    use le_grids, only: nlambda, ng2, lmax, forbid
    use run_parameters, only: tunits
    use theta_grid, only: ntgrid

    implicit none

    !  Arguments

    real, intent(in) :: t, dt
    complex, dimension(ntgriddgmin:ntgriddg,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(in) :: gmodal
    complex, dimension(ntgriddgmin:ntgriddg,1:2,g_lo%llim_proc:g_lo%ulim_proc), &
         intent(out) :: dgdt_modal

    !  Local variables

    integer :: ip, ip2
    integer :: i, j, kk, element, ig, igg, ik, is, il, ie, isgn, iglo
    integer :: lb1, ub1, lb3, ub3, lbp, ubp, lb1a

    complex, parameter :: zi = (0.0,1.0)
    logical :: trapped

    complex, allocatable, dimension(:,:,:) :: g, g_dg
    complex, allocatable, dimension(:,:,:) :: dgdt, dgdt_tmp
    complex, allocatable, dimension(:,:,:) :: dgdt_mod
    complex, allocatable, dimension(:) :: iwg, iwg_modal
    complex, allocatable, dimension(:) :: gplusf, gplusf_modal
    complex, allocatable, dimension(:) :: mgplusf, mgplusf_modal
    complex, allocatable, dimension(:) :: vterm, vterm_modal
    complex, allocatable, dimension(:) :: src, smodal
    complex, allocatable, dimension(:,:,:) :: src3, fluxfn

    real, dimension(3,6) :: mat_p3p, mat_p3m

    complex, dimension(2*p) :: ftemp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Advection term matrix, p=3, positive v

    mat_p3p = real(reshape( &
         source = (/ 1,-3,5, 1,-3,5, 1,-3,5, -1,3,-5, -1,-3,5, -1,-3,-5 /), &
         shape = (/ 3,6 /) ))

    !  Advection term matrix, p=3, negative v

    mat_p3m = real(reshape( &
         source = (/ 1,3,5, -1,3,5, 1,-3,5, -1,-3,-5, 1,3,5, -1,-3,-5 /), &
         shape = (/ 3,6 /) ))

    lb1 = -ntgrid ; ub1 = ntgriddg
    lb1a = ntgriddgmin
    lb3 = g_lo%llim_proc ; ub3 = g_lo%ulim_proc

    allocate( &
         g(-ntgrid:ntgrid,2,lb3:ub3), & ! original theta range
         g_dg(lb1a:ub1,2,lb3:ub3), &
         dgdt(lb1:ub1,2,lb3:ub3), &
         dgdt_tmp(lb1:ub1,2,lb3:ub3), &
         dgdt_mod(lb1:ub1,2,lb3:ub3), &
         iwg(lb1:ub1), &
         iwg_modal(lb1:ub1), &
         gplusf(lb1a:ub1), &
         gplusf_modal(lb1a:ub1), &
         mgplusf(lb1:ub1), &
         mgplusf_modal(lb1:ub1), &
         vterm(lb1:ub1), &
         vterm_modal(lb1:ub1), &
         fluxfn(lb1a:ub1,2,lb3:ub3), &
         src(lb1:ub1), &
         src3(lb1:ub1,2,lb3:ub3), &
         smodal(lb1:ub1) )

    g_dg = gmodal

    !  g halo elements need to be populated, as these are used in the
    !  advection term, which takes information from upwind finite elements

    if (flux_tube) then
       call pass_fluxtube_halos(g_dg,ntgrid,p,lb3,ub3)

       !  Force periodic boundary conditions for the ky=0 zonal mode
       !  (for which all theta points are on the same MPI process)

       do iglo = lb3, ub3
          ik = ik_idx(g_lo,iglo)
          if (ik == 1) then
             g_dg(lb1-p:lb1-1,1,iglo) = g_dg(ntgrid-p:ntgrid-1,1,iglo)
             g_dg(ntgrid:ntgrid+p-1,2,iglo) = g_dg(lb1:lb1+p-1,2,iglo)
          end if
       end do
    end if

    !  Obtain g in nodal form

    call modal2nodal(g_dg,lb1a,ub1,1,2,lb3,ub3,ne2)
    g = g_dg(-ntgrid:ntgrid,:,:)

    !  Evaluate field equations

    call getfieldexp(g,phitmp,apartmp,bpartmp,'i')

    !  Calculate source term and flux function (nodal)

    src3 = 0.0 ; fluxfn = 0.0
    call get_source_term_exp(phitmp,apartmp,bpartmp,istep_dg, &
         fluxfn(lb1:ntgrid,:,:),src3(lb1:ntgrid,:,:))

    !  Flux function halos need to be populated, as these are used in the
    !  advection term, which takes information from upwind finite elements

    if (flux_tube) call pass_fluxtube_halos(fluxfn,ntgrid,p,lb3,ub3)

    do iglo = lb3, ub3
       ik = ik_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       !  Force periodic boundary conditions for the ky=0 zonal mode
       !  (for which all theta points are on the same MPI process)

       if ((flux_tube).and.(ik == 1)) then
          fluxfn(lb1-p:lb1-1,1,iglo) = fluxfn(ntgrid-p:ntgrid-1,1,iglo)
          fluxfn(ntgrid:ntgrid+p-1,2,iglo) = fluxfn(lb1:lb1+p-1,2,iglo)
       end if

       !  Loop over sign of the parallel velocity
       !  isgn=1, v >= 0 ; isgn=2, v < 0

       do isgn = 1,2

          !  Evaluate -i*wd*g term
          !  Have to divide wd by dt here

          iwg = 0.0
          do ig = -ntgrid,ntgrid
             if (forbid(ig,il)) then
                iwg(ig) = 0.0
             else
                iwg(ig) = -zi * wd(ig,isgn,iglo) * g(ig,isgn,iglo) / (dt*tunits(ik))
             end if
          end do

          iwg_modal = iwg
          call nodal2modal(iwg_modal,lb1,ub1,ne)

          !  Add g to F

          gplusf(:) = g_dg(:,isgn,iglo) + fluxfn(:,isgn,iglo)
          gplusf_modal = gplusf

          if (.not.flux_tube) then
             call nodal2modal(gplusf_modal,lb1,ub1,ne2)
          else
             call nodal2modal(gplusf_modal,lb1a,ub1,ne2)
          end if

          !  Perform matrix multiplication to calculate d(g+F)/dz
          !  For speed, ought to store and use transpose...

          mgplusf_modal(:) = 0.0

          if (isgn == 1) then  !  Positive parallel velocity

             if (.not.flux_tube) then

                do element = 1,ne
                   !  Fill ftemp with correct elements from (g+F)modal
                   if (element /= 1) then
                      kk = lb1 + p*(element-2)  !  take values from element to left
                      ftemp(1:p) = gplusf_modal(kk:kk+p-1)
                   else
                      ftemp(1:p) = 0.0  !  'ghost' element to left contains zeroes
                   end if
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1)

                   do ip = 1,p
                      kk = lb1 + p*(element-1) + ip-1
                      do ip2 = 1,2*p
                         mgplusf_modal(kk) = mgplusf_modal(kk) &
                              + mat_p3p(ip,ip2)*ftemp(ip2)
                      end do
                   end do
                end do

             else  !  flux-tube scenario

                do element = 1,ne
                   !  Fill ftemp with correct elements from (g+F)modal
                   !  Element '0' contains halo values
                   kk = lb1 + p*(element-2)  !  take values from element to left
                   ftemp(1:p) = gplusf_modal(kk:kk+p-1)
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1)

                   do ip = 1,p
                      kk = lb1 + p*(element-1) + ip-1
                      do ip2 = 1,2*p
                         mgplusf_modal(kk) = mgplusf_modal(kk) &
                              + mat_p3p(ip,ip2)*ftemp(ip2)
                      end do
                   end do
                end do

             end if

          else  !  Negative parallel velocity

             if (.not.flux_tube) then

                do element = 1,ne
                   !  Fill ftemp with correct elements from (g+F)modal
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(1:p) = gplusf_modal(kk:kk+p-1)
                   if (element /= ne) then
                      kk = lb1 + p*element  !  take values from element to right
                      ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1)
                   else
                      ftemp(p+1:2*p) = 0.0  !  'ghost' element to right contains zeros
                   end if

                   do ip = 1,p
                      kk = lb1 + p*(element-1) + ip-1
                      do ip2 = 1,2*p
                         mgplusf_modal(kk) = mgplusf_modal(kk) &
                              + mat_p3m(ip,ip2)*ftemp(ip2)
                      end do
                   end do
                end do

             else  !  flux-tube scenario

                do element = 1,ne
                   !  Fill ftemp with correct elements from (g+F)modal
                   !  Element 'ne+1' contains halo values
                   kk = lb1 + p*(element-1)  !  take values from this element
                   ftemp(1:p) = gplusf_modal(kk:kk+p-1)
                   kk = lb1 + p*element  !  take values from element to right
                   ftemp(p+1:2*p) = gplusf_modal(kk:kk+p-1)

                   do ip = 1,p
                      kk = lb1 + p*(element-1) + ip-1
                      do ip2 = 1,2*p
                         mgplusf_modal(kk) = mgplusf_modal(kk) &
                              + mat_p3m(ip,ip2)*ftemp(ip2)
                      end do
                   end do
                end do

             end if

          end if

          !  Multiply this term by v/dz, which must be done in nodal form
          !  Array v = (T/q)*vp, where vp = dt/h*(parallel velocity),
          !  i.e. vp is non-dimensional.
          !  h is the finite difference mesh separation, but dz must be
          !  element size p*h, not h.
          !  We must also remove the numerator dt, so we need to divide v by p*dt

          mgplusf = mgplusf_modal
          call modal2nodal(mgplusf,lb1,ub1,ne)

          do ig = lb1,ub1
             if (ig > ntgrid) then
                vterm(ig) = 0.0
             else if (forbid(ig,il)) then
                vterm(ig) = 0.0
             else
                vterm(ig) = 1.0/(p*dt*tunits(ik)) * v(ig,isgn,iglo) &
                     * mgplusf(ig)
             end if
          end do

          vterm_modal = vterm
          call nodal2modal(vterm_modal,lb1,ub1,ne)

          !  Convert source term to modal form.
          !  src3 is actually (2.dt.S), therefore need to divide by 2.dt

          src = 0.0
          src(:) = src3(:,isgn,iglo)/(2.0*dt*tunits(ik))
          smodal = src
          call nodal2modal(smodal,lb1,ub1,ne)

          !  Sum all the modal terms to get dg/dt (modal)

          do i = lb1,ub1 ! modes
             dgdt_mod(i,isgn,iglo) = &
                  iwg_modal(i) &
                  + vterm_modal(i) &
                  + smodal(i)
          end do

       end do  !  isgn loop

    end do  !  iglo loop

    dgdt_tmp = dgdt_mod(:,:,:)
    call modal2nodal(dgdt_tmp,lb1,ub1,1,2,lb3,ub3,ne)

    !  Apply relevant boundary conditions to passing/trapped particles

    dgdt = 0.0

    do iglo = lb3,ub3
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       trapped = (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax)

       if (trapped) then
          lbp = -9999 ; ubp = -9999
          do ig = lb1,ntgrid-1

             !+PJK prototype stuff for totally trapped particles
             if (il == ittp(ig)) then
                if (forbid(ig,il)) then
                   dgdt(ig,:,iglo) = 0.0
                else
                   dgdt(ig,1,iglo) = dgdt_tmp(ig,1,iglo)
                   dgdt(ig,2,iglo) = dgdt(ig,1,iglo)
                end if
                cycle
             end if
             !-PJK

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

                do igg = lbp+1,ubp-1
                   dgdt(igg,:,iglo) = dgdt_tmp(igg,:,iglo)
                end do

                !  Boundary condition at bounce points:
                !  dgdt(v > 0) = dgdt(v < 0)

                dgdt(lbp,1,iglo) = dgdt_tmp(lbp,2,iglo)
                dgdt(lbp,2,iglo) = dgdt_tmp(lbp,2,iglo)

                dgdt(ubp,1,iglo) = dgdt_tmp(ubp,1,iglo)
                dgdt(ubp,2,iglo) = dgdt_tmp(ubp,1,iglo)

                lbp = -9999 ; ubp = -9999
             end if
          end do

       else  !  nothing trapped at this il
          dgdt(:,:,iglo) = dgdt_tmp(:,:,iglo)
       end if

    end do ! iglo

    !  Boundary conditions
!+PJKFT Check this is correct for FT; actually only applies for left-most
!+PJKFT or right-most procs
    dgdt(lb1,1,:) = 0.0
    dgdt(ntgrid:ntgriddg,2,:) = 0.0
!-PJKFT

    !  Finally, convert dg/dt back to modal space

    if (.not.flux_tube) then
       dgdt_modal = dgdt
    else
       dgdt_modal = 0.0
       dgdt_modal(lb1:ub1,:,:) = dgdt
    end if
    call nodal2modal(dgdt_modal,lb1a,ub1,1,2,lb3,ub3,ne2)

  end subroutine dydt_dggs2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rk_adaptive_complex(dydt,y,ynew1,ynew2,nx,t,dt)

    !  Runge-Kutta solver

    use rk_schemes, rk => rkf3

    implicit none

    !  Arguments

    !  Two output values for ynew: one for each of the two orders of the
    !  chosen adaptive RK scheme

    external :: dydt  !  subroutine used to calculate dy/dt at a given x, t

    integer, intent(in) :: nx
    complex, dimension(nx), intent(in) :: y       !  y_old(x) at time t
    complex, dimension(nx), intent(out) :: ynew1  !
    complex, dimension(nx), intent(out) :: ynew2
    real, intent(in) :: t
    real, intent(in) :: dt

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

       call dydt(ytest,ttest,dt,ydot)
       k(:,i) = ydot
    end do

    ynew1 = y
    ynew2 = y

    do i = 1,rk%stages
       ynew1 = ynew1 + dt*rk%bi(1,i)*k(:,i)
       ynew2 = ynew2 + dt*rk%bi(2,i)*k(:,i)
    end do

  end subroutine rk_adaptive_complex

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pass_fluxtube_halos(array, ntgrid, p, lb3, ub3)

    !  Pass the contents of halo finite elements between processors,
    !  during flux-tube runs
    !
    !  Code adapted from David Dickinson's ensure_single_val_g routine
    !
    !  The input/output argument array must have dimensions
    !  array(-ntgrid-p:ntgrid+p-1, 1:2, g_lo%llim_proc:g_lo%ulim_proc)
    !
    !  P J Knight
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use redistribute, only: fill, delete_redist
    use dist_fn, only: init_pass_dg, pass_right, pass_left

    implicit none

    !  Arguments

    integer, intent(in) :: ntgrid, p, lb3, ub3
    complex, dimension (-ntgrid-p:ntgrid+p-1,1:2,lb3:ub3), &
         intent(inout) :: array

    logical :: initd = .false.

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Initialise left and right communicators
    !  (need only do this at the start of a run)

    if (.not.initd) then
       call init_pass_dg(pass_right,'r')
       call init_pass_dg(pass_left,'l')
       initd = .true.
    end if

    !  Pass to the right and left

    call fill(pass_right,array,array)
    call fill(pass_left,array,array)

    !  Deallocate communicator
    !  (currently not done; should move these into a destructor routine)
    !call delete_redist(pass_right)
    !call delete_redist(pass_left)

  end subroutine pass_fluxtube_halos

end module fields_explicit
