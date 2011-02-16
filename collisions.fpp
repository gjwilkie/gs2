# include "define.inc"

module collisions
  
  use redistribute, only: redist_type

  implicit none

  public :: init_collisions, finish_collisions
  public :: solfp1
  public :: reset_init
  public :: dtot, fdf, fdb, lorentz_map, vnmult, vnfac
  public :: ncheck, vnslow, vary_vnew
  public :: etol, ewindow, etola, ewindowa
  public :: init_lorentz, init_ediffuse
  public :: init_lorentz_conserve, init_diffuse_conserve
  public :: init_lorentz_error, collision_model_switch
  public :: g_adjust ! CMR this routine would be better placed in dist_fn.f90

  private

  ! knobs
!> only needed for krook
  real :: vncoef, absom
  integer :: ivnew
  logical :: conserve_number, conserve_momentum
!<
  logical :: const_v, conserve_moments
  logical :: conservative, resistivity
  integer :: collision_model_switch
  integer :: lorentz_switch, ediff_switch
  logical :: adjust
  logical :: heating
  logical :: hyper_colls
  logical :: ei_coll_only
  logical :: test

  integer, parameter :: collision_model_lorentz = 1      ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_krook = 2
  integer, public, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_krook_test = 4
  integer, parameter :: collision_model_lorentz_test = 5 ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_full = 6
  integer, parameter :: collision_model_ediffuse = 7

  integer, parameter :: lorentz_scheme_default = 1
  integer, parameter :: lorentz_scheme_old = 2

  integer, parameter :: ediff_scheme_default = 1
  integer, parameter :: ediff_scheme_old = 2

  real, dimension (2), save :: vnmult = 0.0
  integer :: ncheck
  logical :: vary_vnew
  real :: vnfac, vnslow
  real :: etol, ewindow, etola, ewindowa

  real, dimension (:,:,:), allocatable :: dtot
  ! (-ntgrid:ntgrid,nlambda,max(ng2,nlambda-ng2)) lagrange coefficients for derivative error estimate

  real, dimension (:,:), allocatable :: fdf, fdb
  ! (-ntgrid,ntgrid,nlambda) finite difference coefficients for derivative error estimate

  real, dimension (:,:,:), allocatable :: vnew, vnew_s, vnew_D, vnew_E, delvnew
  ! (naky,negrid,nspec) replicated

  real, dimension (:,:,:), allocatable :: vpdiff
  ! (-ntgrid:ntgrid,2,nlambda) replicated

  complex, dimension (:,:,:), allocatable :: gle, glec
  ! (2*nlambda+1,negrid+1,le_lo%llim_proc:le_lo%ulim_alloc)

  ! only for hyper-diffusive collisions
  real, dimension (:,:,:,:), allocatable :: vnewh
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec) replicated

  ! only for krook
  real, dimension (:), allocatable :: vnewfe
  ! (-*- g_layout -*-)

  ! only for krook
  real, dimension (:,:), allocatable :: aintnorm
  ! (-ntgrid:ntgrid, -*- geint_layout -*-)

  ! for conservation of momentum and energy
!  complex, dimension (:,:,:,:), allocatable :: parfac, perpfac, efac
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec)

  ! only for momentum conservation due to Lorentz operator (8.06)
  complex, dimension(:,:,:), allocatable :: s0, w0, z0
# ifdef USE_LE_LAYOUT
  real, dimension (:,:,:), allocatable :: s0le, w0le, z0le
  real, dimension (:,:,:), allocatable :: aj0le, aj1le
# endif

  ! needed for momentum and energy conservation due to energy diffusion (3.08)
  complex, dimension(:,:,:), allocatable :: bs0, bw0, bz0
# ifdef USE_LE_LAYOUT
  complex, dimension (:,:,:), allocatable :: bs0le, bw0le, bz0le
# endif

  ! only for original parallel mom conservation (not used nowadays)
  real, dimension (:,:,:), allocatable :: sq
  ! (-ntgrid:ntgrid,nlambda,2) replicated

  real :: cfac
# ifdef USE_LE_LAYOUT
  ! only for lorentz
  real, dimension (:,:,:), allocatable :: c1le, betaale, qle, d1le, h1le
  ! only for energy diffusion
  real, dimension (:,:,:), allocatable :: ec1le, ebetaale, eqle
# else
  ! only for lorentz
  real, dimension (:,:), allocatable :: c1, betaa, ql, d1, h1
  ! only for energy diffusion
  real, dimension (:,:), allocatable :: ec1, ebetaa, eql
# endif


  ! momentum conservation
  complex, dimension (:,:), allocatable :: g3int
  ! (-ntgrid:ntgrid, -*- gint_layout -*-)

  type (redist_type), save :: lorentz_map
  type (redist_type), save :: ediffuse_map
  type (redist_type), save :: g2le

  logical :: hypermult
  logical :: initialized = .false.
  logical :: ediffinit = .false.
  logical :: lzinit = .false., leinit = .false.
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.

contains

  subroutine init_collisions
    use species, only: init_species, nspec, spec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0
    use le_grids, only: init_le_grids, nlambda, negrid 
    use run_parameters, only: init_run_parameters
    use gs2_layouts, only: init_dist_fn_layouts, init_gs2_layouts

    implicit none

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts
    call init_species

    hyper_colls = .false.
    if (any(spec%nu_h > epsilon(0.0))) hyper_colls = .true.

    call init_theta_grid
    call init_kt_grids
    call init_le_grids (accelerated_x, accelerated_v)
    call init_run_parameters
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
    call read_parameters
    call init_map  ! ported from AGK July 24, 2009 -- MAB
    call init_arrays

  end subroutine init_collisions

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (8), parameter :: modelopts = &
         (/ text_option('default', collision_model_full), &
            text_option('lorentz', collision_model_lorentz), &
            text_option('ediffuse', collision_model_ediffuse), &
            text_option('krook', collision_model_krook), &
            text_option('krook-test', collision_model_krook_test), &
            text_option('lorentz-test', collision_model_lorentz_test), &
            text_option('none', collision_model_none), &
            text_option('collisionless', collision_model_none) /)
    type (text_option), dimension (2), parameter :: schemeopts = &
         (/ text_option('default', lorentz_scheme_default), &
            text_option('old', lorentz_scheme_old) /)
    type (text_option), dimension (2), parameter :: eschemeopts = &
         (/ text_option('default', ediff_scheme_default), &
            text_option('old', ediff_scheme_old) /)
    character(20) :: collision_model, lorentz_scheme, ediff_scheme
    namelist /collisions_knobs/ collision_model, conserve_moments, &
         heating, adjust, const_v, cfac, hypermult, ei_coll_only, &
         lorentz_scheme, ediff_scheme, resistivity, conservative, test, &
! following only needed for adaptive collisionality
         vnfac, etol, ewindow, ncheck, vnslow, vary_vnew, etola, ewindowa, &
! following only needed for krook
         vncoef, absom, ivnew, conserve_number, conserve_momentum
    integer :: ierr, in_file
    logical :: exist

    if (proc0) then
       collision_model = 'default'
       conserve_moments = .true.  ! DEFAULT CHANGED TO REFLECT IMPROVED MOMENTUM AND ENERGY CONSERVATION 7/08
       hypermult = .false.
       cfac = 1.   ! DEFAULT CHANGED TO INCLUDE CLASSICAL DIFFUSION: APRIL 18, 2006
!> following only used with adaptive collisionality
       vary_vnew = .false.
       vnfac = 1.1
       vnslow = 0.9
       etol = 2.e-2
       ewindow = 1.e-2
       etola = 2.e-2
       ewindowa = 1.e-2
       ncheck = 100
!<
       adjust = .true.
       lorentz_scheme = 'default'
       ediff_scheme = 'default'
!> following only needed for krook
       vncoef = 0.6
       absom = 0.5
       ivnew = 0
       conserve_number = .true.
       conserve_momentum = .true.  ! DEFAULT CHANGED TO REFLECT IMPROVED MOMENTUM CONSERVATION, 8/06
!>
       conservative = .true.
       resistivity = .true.
       const_v = .false.
       heating = .false.
       test = .false.
       ei_coll_only = .false.
       in_file = input_unit_exist ("collisions_knobs", exist)
!       if (exist) read (unit=input_unit("collisions_knobs"), nml=collisions_knobs)
       if (exist) read (unit=in_file, nml=collisions_knobs)

       ierr = error_unit()
       call get_option_value &
            (collision_model, modelopts, collision_model_switch, &
            ierr, "collision_model in collisions_knobs")

       call get_option_value &
            (lorentz_scheme, schemeopts, lorentz_switch, &
            ierr, "lorentz_scheme in collisions_knobs")

       call get_option_value &
            (ediff_scheme, eschemeopts, ediff_switch, &
            ierr, "ediff_scheme in collisions_knobs")
    end if

    call broadcast (hypermult)
    call broadcast (cfac)
!> only need for adaptive collisionality
    call broadcast (vnfac)
    call broadcast (vnslow)
    call broadcast (vary_vnew)
    call broadcast (etol)
    call broadcast (ewindow)
    call broadcast (etola)
    call broadcast (ewindowa)
    call broadcast (ncheck)
!<
!> only needed for krook
    call broadcast (vncoef)
    call broadcast (absom)
    call broadcast (ivnew)
    call broadcast (conserve_number)
    call broadcast (conserve_momentum)
!<
    call broadcast (conservative)
    call broadcast (conserve_moments)
    call broadcast (resistivity)
    call broadcast (const_v)
    call broadcast (collision_model_switch)
    call broadcast (lorentz_switch)
    call broadcast (ediff_switch)
    call broadcast (heating)
    call broadcast (test)
    call broadcast (adjust)
    call broadcast (ei_coll_only)
  end subroutine read_parameters

  subroutine init_map

    use mp, only: finish_mp, proc0
    use redistribute, only: report_map_property

# ifdef USE_LE_LAYOUT
    call init_g2le_redistribute ! new map
    if (test) call check_g2le
# else
    select case (collision_model_switch)
    case (collision_model_full)
       ! init_lorentz_layout is called in redistribute
       call init_lorentz_redistribute
       if (test) then
          if (proc0) print *, '=== Lorentz map property ==='
          call report_map_property (lorentz_map)
       end if
       ! init_ediffuse_layout is called in redistribute
       call init_ediffuse_redistribute
       if (test) then
          if (proc0) print *, '=== Ediffuse map property ==='
          call report_map_property (ediffuse_map)
       end if
    case (collision_model_lorentz,collision_model_lorentz_test)
       ! init_lorentz_layout is called in redistribute
       call init_lorentz_redistribute
    case (collision_model_ediffuse)
       ! init_ediffuse_layout is called in redistribute
       call init_ediffuse_redistribute
    end select
# endif
    
    if (test) then
       if (proc0) print *, 'init_map done'
!!$    call finish_mp
!!$    stop
    end if

  end subroutine init_map

  subroutine init_arrays
    use species, only: nspec
    use le_grids, only: negrid
    use gs2_layouts, only: g_lo
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: c_rate

    implicit none

    real, dimension (negrid,nspec) :: hee

    if (.not. allocated(c_rate)) then
       allocate (c_rate(-ntgrid:ntgrid, ntheta0, naky, nspec, 3))
       c_rate = 0.
    end if

    if (collision_model_switch == collision_model_none) return

    call init_vnew (hee)
    if (all(abs(vnew(:,1,:)) <= 2.0*epsilon(0.0))) then
       collision_model_switch = collision_model_none
       return
    end if
    if (conserve_momentum .and. &
         collision_model_switch == collision_model_krook) call init_g3int

    select case (collision_model_switch)
    case (collision_model_full)
       call init_lorentz
       call init_ediffuse
       if (conserve_moments) then
          call init_lorentz_conserve
          call init_diffuse_conserve
       end if
    case (collision_model_lorentz,collision_model_lorentz_test)
       call init_lorentz
       if (conserve_moments) call init_lorentz_conserve
    case (collision_model_ediffuse)
       call init_ediffuse
       if (conserve_moments) call init_diffuse_conserve
    case (collision_model_krook,collision_model_krook_test)
       call init_krook (hee)
    end select

  end subroutine init_arrays

  subroutine init_lorentz_conserve

! Precompute three quantities needed for momentum and energy conservation:
! z0, w0, s0
    
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec, electron_species
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: e, al, integrate_moment, negrid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, kperp2, vpa
    use run_parameters, only: tunits
    use constants, only: pi
# ifdef USE_LE_LAYOUT
    use le_grids, only: nlambda, ng2
    use gs2_layouts, only: le_lo
    use redistribute, only: gather
# endif

    implicit none
    
    logical :: init_flag = .true.
    complex, dimension (1,1,1) :: dum1 = 0., dum2 = 0.
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: duinv, dtmp
    real, dimension (:,:,:,:), allocatable :: vns
    integer :: ie, il, ik, is, isgn, iglo, all, it
# ifdef USE_LE_LAYOUT
    integer :: nxi
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
# endif
    
! TO DO: 
! tunits not included anywhere yet

    if (.not. allocated(z0)) then
       allocate (z0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (w0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (s0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    end if
    
    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (dtmp(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (vns(naky,negrid,nspec,3))
    
    vns(:,:,:,1) = vnmult(1)*vnew_D
    vns(:,:,:,2) = vnmult(1)*vnew_s
    vns(:,:,:,3) = 0.0

    if (resistivity .and. nspec > 1) then
       do is = 1, nspec
          if (spec(is)%type /= electron_species) cycle
          do ik = 1, naky
             vns(ik,:,is,3) = vnmult(1)*spec(is)%vnewk*tunits(ik)/e(:,is)**1.5
          end do
       end do
    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get z0 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! u0 = -2 nu_D^{ei} vpa J0 dt f0
          if (conservative) then
             z0(:,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpdiff(:,isgn,il) &
                  * sqrt(e(ie,is))*aj0(:,iglo)
          else
             z0(:,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpa(:,isgn,iglo)*aj0(:,iglo)
          end if
       end do
    end do
    
    call solfp_lorentz (z0,dum1,dum2,init=init_flag)   ! z0 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*z0(:,isgn,iglo)
       end do
    end do
    
    call integrate_moment (gtmp, dtmp, all) ! v0z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine z0 = z0 / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          z0(:,isgn,iglo) = z0(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! du == int (E nu_s f_0);  du = du(z, kx, ky, s)
    ! duinv = 1/du
    if (conservative) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             gtmp(:,isgn,iglo)  = vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
                  * vpdiff(:,isgn,il)*sqrt(e(ie,is))
!             gtmp(:,isgn,iglo)  = vns(ik,ie,is,2)*e(ie,is)
          end do
       end do

       all = 1
       call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
    else
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             gtmp(:,isgn,iglo)  = vns(ik,ie,is,1)*vpa(:,isgn,iglo)**2
          end do
       end do

       all = 1
       call integrate_moment (gtmp, duinv, all)  ! not 1/du yet

!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(1)*spec(is)%vnewk*sqrt(2./pi)
!       end do
    end if

    where (cabs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                        ! duinv=0 iff vnew=0 so ok to keep duinv=0.
       duinv = 1./duinv  ! now it is 1/du
    end where

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get s0 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! u1 = -3 nu_s vpa dt J0 f_0 / du
          if (conservative) then
             s0(:,isgn,iglo) = -vns(ik,ie,is,1)*vpdiff(:,isgn,il)*sqrt(e(ie,is)) &
                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
          else
!             s0(:,isgn,iglo) = -3.0*vns(ik,ie,is,2)*vpa(:,isgn,iglo) &
             s0(:,isgn,iglo) = -vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
          end if
       end do
    end do

    call solfp_lorentz (s0,dum1,dum2,init=init_flag)    ! s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0s0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*s0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v0s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 - v0s0 * z0 / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          s0(:,isgn,iglo) = s0(:,isgn,iglo) - dtmp(:,it,ik,is) &
               * z0(:,isgn,iglo)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1s0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nu_D vpa J0
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*sqrt(e(ie,is))*vpdiff(:,isgn,il) &
                  * aj0(:,iglo)*s0(:,isgn,iglo)
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * s0(:,isgn,iglo)
          end if
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v1s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 / (1 + v0s0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          s0(:,isgn,iglo) = s0(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get w0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! u2 = -3 dt J1 vperp vus a f0 / du
          if (conservative) then
             w0(:,isgn,iglo) = -vns(ik,ie,is,1)*e(ie,is)*al(il)*aj1(:,iglo) &
                  * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                  / bmag
          else
!             w0(:,isgn,iglo) = -3.*vns(ik,ie,is,2)*e(ie,is)*al(il)*aj1(:,iglo) &
             w0(:,isgn,iglo) = -vns(ik,ie,is,1)*e(ie,is)*al(il)*aj1(:,iglo) &
                  * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                  / bmag
          end if
       end do
    end do

    call solfp_lorentz (w0,dum1,dum2,init=init_flag)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0w0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*w0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v0w0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w0 - v0w0 * z0 / (1 + v0z0) (this is w1 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          w0(:,isgn,iglo) = w0(:,isgn,iglo) - z0(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v1w1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nud vpa J0 f0
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*sqrt(e(ie,is))*vpdiff(:,isgn,il) &
                  * aj0(:,iglo)*w0(:,isgn,iglo)
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * w0(:,isgn,iglo)
          end if
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v1w1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w1 - v1w1 * s1 / (1 + v1s1) (this is w2 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          w0(:,isgn,iglo) = w0(:,isgn,iglo) - s0(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v2w2

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = nud vperp J1 f0 
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*e(ie,is)*al(il)*aj1(:,iglo) &
               * w0(:,isgn,iglo)
       end do
    end do
    
    call integrate_moment (gtmp, dtmp, all)   ! v2w2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w2 / (1 + v2w2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          w0(:,isgn,iglo) = w0(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

    deallocate (gtmp, duinv, dtmp, vns)

# ifdef USE_LE_LAYOUT

    nxi = max(2*nlambda-1,2*ng2)
    allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))

    ! first set s0le, w0le & z0le
    z_big = 0.0
    z_big = cmplx(real(s0), real(w0))

    call gather (g2le, z_big, ctmp)

    if (.not. allocated(z0le)) then
       allocate (s0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (w0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (z0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end if

    s0le = real(ctmp)
    w0le = aimag(ctmp)

    ! set z0le
    call gather (g2le, z0, ctmp)

    z0le = real(ctmp)

    ! next set aj0le & aj1l
    z_big(:,1,:) = cmplx(aj0,aj1)
    z_big(:,2,:) = z_big(:,1,:)

    call gather (g2le, z_big, ctmp)

    if (.not. allocated(aj0le)) then
       allocate (aj0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (aj1le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end if

    aj0le = real(ctmp)
    aj1le = aimag(ctmp)

    deallocate (ctmp, z_big)

# endif
    
  end subroutine init_lorentz_conserve

  subroutine init_diffuse_conserve

! Precompute three quantities needed for momentum and energy conservation:
! bz0, bw0, bs0
    
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: e, al, integrate_moment, negrid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, kperp2, vpa
    use run_parameters, only: tunits
    use constants, only: pi
# ifdef USE_LE_LAYOUT
    use le_grids, only: nlambda, ng2
    use gs2_layouts, only: le_lo
    use redistribute, only: gather
# endif

    implicit none

    logical :: init = .true.
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: duinv, dtmp
    real, dimension (:,:,:,:), allocatable :: vns
    integer :: ie, il, ik, is, isgn, iglo, all, it
# ifdef USE_LE_LAYOUT
    integer :: nxi
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
# endif

! TO DO: 
! tunits not included anywhere yet

!    if (first) then
    if (.not. allocated(bz0)) then
       allocate (bz0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (bw0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (bs0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
!       first = .false.
    end if

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (dtmp(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (vns(naky,negrid,nspec,2))
    
    vns(:,:,:,1) = vnmult(2)*delvnew
    vns(:,:,:,2) = vnmult(2)*vnew_s

    ! first obtain 1/du
!    if (conservative) then
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          gtmp(:,isgn,iglo) = e(ie,is)*vnmult(2)*vnew_E(ik,ie,is)
       end do
    end do
       
    all = 1
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
!    else
!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(2)*spec(is)%vnewk*sqrt(2./pi)
!       end do
!    end if

    where (cabs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                        ! duinv=0 iff vnew=0 so ok to keep duinv=0.
       duinv = 1./duinv  ! now it is 1/du
    end where

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get z0 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! u0 = -nu_E E dt J0 f_0 / du
           bz0(:,isgn,iglo) = -code_dt*vnmult(2)*vnew_E(ik,ie,is) &
               * aj0(:,iglo)*duinv(:,it,ik,is)
       end do
    end do
    
    call solfp_ediffuse (bz0,init)   ! s0 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0 f_0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * bz0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all) ! v0z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine z0 = z0 / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          bz0(:,isgn,iglo) = bz0(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! redefine dq = du (for momentum-conserving terms)
    ! du == int (E nu_s f_0);  du = du(z, kx, ky, s)
    ! duinv = 1/du
!    if (conservative) then
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          gtmp(:,isgn,iglo)  = vns(ik,ie,is,1)*vpa(:,isgn,iglo)**2
!             gtmp(:,isgn,iglo)  = vns(ik,ie,is,2)*e(ie,is)

       end do
    end do
       
    all = 1
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
!    else
!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(2)*spec(is)%vnewk*sqrt(2./pi)
!       end do
!    end if

    where (cabs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                        ! duinv=0 iff vnew=0 so ok to keep duinv=0.
       duinv = 1./duinv  ! now it is 1/du
    end where

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get s0 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! u1 = -3 nu_s vpa dt J0 f_0 / du
!          if (conservative) then
          bs0(:,isgn,iglo) = -vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
               * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
!          else
!             bs0(:,isgn,iglo) = -3.0*vns(ik,ie,is,2)*vpa(:,isgn,iglo) &
!                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
!          end if
       end do
    end do

    call solfp_ediffuse (bs0,init)    ! s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0s0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * bs0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v0s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 - v0s0 * z0 / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bs0(:,isgn,iglo) = bs0(:,isgn,iglo) - dtmp(:,it,ik,is) &
               * bz0(:,isgn,iglo)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1s0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = (nu_s - nu_D) vpa J0
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
               * bs0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v1s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 / (1 + v0s0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bs0(:,isgn,iglo) = bs0(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get w0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! u0 = -3 dt J1 vperp vus a f0 / du
!          if (conservative) then
          bw0(:,isgn,iglo) = -vns(ik,ie,is,1)*e(ie,is)*al(il)*aj1(:,iglo) &
               * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
               / bmag
!          else
!             bw0(:,isgn,iglo) = -3.*vns(ik,ie,is,2)*e(ie,is)*al(il)*aj1(:,iglo) &
!                  * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
!                  / bmag
!          end if
       end do
    end do

    call solfp_ediffuse (bw0,init)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0w0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * bw0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v0w0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w0 - v0w0 * z0 / (1 + v0z0) (this is w1 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bw0(:,isgn,iglo) = bw0(:,isgn,iglo) - bz0(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v1w1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = (nus-nud) vpa J0 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
               * bw0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)    ! v1w1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w1 - v1w1 * s1 / (1 + v1s1) (this is w2 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bw0(:,isgn,iglo) = bw0(:,isgn,iglo) - bs0(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v2w2

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = (nus-nud) vperp J1 f0 
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*e(ie,is)*al(il)*aj1(:,iglo) &
                  * bw0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, dtmp, all)   ! v2w2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w2 / (1 + v2w2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bw0(:,isgn,iglo) = bw0(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

    deallocate (gtmp, duinv, dtmp, vns)

# ifdef USE_LE_LAYOUT

    nxi = max(2*nlambda-1,2*ng2)

    if (.not. allocated(bs0le)) then
       allocate (bs0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (bz0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (bw0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end if
    call gather (g2le, bs0, bs0le)
    call gather (g2le, bz0, bz0le)
    call gather (g2le, bw0, bw0le)

    if (.not. allocated(aj0le)) then
       allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
       ! next set aj0le & aj1l
       z_big(:,1,:) = cmplx(aj0,aj1)
       z_big(:,2,:) = z_big(:,1,:)
       call gather (g2le, z_big, ctmp)
       allocate (aj0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (aj1le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       aj0le = real(ctmp)
       aj1le = aimag(ctmp)
       deallocate (ctmp, z_big)
    end if

# endif

  end subroutine init_diffuse_conserve

  subroutine init_vnew (hee)
    use species, only: nspec, spec, electron_species, has_electron_species
    use le_grids, only: negrid, e, w
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use run_parameters, only: zeff, tunits
    use dist_fn_arrays, only: kperp2
    use constants
    use spfunc, only: erf => erf_ext

    real, dimension (:,:), intent (out) :: hee
    real,dimension (negrid,nspec)::heevth, hsg, hsgvth
    integer :: ik, ie, is, it, ig
    real :: v, k4max
    real :: vl, vr, dv2l, dv2r
!    real :: erf ! this is needed for PGI: RN


    do is = 1, nspec
       do ie = 1, negrid
!          v = sqrt(e(ie,is))
!          hee(ie,is) = 1.0/sqrt(pi)/v*exp(-e(ie,is)) &  ! hee is vnew_D from MAB notes
!               + (1.0 - 0.5/e(ie,is)) &
!               *(1.0 - 1.0/(1.0          + v &          ! this line on is erf(v)
!               *(0.0705230784 + v &
!               *(0.0422820123 + v &
!               *(0.0092705272 + v &
!               *(0.0001520143 + v &
!               *(0.0002765672 + v &
!               *(0.0000430638)))))))**16)
          hee(ie,is) = exp(-e(ie,is))/sqrt(pi*e(ie,is)) &
               + (1.0 - 0.5/e(ie,is))*erf(sqrt(e(ie,is)))

!>MAB
! hsg is the G of Hirshman and Sigmar
! added to allow for momentum conservation with energy diffusion
!          hsg(ie,is) = -1.0/sqrt(pi)/v*exp(-e(ie,is)) &
!               + 0.5/e(ie,is) &
!               *(1.0 - 1.0/(1.0 + v &  ! this line on is erf(v)
!               *(0.0705230784 + v &
!               *(0.0422820123 + v &
!               *(0.0092705272 + v &
!               *(0.0001520143 + v &
!               *(0.0002765672 + v &
!               *(0.0000430638)))))))**16)
          hsg(ie,is) = hsg_func(sqrt(e(ie,is)))
!<MAB
       end do
    end do

!heevth is hee but using only thermal velocity (energy independent)

    do is = 1, nspec
       do ie = 1, negrid
!          v = 1           
!          heevth(ie,is) = 1.0/sqrt(pi)/v*exp(-v**2) &
!               + (1.0 - 0.5/v**2) &
!               *(1.0 - 1.0/(1.0          + v &
!               *(0.0705230784 + v &
!               *(0.0422820123 + v &
!               *(0.0092705272 + v &
!               *(0.0001520143 + v &
!               *(0.0002765672 + v &
!               *(0.0000430638)))))))**16)
          heevth(ie,is) = exp(-1.0)/sqrt(pi) &
               + 0.5*erf(1.0)

!>MAB
! hsg is the G of Helander and Sigmar
! added to allow for momentum conservation with energy diffusion
!          hsgvth(ie,is) = -1.0/sqrt(pi)/v*exp(-v**2) &
!               + 0.5/v**2 &
!               *(1.0 - 1.0/(1.0 + v &
!               *(0.0705230784 + v &
!               *(0.0422820123 + v &
!               *(0.0092705272 + v &
!               *(0.0001520143 + v &
!               *(0.0002765672 + v &
!               *(0.0000430638)))))))**16)
          hsgvth(ie,is) = hsg_func(1.0)
!<MAB
       end do
    end do                                                                  

!    do is = 1, nspec
!       if (spec(is) % nustar < 0.)
!
!    end do

    if(.not.allocated(vnew)) then
       allocate (vnew(naky,negrid,nspec))
       allocate (vnew_s(naky,negrid,nspec))
       allocate (vnew_D(naky,negrid,nspec))
       allocate (vnew_E(naky,negrid,nspec))
       allocate (delvnew(naky,negrid,nspec))
    end if
    if(.not.allocated(vnewh)) allocate (vnewh(-ntgrid:ntgrid,ntheta0,naky,nspec))

    do is = 1, nspec
       if (spec(is)%type == electron_species) then
          do ie = 1, negrid
             do ik = 1, naky
                if (const_v) then
                   vnew(ik,ie,is) = spec(is)%vnewk &
                        *(zeff+heevth(ie,is))*0.5*tunits(ik)                   
                   vnew_s(ik,ie,is) = spec(is)%vnewk &
                        *hsgvth(ie,is)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = spec(is)%vnewk &
                        *heevth(ie,is)*tunits(ik)                   
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                           - 2.0*vnew_D(ik,ie,is)
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *(zeff + hee(ie,is))*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(e(ie,is)) &
                        *hsg(ie,is)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *hee(ie,is)*tunits(ik)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = e(ie,is)*(vnew_s(ik,ie,is)*(2.0-0.5/e(ie,is)) &
                           - 2.0*vnew_D(ik,ie,is))
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                end if
                if (ei_coll_only) then
                   vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *zeff*0.5*tunits(ik)
                   vnew_s(ik,ie,is)=0.
                   vnew_D(ik,ie,is)=0.
                   vnew_E(ik,ie,is)=0.
                   delvnew(ik,ie,is)=0.
                endif
             end do
          end do
       else
          do ie = 1, negrid
             do ik = 1, naky
                if (const_v) then
                   vnew(ik,ie,is) = spec(is)%vnewk &
                        *heevth(ie,is)*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk &
                        *hsgvth(ie,is)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                           - 2.0*vnew_D(ik,ie,is)
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *hee(ie,is)*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(e(ie,is)) &
                        *hsg(ie,is)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = e(ie,is)*(vnew_s(ik,ie,is)*(2.0-0.5/e(ie,is)) &
                           - 2.0*vnew_D(ik,ie,is))
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                end if
                if (ei_coll_only) then
                   vnew(ik,ie,is) = 0.
                   vnew_s(ik,ie,is)=0.
                   vnew_D(ik,ie,is)=0.
                   vnew_E(ik,ie,is)=0.
                   delvnew(ik,ie,is)=0.
                endif
             end do
          end do
       end if

       if (conservative) then

          do ie = 2, negrid-1
             vr = 0.5*(sqrt(e(ie+1,is)) + sqrt(e(ie,is)))
             vl = 0.5*(sqrt(e(ie,is)) + sqrt(e(ie-1,is)))
             dv2r = (e(ie+1,is) - e(ie,is)) / (sqrt(e(ie+1,is)) - sqrt(e(ie,is)))
             dv2l = (e(ie,is) - e(ie-1,is)) / (sqrt(e(ie,is)) - sqrt(e(ie-1,is)))
          
             vnew_E(:,ie,is) = spec(is)%vnewk*tunits*(vl*exp(-vl**2)*dv2l*hsg_func(vl) &
                  - vr*exp(-vr**2)*dv2r*hsg_func(vr)) / (sqrt(pi)*w(ie,is))
             delvnew(:,ie,is) = spec(is)%vnewk*tunits*(vl*exp(-vl**2)*hsg_func(vl) &
                     - vr*exp(-vr**2)*hsg_func(vr)) / (sqrt(pi*e(ie,is))*w(ie,is))
          end do

          ! boundary at v = 0
          vr = 0.5*(sqrt(e(2,is)) + sqrt(e(1,is)))
          dv2r = (e(2,is) - e(1,is)) / (sqrt(e(2,is)) - sqrt(e(1,is)))

          vnew_E(:,1,is) = -spec(is)%vnewk*tunits*vr*exp(-vr**2)*hsg_func(vr)*dv2r &
               / (sqrt(pi)*w(1,is))
          delvnew(:,1,is) = -spec(is)%vnewk*tunits*vr*exp(-vr**2)*hsg_func(vr) &
               / (sqrt(pi*e(1,is))*w(1,is))

          ! boundary at v -> infinity
          vl = 0.5*(sqrt(e(negrid,is)) + sqrt(e(negrid-1,is)))
          dv2l = (e(negrid,is) - e(negrid-1,is)) / (sqrt(e(negrid,is)) - sqrt(e(negrid-1,is)))

          vnew_E(:,negrid,is) = spec(is)%vnewk*tunits*vl*exp(-vl**2)*hsg_func(vl)*dv2l &
               / (sqrt(pi)*w(negrid,is))
          delvnew(:,negrid,is) = spec(is)%vnewk*tunits*vl*exp(-vl**2)*hsg_func(vl) &
               / (sqrt(pi*e(negrid,is))*w(negrid,is))

       end if

       ! add hyper-terms inside collision operator
!BD: Warning!
!BD: For finite magnetic shear, this is different in form from what appears in hyper.f90 
!BD: because kperp2 /= akx**2 + aky**2;  there are cross terms that are dropped in hyper.f90
!BD: Warning!
!BD: Also: there is no "grid_norm" option here and the exponent is fixed to 4 for now
       if (hyper_colls) then
          k4max = (maxval(kperp2))**2 
          do ik = 1, naky
             do it = 1, ntheta0
                do ig=-ntgrid,ntgrid
                   vnewh(ig,it,ik,is) = spec(is)%nu_h * kperp2(ig,it,ik)**2/k4max
                end do
             end do
          end do
       else
          vnewh = 0.
       end if
    end do
    
  end subroutine init_vnew

  function hsg_func (vel)

    use constants, only: pi
    use spfunc, only: erf => erf_ext

    implicit none

    real, intent (in) :: vel
    real :: hsg_func

    hsg_func = 0.5*erf(vel)/vel**2-exp(-vel**2)/(sqrt(pi)*vel)

  end function hsg_func

  subroutine init_g3int
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: nlambda, al, lintegrate
    use gs2_layouts, only: g_lo, gint_lo, il_idx
    implicit none
    complex, dimension (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc) :: g
    real :: x
    integer :: il, ig, iglo

    if (.not. allocated(sq)) allocate (sq(-ntgrid:ntgrid,nlambda,2))
    do il = 1, nlambda
       do ig = -ntgrid, ntgrid
          x = sqrt(max(0.0, 1.0 - al(il)*bmag(ig)))
          sq(ig,il,1) =  x
          sq(ig,il,2) = -x
       end do
    end do

    if (.not.allocated(g3int)) &
         allocate (g3int(-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc))
    g3int = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       x = al(il_idx(g_lo,iglo))
       g(:,1,iglo) = max(0.0, 1.0 - x*bmag(:))
       g(:,2,iglo) = max(0.0, 1.0 - x*bmag(:))
    end do
    call lintegrate (g, g3int)
  end subroutine init_g3int

  subroutine init_krook (hee)
    use species, only: spec, electron_species, ion_species
    use theta_grid, only: ntgrid, eps
    use le_grids, only: integrate, ng2, al, e
    use gs2_layouts, only: g_lo, geint_lo, ik_idx, il_idx, ie_idx, is_idx
    use run_parameters, only: zeff, tunits
    use gs2_time, only: code_dt
    implicit none
    real, dimension (:,:), intent (in) :: hee
    complex, dimension (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc) :: g
    complex, dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc) &
         :: geint
    integer :: iglo, ik, il, ie, is, ige
    real :: zeff1, vep, vhat, delta00

    if (.not.allocated(aintnorm)) &
         allocate (aintnorm(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
    aintnorm = 0. ; g = 1.0
    call integrate (g, geint)

    do ige = geint_lo%llim_proc, geint_lo%ulim_proc
       aintnorm(:,ige) = 1.0/real(geint(:,ige))
    end do

    if (.not.allocated(vnewfe)) allocate (vnewfe(g_lo%llim_proc:g_lo%ulim_alloc))
    vnewfe = 0.

    if (collision_model_switch == collision_model_krook_test) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          vnewfe(iglo) = abs(spec(is)%vnewk)*tunits(ik)*code_dt
       end do
    else
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)

          if (spec(is)%type == electron_species &
               .or. spec(is)%type == ion_species) &
          then
             zeff1 = zeff
          else
             zeff1 = 0.0
          end if

          vep = abs(spec(is)%vnewk)*tunits(ik)
          vhat = sqrt(e(ie,is))
          if (ivnew > 0) then
             delta00 = absom/((vep + 2.0*epsilon(0.0))*zeff)
             if (ivnew >= 2) delta00 = (delta00*eps/37.2)**.333333333
             vnewfe(iglo) = code_dt*vep*eps*(zeff1 + hee(ie,is)) &
                  /(vhat**3*((1.0-eps-al(il))**2 + 1e-8)) &
                  *(.111*delta00+1.31)/(11.79*delta00+1.0)
          else
             vnewfe(iglo) = 0.00941/((1.0 - eps - al(il))**2 + 1e-8)
             if (il > ng2) vnewfe(iglo) = vnewfe(iglo) + vncoef/eps**2
             vnewfe(iglo) = vnewfe(iglo)*code_dt*vep*eps &
                  *(zeff1 + hee(ie,is))/vhat**3
          end if
       end do
    end if
  end subroutine init_krook

  subroutine init_ediffuse (vnmult_target)
    use constants, only: pi
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, ng2, e, ecut, integrate_moment, al, w
    use le_grids, only: vgrid, forbid
    use egrid, only: zeroes, x0, energy
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
# ifdef USE_LE_LAYOUT
    use le_grids, only: nlambda
    use gs2_layouts, only: le_lo
# else
    use gs2_layouts, only: e_lo, il_idx
# endif
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    use dist_fn_arrays, only: kperp2
    use spfunc, only: erf => erf_ext

    implicit none
    
    real, intent (in), optional :: vnmult_target

    integer :: ie, is, iglo, ik, il, ig, it
    real, dimension (:), allocatable :: aa, bb, cc, xe, el
    real :: vn, xe0, xe1, xe2, xer, xel, er, fac, ee, capgl, capgr, slb1
!    real :: erf ! this is needed for PGI: RN
# ifdef USE_LE_LAYOUT
    integer :: ile, ixi, nxi
# else
    integer :: ielo
# endif

    allocate (aa(negrid), bb(negrid), cc(negrid))
    allocate (xe(negrid))
    allocate (el(negrid))

    ! want to use x variables instead of e because we want conservative form
    ! for the x-integration
    xe(1:negrid-1) = zeroes
    xe(negrid) = x0

    if (present(vnmult_target)) then
       vnmult(2) = max (vnmult_target, 1.0)
    end if

    if (.not. vgrid) then
       do ie = 2, negrid
          xel = (xe(ie-1)+xe(ie))*0.5
          el(ie) = energy(xel,ecut)
       end do
    end if

# ifdef USE_LE_LAYOUT

    nxi = max(2*nlambda-1, 2*ng2)
    if (.not.allocated(gle)) then
       allocate (gle(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       gle = 0.0
    end if
    if (.not.allocated(ec1le)) then
       allocate (ec1le   (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (ebetaale(nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (eqle    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
       vnmult(2) = max(1.0, vnmult(2))
    end if
    ec1le = 0.0 ; ebetaale = 0.0 ; eqle = 0.0

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ik = ik_idx(le_lo, ile)
       it = it_idx(le_lo, ile)
       is = is_idx(le_lo, ile)
       ig = ig_idx(le_lo, ile)
       do ixi=1, nxi
          il = min(ixi, nxi+1-ixi)
          if (forbid(ig,il)) cycle
          call get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, el, xe)
          ec1le(ixi,:,ile) = cc
          ebetaale(ixi,1,ile) = 1.0 / bb(1)
          do ie = 1, negrid-1
             eqle(ixi,ie+1,ile) = aa(ie+1) * ebetaale(ixi,ie,ile)
             ebetaale(ixi,ie+1,ile) = 1.0 / ( bb(ie+1) - eqle(ixi,ie+1,ile) &
                  * ec1le(ixi,ie,ile) )
          end do
       end do
    end do

# else

    if (.not.allocated(ec1)) then
       allocate (ec1    (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (ebetaa (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (eql    (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       vnmult(2) = max(1.0, vnmult(2))
    endif
    ec1 = 0.0 ; ebetaa = 0.0 ; eql = 0.0

    do ielo = e_lo%llim_proc, e_lo%ulim_proc
       is = is_idx(e_lo, ielo)
       ik = ik_idx(e_lo, ielo)
       il = il_idx(e_lo, ielo)
       ig = ig_idx(e_lo, ielo)
       it = it_idx(e_lo, ielo)
       if (forbid(ig,il)) cycle
       call get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, el, xe)
       ec1(:,ielo) = cc

! fill in the arrays for the tridiagonal
       ebetaa(1,ielo) = 1.0/bb(1)
       do ie = 1, negrid-1
          eql(ie+1,ielo) = aa(ie+1)*ebetaa(ie,ielo)
          ebetaa(ie+1,ielo) = 1.0/(bb(ie+1)-eql(ie+1,ielo)*ec1(ie,ielo))
       end do

    end do

# endif

    deallocate(aa, bb, cc, xe, el)

  end subroutine init_ediffuse

  subroutine get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, el, xe)

    use species, only: spec
    use run_parameters, only: tunits
    use theta_grid, only: bmag
    use le_grids, only: al, negrid, w, vgrid
    use spfunc, only: erf => erf_ext
    use constants, only: pi
    use dist_fn_arrays, only: kperp2
    use gs2_time, only: code_dt

    implicit none

    integer, intent (in) :: ig, ik, it, il, is
    real, dimension (:), intent (in) :: el, xe
    real, dimension (:), intent (out) :: aa, bb, cc

    integer :: ie
    real :: vn, slb1, xe0, xe1, xe2, xel, xer
    real :: capgr, capgl, ee

    vn = vnmult(2)*spec(is)%vnewk*tunits(ik)

    slb1 = sqrt(max(0.0,1.0 - bmag(ig)*al(il)))     ! xi_j

    select case (ediff_switch)

    case (ediff_scheme_default)
       
       do ie = 2, negrid-1
          xe0 = xe(ie-1)
          xe1 = xe(ie)
          xe2 = xe(ie+1)
          
          xel = (xe0 + xe1)*0.5
          xer = (xe1 + xe2)*0.5
          
          if (vgrid) then
             capgr = 2.0*exp(-xer**2)*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/(xer*sqrt(pi))
             capgl = 2.0*exp(-xel**2)*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/(xel*sqrt(pi))
          else
             capgr = 8.0*xer*sqrt(el(ie+1))*exp(-2.0*el(ie+1))/pi
             capgl = 8.0*xel*sqrt(el(ie))*exp(-2.0*el(ie))/pi
          end if
          
          ee = 0.125*(1.-slb1**2)*vnew_s(ik,ie,is) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(ie) = -0.25*vn*code_dt*capgr/(w(ie,is)*(xe2 - xe1))
          aa(ie) = -0.25*vn*code_dt*capgl/(w(ie,is)*(xe1 - xe0))
          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*code_dt
          
             ! coefficients for entropy heating calculation
!             if (heating) then
!                dd(ie) =vnc*(0.25*capgr/(w(ie,is)*(xe2-xe1)) + ee)
!                hh(ie) =vnh*(0.25*capgr/(w(ie,is)*(xe2-xe1)) + ee)
!             end if
       end do

       ! boundary at v = 0
       xe1 = xe(1)
       xe2 = xe(2)
       
       xer = (xe1 + xe2)*0.5
       
       if (vgrid) then
          capgr = 2.0*exp(-xer**2)*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer/sqrt(pi)
       else
          capgr = 8.0*xer*sqrt(el(2))*exp(-2.0*el(2))/pi
       end if
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,1,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -0.25*vn*code_dt*capgr/(w(1,is)*(xe2 - xe1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*code_dt
       
!          if (heating) then
!             dd(1) =vnc*(0.25*capgr/(w(1,is)*(xe2-xe1)) + ee)
!             hh(1) =vnh*(0.25*capgr/(w(1,is)*(xe2-xe1)) + ee)
!          end if

       ! boundary at v = infinity
       xe0 = xe(negrid-1)
       xe1 = xe(negrid)
       
       xel = (xe1 + xe0)*0.5
       
       if (vgrid) then
          capgl = 2.0*exp(-xel**2)*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel/sqrt(pi)
       else
          capgl = 8.0*xel*sqrt(el(negrid))*exp(-2.0*el(negrid))/pi
       end if
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,negrid,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(negrid) = 0.0
       aa(negrid) = -0.25*vn*code_dt*capgl/(w(negrid,is)*(xe1 - xe0))
       bb(negrid) = 1.0 - aa(negrid) + ee*code_dt
       
!          if (heating) then
!             dd(negrid) =vnc*ee
!             hh(negrid) =vnh*ee
!          end if

    case (ediff_scheme_old)
       
       ! non-conservative scheme
       
       do ie = 2, negrid-1
          xe0 = xe(ie-1)
          xe1 = xe(ie)
          xe2 = xe(ie+1)
          
          xel = (xe0 + xe1)*0.5
          xer = (xe1 + xe2)*0.5
          
          if (vgrid) then
             capgr = 0.5*exp(xe1**2-xer**2)/xe1**2*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer
             capgl = 0.5*exp(xe1**2-xel**2)/xe1**2*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel
          else
             capgr = 8.0*xer*sqrt(el(ie+1))*exp(-2.0*el(ie+1))/pi
             capgl = 8.0*xel*sqrt(el(ie))*exp(-2.0*el(ie))/pi
          end if
          
          ee = 0.125*(1.-slb1**2)*vnew_s(ik,ie,is) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(ie) = -vn*code_dt*capgr/((xer-xel)*(xe2 - xe1))
          aa(ie) = -vn*code_dt*capgl/((xer-xel)*(xe1 - xe0))
          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*code_dt
          
          ! coefficients for entropy heating calculation
!             if (heating) then
!                dd(ie) =vnc*(0.25*capgr/(w(ie,is)*(xe2-xe1)) + ee)
!                hh(ie) =vnh*(0.25*capgr/(w(ie,is)*(xe2-xe1)) + ee)
!             end if
       end do

       ! boundary at xe = 0
       xe1 = xe(1)
       xe2 = xe(2)
       
       xer = (xe1 + xe2)*0.5
       
       if (vgrid) then
          capgr = 0.5*exp(xe1**2-xer**2)/xe1**2*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer
       else
          capgr = 8.0*xer*sqrt(el(2))*exp(-2.0*el(2))/pi
       end if
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,1,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -vn*code_dt*capgr/(xer*(xe2 - xe1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*code_dt
       
!          if (heating) then
!             dd(1) =vnc*(0.25*capgr/(w(1,is)*(xe2-xe1)) + ee)
!             hh(1) =vnh*(0.25*capgr/(w(1,is)*(xe2-xe1)) + ee)
!          end if

       ! boundary at xe = 1
       xe0 = xe(negrid-1)
       xe1 = xe(negrid)
       
       xel = (xe1 + xe0)*0.5
       
       if (vgrid) then
          capgl = 0.5*exp(xe1**2-xel**2)/xe1**2*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel
       else
          capgl = 8.0*xel*sqrt(el(negrid))*exp(-2.0*el(negrid))/pi
       end if
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,negrid,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(negrid) = 0.0
       aa(negrid) = -vn*code_dt*capgl/((1.0-xel)*(xe1 - xe0))
       bb(negrid) = 1.0 - aa(negrid) + ee*code_dt

!          if (heating) then
!             dd(negrid) =vnc*ee
!             hh(negrid) =vnh*ee
!          end if

    end select

  end subroutine get_ediffuse_matrix

  subroutine init_ediffuse_redistribute
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, ng2
    use gs2_layouts, only: init_ediffuse_layouts
    use gs2_layouts, only: g_lo, e_lo, ie_idx, ik_idx, it_idx, il_idx, is_idx
    use gs2_layouts, only: idx_local, proc_id, idx
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, ik, it, il, ie, is, ielo
    integer :: n, ip
!    logical :: ediffinit = .false.

    if (ediffinit) return
    
!DD, March 2009: Nullify pointers on initialisation, so do not pass association
!                test when not allocated. 
!  Problem arose on Pascali compiler (York) with  to_list(ip)%third and fourth
    do ip=0, (nproc-1)
       nullify(to_list(ip)%first,from_list(ip)%first,to_list(ip)%second,from_list(ip)%second,to_list(ip)%third,from_list(ip)%third,to_list(ip)%fourth,from_list(ip)%fourth)
    end do
!<DD
	
!DD, March 2009: Nullify pointers on initialisation, so do not pass association
!                test when not allocated. 
!  Problem arose on Pascali compiler (York) with  to_list(ip)%third and fourth
    do ip=0, (nproc-1)
       nullify(to_list(ip)%first,from_list(ip)%first,to_list(ip)%second,from_list(ip)%second,to_list(ip)%third,from_list(ip)%third,to_list(ip)%fourth,from_list(ip)%fourth)
    end do
!<DD

    call init_ediffuse_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             ielo = idx(e_lo,ig,isign,ik,it,il,is)
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(e_lo,ielo)) = nn_from(proc_id(e_lo,ielo)) + 1
             if (idx_local(e_lo,ielo)) &
                  nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
          end do
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             ielo = idx(e_lo,ig,isign,ik,it,il,is)

             if (idx_local(g_lo,iglo)) then
                ip = proc_id(e_lo,ielo)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end if
             if (idx_local(e_lo,ielo)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = ie
                to_list(ip)%second(n) = ielo
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = e_lo%llim_proc

    to_high(1) = negrid+1
    to_high(2) = e_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc
    

    call init_redist (ediffuse_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

    ediffinit = .true.

  end subroutine init_ediffuse_redistribute

  subroutine init_lorentz (vnmult_target)
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, al, jend, ng2, e, wl
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: kperp2
# ifdef USE_LE_LAYOUT
    use gs2_layouts, only: le_lo
# else
    use gs2_layouts, only: lz_lo
# endif
    use gs2_layouts, only: ig_idx, ik_idx, ie_idx, is_idx, it_idx

    implicit none

    real, intent (in), optional :: vnmult_target

    integer :: ig, il, it, ik, ie, is, je, te, te2, teh, nxi
    real, dimension (:), allocatable :: aa, bb, cc
    real, dimension (:), allocatable :: dd, hh
    real :: slb0, slb1, slb2, slbl, slbr, vn, ee, vnh, vnc
# ifdef USE_LE_LAYOUT
    integer :: ile
# else
    integer :: ilz
# endif

    nxi = max(2*nlambda-1, 2*ng2)

    allocate (aa(nxi+1), bb(nxi+1), cc(nxi+1), dd(nxi+1), hh(nxi+1))

    call init_vpdiff

# ifdef USE_LE_LAYOUT

    if (.not. allocated(gle)) then
       allocate (gle(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       gle = 0.0
       if (heating) then
          allocate (glec(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
          glec = 0.0
       end if
    end if
    if (.not.allocated(c1le)) then
       allocate (c1le    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (betaale (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (qle     (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
       if (heating) then
          allocate (d1le    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (h1le    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          d1le = 0.0 ; h1le = 0.0
       end if
       vnmult(1) = max(1.0, vnmult(1))
    endif
    c1le = 0.0 ; betaale = 0.0 ; qle = 0.0

    if (present(vnmult_target)) then
       vnmult(1) = max (vnmult_target, 1.0)
    end if

    do ile = le_lo%llim_proc, le_lo%ulim_proc

       ig = ig_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       is = is_idx(le_lo,ile)
       je = jend(ig)
       if (je <= ng2+1) then
          te2 = 2*ng2
       else
          te2 = 2*je-1
       end if
       do ie = 1, negrid

          call get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)

          c1le(:,ie,ile) = cc
          if (allocated(d1le)) then
             d1le(:,ie,ile) = dd
             h1le(:,ie,ile) = hh
          end if

          qle(1,ie,ile) = 0.0
          betaale(1,ie,ile) = 1.0/bb(1)
          do il = 1, te2-1
             qle(il+1,ie,ile) = aa(il+1) * betaale(il,ie,ile)
             betaale(il+1,ie,ile) = 1.0 / ( bb(il+1) - qle(il+1,ie,ile) &
                  * c1le(il,ie,ile) )
          end do
          qle(te2+1,ie,ile) = 0.0
          betaale(te2+1,ie,ile) = 0.0

       end do
    end do

# else

    if (.not.allocated(c1)) then
       allocate (c1(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (betaa(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (ql(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
       if (heating) then
          allocate (d1   (nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
          allocate (h1   (nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
          d1 = 0.0 ; h1 = 0.0
       end if
       vnmult(1) = max(1.0, vnmult(1))
    end if
       
    c1 = 0.0 ; betaa = 0.0 ; ql = 0.0

    if (present(vnmult_target)) then
       vnmult(1) = max (vnmult_target, 1.0)
    end if

    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       is = is_idx(lz_lo,ilz)
       ik = ik_idx(lz_lo,ilz)
       it = it_idx(lz_lo,ilz)
       ie = ie_idx(lz_lo,ilz)
       ig = ig_idx(lz_lo,ilz)
       je = jend(ig)
       if (je <= ng2+1) then
          te2 = 2*ng2
       else
          te2 = 2*je-1
       end if
       ! tested get_lorentz_matrix -- MAB
       call get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)
       c1(:,ilz) = cc
       if (allocated(d1)) then
          d1(:,ilz) = dd
          h1(:,ilz) = hh
       end if

       betaa(1,ilz) = 1.0/bb(1)
       do il = 1, te2-1
          ql(il+1,ilz) = aa(il+1)*betaa(il,ilz)
          betaa(il+1,ilz) = 1.0/(bb(il+1)-ql(il+1,ilz)*c1(il,ilz))
       end do

       ql(1,ilz) = 0.0
       ql(te2+1:,ilz) = 0.0
       c1(te2+1:,ilz) = 0.0
       betaa(te2+1:,ilz) = 0.0

    end do

# endif

    deallocate (aa, bb, cc, dd, hh)

  end subroutine init_lorentz

  subroutine get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)

    use species, only: spec
    use le_grids, only: nlambda, al, e, ng2
    use le_grids, only: wl, jend, al
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: kperp2
    use theta_grid, only: bmag, ntgrid
    use run_parameters, only: tunits

    implicit none

    real, dimension (:), intent (out) :: aa, bb, cc, dd, hh
    integer, intent (in) :: ig, ik, it, ie, is

    integer :: il, je, te, te2, teh
    real :: slb0, slb1, slb2, slbl, slbr, vn, vnh, vnc, ee

    je = jend(ig)
    if (je <= ng2+1) then
       te = ng2
       te2 = 2*ng2
       teh = ng2
    else
       te = je
       te2 = 2*je-1
       teh = je-1
    end if
    if (collision_model_switch == collision_model_lorentz_test) then
       vn = vnmult(1)*abs(spec(is)%vnewk)*tunits(ik)
       vnc = 0.
       vnh = 0.
    else
       if (hypermult) then
          vn = vnmult(1)*vnew(ik,ie,is)*(1.+vnewh(ig,it,ik,is))
          vnc = vnmult(1)*vnew(ik,ie,is)
          vnh = vnewh(ig,it,ik,is)*vnc
       else
          vn = vnmult(1)*vnew(ik,ie,is)+vnewh(ig,it,ik,is)
          vnc = vnmult(1)*vnew(ik,ie,is)
          vnh = vnewh(ig,it,ik,is)
       end if
    end if

    aa = 0.0 ; bb = 0.0 ; cc = 0.0 ; dd = 0.0 ; hh = 0.0

    select case (lorentz_switch)
       
    case (lorentz_scheme_default)          
       
       do il = 2, te-1
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
          slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1))) 
          
          slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
          slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
          
          ee = 0.5*e(ie,is)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(il) = 2.0*vn*code_dt*(1.0 - slbr**2)/(wl(ig,il)*(slb2 - slb1))
          aa(il) = 2.0*vn*code_dt*(1.0 - slbl**2)/(wl(ig,il)*(slb1 - slb0))
          bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
          
          ! coefficients for entropy heating calculation
          if (heating) then
             dd(il) =vnc*(-2.0*(1.0-slbr*slbr)/(wl(ig,il)*(slb2-slb1)) + ee)
             hh(il) =vnh*(-2.0*(1.0-slbr*slbr)/(wl(ig,il)*(slb2-slb1)) + ee)
          end if
       end do
       
       ! boundary at xi = 1
       slb1 = sqrt(abs(1.0-bmag(ig)*al(1)))
       slb2 = sqrt(abs(1.0-bmag(ig)*al(2)))
       
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*e(ie,is)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = 2.0*vn*code_dt*(1.0 - slbr**2)/(wl(ig,1)*(slb2 - slb1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*vn*code_dt
       
       if (heating) then
          dd(1) =vnc*(-2.0*(1.0-slbr**2)/(wl(ig,1)*(slb2-slb1)) + ee)
          hh(1) =vnh*(-2.0*(1.0-slbr**2)/(wl(ig,1)*(slb2-slb1)) + ee)
       end if
       
       ! boundary at xi = 0
       il = te
       slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
       if (te == ng2) then
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
          slb2 = -slb1
       else
          slb1 = 0.0
          slb2 = -slb0
       end if
       
       slbl = (slb1 + slb0)*0.5
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*e(ie,is)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(il) = 2.0*vn*code_dt*(1.0 - slbr**2)/(wl(ig,il)*(slb2 - slb1))
       aa(il) = 2.0*vn*code_dt*(1.0 - slbl**2)/(wl(ig,il)*(slb1 - slb0))
       bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
       
       if (heating) then
          dd(il) =vnc*(-2.0*(1.0-slbr**2)/(wl(ig,il)*(slb2-slb1)) + ee)
          hh(il) =vnh*(-2.0*(1.0-slbr**2)/(wl(ig,il)*(slb2-slb1)) + ee)
       end if
       
    case (lorentz_scheme_old)
       
       do il = 2, te-1
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
          slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1))) 
          
          slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
          slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
          
          ee = 0.5*e(ie,is)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
          aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
          
          ! coefficients for entropy heating calculation
          if (heating) then
             dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
             hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          end if
       end do
       
       ! boundary at xi = 1
       slb0 = 1.0
       slb1 = sqrt(abs(1.0-bmag(ig)*al(1)))
       slb2 = sqrt(abs(1.0-bmag(ig)*al(2)))
       
       slbl = (slb1 + slb0)*0.5
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*e(ie,is)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -vn*code_dt*(-1.0 - slbr)/(slb2-slb1)
       aa(1) = 0.0
       bb(1) = 1.0 - (aa(1) + cc(1)) + ee*vn*code_dt
       
       if (heating) then
          dd(1) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(1) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
       end if
       
       ! boundary at xi = 0
       il = te
       slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
       if (te == ng2) then
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
          slb2 = -slb1
       else
          slb1 = 0.0
          slb2 = -slb0
       end if
       
       slbl = (slb1 + slb0)*0.5
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*e(ie,is)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
       aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
       bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
       
       if (heating) then
          dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
       end if
       
    end select

    ! assuming symmetry in xi, fill in the rest of the arrays.
    aa(te+1:te2) = cc(teh:1:-1)
    bb(te+1:te2) = bb(teh:1:-1)
    cc(te+1:te2) =aa(teh:1:-1)
       
    if (heating) then
       dd(te+1:te2) = dd(teh:1:-1)
       hh(te+1:te2) = hh(teh:1:-1)
    end if

  end subroutine get_lorentz_matrix

  subroutine init_lorentz_error

    use le_grids, only: jend, al, ng2, nlambda
    use theta_grid, only: ntgrid, bmag
    implicit none
    
    integer :: je, ig, il, ip, ij, im
    real :: slb0, slb1, slb2, slbr, slbl
    real, dimension (:), allocatable :: slb
    real, dimension (:,:), allocatable :: dprod
    real, dimension (:,:,:), allocatable :: dlcoef, d2lcoef

    allocate(slb(2*nlambda))
    allocate (dprod(nlambda,5))

    allocate (dlcoef(-ntgrid:ntgrid,nlambda,5))
    allocate (d2lcoef(-ntgrid:ntgrid,nlambda,5))
    allocate (dtot(-ntgrid:ntgrid,nlambda,5))
    allocate (fdf(-ntgrid:ntgrid,nlambda), fdb(-ntgrid:ntgrid,nlambda))

    dlcoef = 1.0; d2lcoef = 0.0; dtot = 0.0
    fdf = 0.0; fdb = 0.0; slb = 0.0

    do ig=-ntgrid,ntgrid
       je = jend(ig)
       
       if (je <= ng2+1) then            ! no trapped particles

! calculation of xi and finite difference coefficients for non-boundary points
          do il=2,ng2-1
             slb(il) = sqrt(abs(1.0-al(il)*bmag(ig)))   ! xi_{j}
             
             slb2 = sqrt(abs(1.0-al(il+1)*bmag(ig)))    ! xi_{j+1}
             slb1 = slb(il)
             slb0 = sqrt(abs(1.0-al(il-1)*bmag(ig)))    ! xi_{j-1}
             
             slbr = (slb2+slb1)*0.5                     ! xi_{j+1/2}
             slbl = (slb1+slb0)*0.5                     ! xi_{j-1/2}

! finite difference coefficients
             fdf(ig,il) = (1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             fdb(ig,il) = (1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          end do

! boundary at xi = 1
          slb(1) = sqrt(abs(1.0-al(1)*bmag(ig)))
          slb0 = 1.0
          slb1 = slb(1)
          slb2 = slb(2)

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

! derivative of [(1-xi**2)*df/dxi] at xi_{j=1} is centered, with upper xi=1 and
! lower xi = xi_{j+1/2}
          fdf(ig,1) = (-1.0-slbr)/(slb2-slb1)
          fdb(ig,1) = 0.0

! boundary at xi = 0
          il = ng2
          slb(il) = sqrt(abs(1.0 - al(il)*bmag(ig)))
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))
          slb1 = slb(il)
          slb2 = -slb1

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

          fdf(ig,il) = (1.0 - slbr*slbr)/(slbr-slbl)/(slb2-slb1)
          fdb(ig,il) = (1.0 - slbl*slbl)/(slbr-slbl)/(slb1-slb0)

          slb(ng2+1:) = -slb(ng2:1:-1)
       else          ! run with trapped particles
          do il=2,je-1
             slb(il) = sqrt(abs(1.0-al(il)*bmag(ig)))
             
             slb2 = sqrt(abs(1.0-al(il+1)*bmag(ig)))
             slb1 = slb(il)
             slb0 = sqrt(abs(1.0-al(il-1)*bmag(ig)))
             
             slbr = (slb2+slb1)*0.5
             slbl = (slb1+slb0)*0.5

             fdf(ig,il) = (1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             fdb(ig,il) = (1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          end do

! boundary at xi = 1
          slb(1) = sqrt(abs(1.0-bmag(ig)*al(1)))
          slb0 = 1.0
          slb1 = slb(1)
          slb2 = slb(2)

          slbr = (slb1 + slb2)/2.0

          fdf(ig,1) = (-1.0 - slbr)/(slb2-slb1)
          fdb(ig,1) = 0.0

! boundary at xi = 0
          il = je
          slb(il) = sqrt(abs(1.0-bmag(ig)*al(il)))
          slb0 = slb(je-1)
          slb1 = 0.
          slb2 = -slb0                                                        

          slbl = (slb1 + slb0)/2.0

          fdf(ig,il) = (1.0 - slbl*slbl)/slb0/slb0
          fdb(ig,il) = fdf(ig,il)

          slb(je+1:2*je-1) = -slb(je-1:1:-1)
       end if

! compute coefficients (dlcoef) multipyling first derivative of h
       do il=3,ng2
          do ip=il-2,il+2
             if (il == ip) then
                dlcoef(ig,il,ip-il+3) = 0.0
                do ij=il-2,il+2
                   if (ij /= ip) dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3) + 1/(slb(il)-slb(ij))
                end do
             else
                do ij=il-2,il+2
                   if (ij /= ip .and. ij /= il) then
                      dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                   end if
                end do
                dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
             end if
             dlcoef(ig,il,ip-il+3) = -2.0*slb(il)*dlcoef(ig,il,ip-il+3)
          end do
       end do

       il = 1
       do ip=il,il+2
          if (il == ip) then
             dlcoef(ig,il,ip) = 0.0
             do ij=il,il+2
                if (ij /= ip) dlcoef(ig,il,ip) = dlcoef(ig,il,ip) + 1./(slb(il)-slb(ij))
             end do
          else
             do ij=il,il+2
                if (ij /= ip .and. ij /= il) then
                   dlcoef(ig,il,ip) = dlcoef(ig,il,ip)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do
             dlcoef(ig,il,ip) = dlcoef(ig,il,ip)/(slb(ip)-slb(il))
          end if
          dlcoef(ig,il,ip) = -2.0*slb(il)*dlcoef(ig,il,ip)
       end do

       il = 2
       do ip=il-1,il+1
          if (il == ip) then
             dlcoef(ig,il,ip-il+2) = 0.0
             do ij=il-1,il+1
                if (ij /= ip) dlcoef(ig,il,ip-il+2) = dlcoef(ig,il,ip-il+2) + 1/(slb(il)-slb(ij))
             end do
          else
             do ij=il-1,il+1
                if (ij /= ip .and. ij /= il) then
                   dlcoef(ig,il,ip-il+2) = dlcoef(ig,il,ip-il+2)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do
             dlcoef(ig,il,ip-il+2) = dlcoef(ig,il,ip-il+2)/(slb(ip)-slb(il))
          end if
          dlcoef(ig,il,ip-il+2) = -2.0*slb(il)*dlcoef(ig,il,ip-il+2)
       end do

       dprod = 2.0

! compute coefficients (d2lcoef) multiplying second derivative of h
       do il=3,ng2
          do ip=il-2,il+2
             if (il == ip) then
                do ij=il-2,il+2
                   if (ij /= ip) then
                      do im=il-2,il+2
                         if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip-il+3) = &
                              d2lcoef(ig,il,ip-il+3) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                      end do
                   end if
                end do
             else
                do ij=il-2,il+2
                   if (ij /= il .and. ij /= ip) then
                      dprod(il,ip-il+3) = dprod(il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                   end if
                end do

                do ij=il-2,il+2
                   if (ij /= ip .and. ij /= il) then
                      d2lcoef(ig,il,ip-il+3) = d2lcoef(ig,il,ip-il+3) + 1./(slb(il)-slb(ij))
                   end if
                end do
                d2lcoef(ig,il,ip-il+3) = dprod(il,ip-il+3) &
                     *d2lcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
             end if
             d2lcoef(ig,il,ip-il+3) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip-il+3)
          end do
       end do

       il = 1
       do ip=il,il+2
          if (il == ip) then
             do ij=il,il+2
                if (ij /= ip) then
                   do im=il,il+2
                      if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                   end do
                end if
             end do
          else
             do ij=il,il+2
                if (ij /= il .and. ij /= ip) then
                   dprod(il,ip) = dprod(il,ip)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do

             do ij=il,il+2
                if (ij /= ip .and. ij /= il) then
                   d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./(slb(il)-slb(ij))
                end if
             end do
             d2lcoef(ig,il,ip) = dprod(il,ip)*d2lcoef(ig,il,ip)/(slb(ip)-slb(il))
          end if
          d2lcoef(ig,il,ip) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip)
       end do

       il = 2
       do ip=il-1,il+1
          if (il == ip) then
             do ij=il-1,il+1
                if (ij /= ip) then
                   do im=il-1,il+1
                      if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip-il+2) &
                           = d2lcoef(ig,il,ip-il+2) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                   end do
                end if
             end do
          else
             do ij=il-1,il+1
                if (ij /= il .and. ij /= ip) then
                   dprod(il,ip-il+2) = dprod(il,ip-il+2)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do

             do ij=il-1,il+1
                if (ij /= ip .and. ij /= il) then
                   d2lcoef(ig,il,ip-il+2) = d2lcoef(ig,il,ip-il+2) + 1./(slb(il)-slb(ij))
                end if
             end do
             d2lcoef(ig,il,ip-il+2) = dprod(il,ip-il+2)*d2lcoef(ig,il,ip-il+2)/(slb(ip)-slb(il))
          end if
          d2lcoef(ig,il,ip-il+2) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip-il+2)
       end do
       
       if (je > ng2+1) then      ! have to handle trapped particles

          do il=ng2+1,je
             do ip=il-2,il+2
                if (il == ip) then
                   dlcoef(ig,il,ip-il+3) = 0.0
                   do ij=il-2,il+2
                      if (ij /= ip) dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3) + 1/(slb(il)-slb(ij))
                   end do
                else
                   do ij=il-2,il+2
                      if (ij /= ip .and. ij /= il) then
                         dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                      end if
                   end do
                   dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
                end if
                dlcoef(ig,il,ip-il+3) = -2.0*slb(il)*dlcoef(ig,il,ip-il+3)
             end do
          end do

          do il=ng2+1,je
             do ip=il-2,il+2
                if (il == ip) then
                   do ij=il-2,il+2
                      if (ij /= ip) then
                         do im=il-2,il+2
                            if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip-il+3) = &
                                 d2lcoef(ig,il,ip-il+3) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                         end do
                      end if
                   end do
                else
                   do ij=il-2,il+2
                      if (ij /= il .and. ij /= ip) then
                         dprod(il,ip-il+3) = dprod(il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                      end if
                   end do
                   
                   do ij=il-2,il+2
                      if (ij /= ip .and. ij /= il) then
                         d2lcoef(ig,il,ip-il+3) = d2lcoef(ig,il,ip-il+3) + 1./(slb(il)-slb(ij))
                      end if
                   end do
                   d2lcoef(ig,il,ip-il+3) = dprod(il,ip-il+3) &
                        *d2lcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
                end if
                d2lcoef(ig,il,ip-il+3) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip-il+3)
             end do
          end do
          
       end if
    end do

    dtot = dlcoef + d2lcoef

    deallocate (slb, dprod, dlcoef, d2lcoef)
  end subroutine init_lorentz_error

  subroutine init_lorentz_redistribute
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, jend, ng2
    use gs2_layouts, only: init_lorentz_layouts
    use gs2_layouts, only: g_lo, lz_lo
    use gs2_layouts, only: idx_local, proc_id
! TT>
!!$    use gs2_layouts, only: gidx2lzidx
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx
! <TT
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ilz
! TT>
    integer :: ik, it, ie, is, je
! <TT
    integer :: n, ip
!    logical :: done = .false.

    if (lzinit) return

    call init_lorentz_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
! TT> expansion of gidx2lzidx
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
! <TT
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
! TT>
!!$             call gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
             ilz = idx(lz_lo, ig, ik, it, ie, is)
! <TT
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(lz_lo,ilz)) = nn_from(proc_id(lz_lo,ilz)) + 1
             if (idx_local(lz_lo,ilz)) &
                  nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
          end do
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
! TT> expansion of gidx2lzidx
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
! <TT
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
! TT>
!!$             call gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
! TT: following just the inline expansion of gidx2lzidx except that
! TT: the arguments of idx function is replaced by variables
             je = jend(ig)
             il = il_idx(g_lo,iglo)
             if (je == 0) then
                if (isign == 2) then
                   il = 2*g_lo%nlambda+1 - il
                end if
             else
                if (il == je) then
                   if (isign == 1) il = 2*je  ! throw this info away
                else if (il > je) then 
                   if (isign == 1) il = il + je
                   if (isign == 2) il = 2*g_lo%nlambda + 1 - il + je
                else
                   if (isign == 2) il = 2*je - il !+ 1
                end if
             end if
             ilz = idx(lz_lo, ig, ik, it, ie, is)
! <TT
!             write(*,*) ig,':',isign,':',iglo,':',il,':',ilz
             if (idx_local(g_lo,iglo)) then
                ip = proc_id(lz_lo,ilz)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end if
             if (idx_local(lz_lo,ilz)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = il
                to_list(ip)%second(n) = ilz
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = lz_lo%llim_proc

    to_high(1) = max(2*nlambda, 2*ng2+1)
    to_high(2) = lz_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc
    call init_redist (lorentz_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)
    call delete_list (to_list)
    call delete_list (from_list)

    lzinit = .true.

  end subroutine init_lorentz_redistribute

  subroutine init_g2le_redistribute
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, jend, ng2
    use gs2_layouts, only: init_le_layouts
    use gs2_layouts, only: g_lo, le_lo
    use gs2_layouts, only: idx_local, proc_id
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, il_idx, idx
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (3) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ile
    integer :: ik, it, ie, is
    integer :: n, ip, je

    if (leinit) return

    call init_le_layouts (ntgrid, naky, ntheta0, nspec)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
!       do isign = 1, 2
       do ig = -ntgrid, ntgrid
          ile = idx (le_lo, ig, ik, it, is)
          if (idx_local(g_lo,iglo)) &
!               nn_from(proc_id(le_lo,ile)) = nn_from(proc_id(le_lo,ile)) + 1
               nn_from(proc_id(le_lo,ile)) = nn_from(proc_id(le_lo,ile)) + 2
          if (idx_local(le_lo,ile)) &
!               nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
               nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 2
       end do
!       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first (nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third (nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first (nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third (nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
!       il = il_idx(g_lo,iglo)
       do isign = 1, 2
!          if (isign==2) il = 2*g_lo%nlambda+1 - il
          do ig = -ntgrid, ntgrid
             il = il_idx(g_lo,iglo)
             je = jend(ig)
             if (je == 0) then
                if (isign == 2) then
                   il = 2*g_lo%nlambda+1 - il
                end if
             else
                if (il == je) then
                   if (isign == 1) il = 2*je  ! throw this info away
                else if (il > je) then 
                   if (isign == 1) il = il + je
                   if (isign == 2) il = 2*g_lo%nlambda + 1 - il + je
                else
                   if (isign == 2) il = 2*je - il !+ 1
                end if
             end if
             ile = idx (le_lo, ig, ik, it, is)
             if (idx_local(g_lo,iglo)) then
                ip = proc_id(le_lo,ile)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n)  = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n)  = iglo
             end if
             if (idx_local(le_lo,ile)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n)  = il
                to_list(ip)%second(n) = ie
                to_list(ip)%third(n)  = ile
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    to_low = le_lo%llim_proc

    to_high(1) = max(2*nlambda, 2*ng2+1)
    to_high(2) = negrid + 1  ! TT: just followed convention with +1.
    ! TT: It may be good to avoid bank conflict.
    to_high(3) = le_lo%ulim_alloc

    call init_redist (g2le, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

    leinit = .true.

  end subroutine init_g2le_redistribute

  subroutine solfp1 (g, gold, g1, phi, apar, bpar, phinew, aparnew, bparnew, diagnostics)
    use gs2_layouts, only: g_lo, it_idx, ik_idx, ie_idx, is_idx
    use theta_grid, only: ntgrid
    use run_parameters, only: tunits, beta
    use run_parameters, only: fphi, fbpar
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    use le_grids, only: e, integrate_moment
    use species, only: nspec, spec, electron_species
    use dist_fn_arrays, only: c_rate, vpa, kperp2, aj0
    use run_parameters, only: ieqzip

    use constants

    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gold, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar, phinew, aparnew, bparnew
    complex, dimension (:,:,:), allocatable :: gc1, gc2, gc3
    integer, optional, intent (in) :: diagnostics

    integer :: ig, it, ik, ie, is, iglo

    ! TMP FOR TESTING -- MAB
!    integer :: t0, t1, t2, t3, t4, t5, tr
!    real :: t1tot = 0., t2tot = 0., t3tot = 0., t4tot = 0., t5tot = 0.

    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, upar, uperp, ttot
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntotnew, uparnew, uperpnew, ttotnew

    if (collision_model_switch == collision_model_none) return

    if (heating .and. present(diagnostics)) then
       allocate (gc1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; gc1 = 0.
       allocate (gc2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; gc2 = 0.
       allocate (gc3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; gc3 = 0.
    else
       allocate (gc1(1,1,1)) ; gc1 = 0.
       allocate (gc2(1,1,1)) ; gc2 = 0.
       allocate (gc3(1,1,1)) ; gc3 = 0.
    end if
    
    ! changed 28/07/08 -- MAB
!    if (adjust) call g_adjust (g, phi, bpar, fphi, fbpar)
    if (adjust) call g_adjust (g, phinew, bparnew, fphi, fbpar)
    if (adjust) call g_adjust (gold, phi, bpar, fphi, fbpar)

    ! added 11/05/08 -- MAB
    if (heating .and. present(diagnostics)) gc3 = g

    select case (collision_model_switch)
    case (collision_model_full)

       ! TMP FOR TESTING -- MAB
!       call system_clock (count=t0, count_rate=tr)

       call solfp_ediffuse (g)

       ! TMP FOR TESTING -- MAB
!       call system_clock (count=t1)
!       t1tot = t1tot + real(t1-t0)/tr

       if (conserve_moments) call conserve_diffuse (g, g1)

       ! TMP FOR TESTING -- MAB
!       call system_clock (count=t2)
!       t2tot = t2tot + real(t2-t1)/tr

       if (resistivity .and. beta > epsilon(0.0) .and. nspec > 1) then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             is = is_idx(g_lo,iglo)
             if (spec(is)%type /= electron_species) cycle
             it = it_idx(g_lo,iglo)
             ik = ik_idx(g_lo,iglo)
             ie = ie_idx(g_lo,iglo)
             do ig = -ntgrid, ntgrid
                g(ig,:,iglo) = g(ig,:,iglo) + ieqzip(it,ik) * &
                     vnmult(1)*spec(is)%vnewk*code_dt &
                     * vpa(ig,:,iglo)*kperp2(ig,it,ik)*aparnew(ig,it,ik)*aj0(ig,iglo) &
                     / (beta*spec(is)%stm*e(ie,is)**1.5)
                ! probably need 1/(spec(is_ion)%z*spec(is_ion)%dens) above
             end do
          end do
       end if

       ! TMP FOR TESTING -- MAB
!       call system_clock (count=t3)
!       t3tot = t3tot + real(t3-t2)/tr
       
       if (heating .and. present(diagnostics)) then
          call solfp_lorentz (g, gc1, gc2, diagnostics)
       else
          call solfp_lorentz (g, gc1, gc2)
       end if

       ! TMP FOR TESTING -- MAB
!       call system_clock (count=t4)
!       t4tot = t4tot + real(t4-t3)/tr

       if (conserve_moments) call conserve_lorentz (g, g1)

       ! TMP FOR TESTING -- MAB
!       call system_clock (count=t5)
!       t5tot = t5tot + real(t5-t4)/tr

!       write (*,*) t1tot, t2tot, t3tot, t4tot, t5tot

    case (collision_model_lorentz,collision_model_lorentz_test)

       if (resistivity .and. beta > epsilon(0.0) .and. nspec > 1) then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             is = is_idx(g_lo,iglo)
             if (spec(is)%type /= electron_species) cycle
             it = it_idx(g_lo,iglo)
             ik = ik_idx(g_lo,iglo)
             ie = ie_idx(g_lo,iglo)
             do ig = -ntgrid, ntgrid
                g(ig,:,iglo) = g(ig,:,iglo) + ieqzip(it,ik) * &
                     vnmult(1)*spec(is)%vnewk*code_dt &
                     * vpa(ig,:,iglo)*kperp2(ig,it,ik)*aparnew(ig,it,ik)*aj0(ig,iglo) &
                     / (beta*spec(is)%stm*e(ie,is)**1.5)
             end do
          end do
       end if
       
       if (heating .and. present(diagnostics)) then
          call solfp_lorentz (g, gc1, gc2, diagnostics)
       else
          call solfp_lorentz (g, gc1, gc2)
       end if

       if (conserve_moments) call conserve_lorentz (g, g1)

    case (collision_model_ediffuse)

       call solfp_ediffuse (g)
       if (conserve_moments) call conserve_diffuse (g, g1)

    case (collision_model_krook,collision_model_krook_test)
       call solfp_krook (g, g1)
    end select
    
    if (heating .and. present(diagnostics)) then

       call integrate_moment (gc1, c_rate(:,:,:,:,1))
       deallocate (gc1)

       if (hyper_colls) call integrate_moment (gc2, c_rate(:,:,:,:,2))
       deallocate (gc2)

! form (h_i+1 + h_i)/2 * C(h_i+1) and integrate.  

       gc3 = 0.5*conjg(g+gold)*(g-gc3)/code_dt

       call integrate_moment (gc3, c_rate(:,:,:,:,3))

       deallocate (gc3)
    end if

    ! changed 28/07/08 -- MAB
!    if (adjust) call g_adjust (g, phi, bpar, -fphi, -fbpar)
    if (adjust) call g_adjust (g, phinew, bparnew, -fphi, -fbpar)
    if (adjust) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
    
  end subroutine solfp1

  subroutine conserve_lorentz (g, g1)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: e, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: ieqzip
# ifdef USE_LE_LAYOUT
    use le_grids, only: nlambda, negrid, e, al, ng2
    use gs2_layouts, only: le_lo, ig_idx
    use redistribute, only: scatter
    use run_parameters, only: tunits
    use theta_grid, only: bmag
# endif    

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns
# ifdef USE_LE_LAYOUT
    real, dimension (:,:,:,:), allocatable :: vpanud
    complex, dimension (:), allocatable :: v0y0, v1y1, v2y2
# else
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y1, v2y2
# endif

    integer :: ig, isgn, iglo, ik, ie, il, is, it, all = 1
# ifdef USE_LE_LAYOUT
    integer :: ile, ixi, nxi
# endif


# ifdef USE_LE_LAYOUT

    allocate (v0y0(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v1y1(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v2y2(le_lo%llim_proc:le_lo%ulim_alloc))

    nxi = max(2*nlambda-1,2*ng2)
    ! Let's work on gle directly instead of g for the moment
    allocate (gtmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc)) ; gtmp = 0.0
    allocate (vpanud(-ntgrid:ntgrid, nxi+1, negrid+1, nspec)) ; vpanud = 0.0

    if (resistivity) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

       isgn = 1
       do is = 1, nspec
          do ie = 1, negrid
             do ixi = 1, nxi
                il = min(ixi, nxi+1-ixi)
                vpanud(:,ixi,ie,is) = sign(isgn, nlambda-ixi)*sqrt((1.0-al(il)*bmag)*e(ie,is))
             end do
          end do
       end do
       
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          is = is_idx(le_lo,ile)
          ig = ig_idx(le_lo,ile)
          gtmp(:,:,ile) = vpanud(ig,:,:,is) * aj0le(:,:,ile) * gle(:,:,ile)
       end do
       call integrate_moment (le_lo, gtmp, v0y0)    ! v0y0

       ! add part of ion-drag term
!        do ile = le_lo%llim_proc, le_lo%ulim_proc
!           ig = ig_idx(le_lo,ile)
!           ik = ik_idx(le_lo,ile)
!           it = it_idx(le_lo,ile)
!           is = is_idx(le_lo,ile)
!           gle(:,:,ile) = gle(:,:,ile) - z0le(:,:,ile) * v0y0(ig,it,ik,is)
!        end do

       do ile = le_lo%llim_proc, le_lo%ulim_proc
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)* z0le(:,:,ile) * v0y0(ile)
       end do

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    if (conservative) then
       do is = 1, nspec
          do ie = 1, negrid
             do ixi = 1, nxi
                il = min(ixi, nxi+1-ixi)
                isgn = 2 - il/ixi
                vpanud(:,ixi,ie,is) = vpdiff(:,isgn,il) * sqrt(e(ie,is)) * vnmult(1)*vnew_D(1,ie,is)/tunits(1)
             end do
          end do
       end do
    else
       isgn = 1
       do is = 1, nspec
          do ie = 1, negrid
             do ixi = 1, nxi
                il = min(ixi, nxi+1-ixi)
                vpanud(:,ixi,ie,is) = sign(isgn, nlambda-ixi)*sqrt((1.0-al(il)*bmag)*e(ie,is))*vnmult(1)*vnew_D(1,ie,is)/tunits(1)
             end do
          end do
       end do
    end if

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       gtmp(:,:,ile) = vpanud(ig,:,:,is) * tunits(ik) * aj0le(:,:,ile) &
            * gle(:,:,ile)
    end do

    call integrate_moment (le_lo, gtmp, v1y1)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

!     do ile = le_lo%llim_proc, le_lo%ulim_proc
!        it = it_idx(le_lo,ile)
!        ik = ik_idx(le_lo,ile)
!        ig = ig_idx(le_lo,ile)
!        is = is_idx(le_lo,ile)
!        gle(:,:,ile) = gle(:,:,ile) - s0le(:,:,ile) * v1y1(ig,it,ik,is)
!     end do

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*s0le(:,:,ile) * v1y1(ile)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       do ixi=1, nxi
          il = min(ixi, nxi+1-ixi)
          ! aj1vp2 = 2 * J1(arg)/arg * vperp^2
          gtmp(ixi,:negrid,ile) = vnmult(1) * vnew_D(ik,:negrid,is) * aj1le(ixi,:negrid,ile) &
               * e(:,is) * al(il) * gle(ixi,:negrid,ile)
       end do
    end do

    call integrate_moment (le_lo, gtmp, v2y2)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

!     do ile = le_lo%llim_proc, le_lo%ulim_proc
!        it = it_idx(le_lo,ile)
!        ik = ik_idx(le_lo,ile)
!        ig = ig_idx(le_lo,ile)
!        is = is_idx(le_lo,ile)
!        gle(:,:,ile) = gle(:,:,ile) - w0le(:,:,ile) * v2y2(ig,it,ik,is)
!     end do
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*w0le(:,:,ile) * v2y2(ile)
    end do

    deallocate (vpanud, v0y0, v1y1, v2y2)

# else

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v1y1(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v2y2(-ntgrid:ntgrid, ntheta0, naky, nspec))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (vns(naky,negrid,nspec))
    vns = vnmult(1)*vnew_D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*g(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v0y0, all)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g(:,isgn,iglo) - ieqzip(it,ik)*v0y0(:,it,ik,is) &
               * z0(:,isgn,iglo)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nud vpa J0 f0
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is)*sqrt(e(ie,is))*vpdiff(:,isgn,il) &
                  * aj0(:,iglo)*g1(:,isgn,iglo)
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * g1(:,isgn,iglo)
          end if
       end do
    end do

    call integrate_moment (gtmp, v1y1, all)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g1(:,isgn,iglo) - ieqzip(it,ik)*v1y1(:,it,ik,is) &
               * s0(:,isgn,iglo)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = nud vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v2y2, all)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g(:,isgn,iglo) = g1(:,isgn,iglo) - ieqzip(it,ik)*v2y2(:,it,ik,is) &
               * w0(:,isgn,iglo)
       end do
    end do

    deallocate (vns, v0y0, v1y1, v2y2)

# endif

# ifdef USE_LE_LAYOUT
    if (collision_model_switch == collision_model_lorentz) call scatter (g2le, gle, g)
# endif

  end subroutine conserve_lorentz

  subroutine conserve_diffuse (g, g1)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: e, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: ieqzip
# ifdef USE_LE_LAYOUT
    use le_grids, only: nlambda, ng2, forbid
    use gs2_layouts, only: le_lo, ig_idx
    use redistribute, only: scatter
    use theta_grid, only: bmag
    use run_parameters, only: tunits
# endif

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns

    integer :: ig, isgn, iglo, ik, ie, il, is, it, all = 1
# ifdef USE_LE_LAYOUT
    integer :: ile, ixi, nxi
    real, dimension (:,:,:,:), allocatable :: vpadelnu
    real, dimension (:,:), allocatable :: vpatmp
    complex, dimension (:), allocatable :: v0y0, v1y1, v2y2

    ! TMP FOR TESTING -- MAB
!    integer :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, tr
!    real :: t1tot=0., t2tot = 0., t3tot = 0., t4tot = 0., t5tot = 0., t6tot = 0., t7tot = 0., t8tot = 0., t9tot = 0., t10tot=0.
# else
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y1, v2y2    
# endif


# ifdef USE_LE_LAYOUT

    allocate (v0y0(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v1y1(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v2y2(le_lo%llim_proc:le_lo%ulim_alloc))

    nxi = max(2*nlambda-1,2*ng2)

    allocate (gtmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (vpadelnu(-ntgrid:ntgrid, nxi+1, negrid+1, nspec)) ; vpadelnu = 0.0
    allocate (vns(naky,negrid,nspec))
    allocate (vpatmp(-ntgrid:ntgrid,nxi)) ; vpatmp = 0.0

    vns = vnmult(2)*vnew_E

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0
    
    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t0, count_rate=tr)

     do ile = le_lo%llim_proc, le_lo%ulim_proc
        is = is_idx(le_lo,ile)
        ik = ik_idx(le_lo,ile)
       
        do ie=1, negrid
           do ixi = 1, nxi
              gtmp(ixi,ie,ile) = vns(ik,ie,is)*aj0le(ixi,ie,ile)*gle(ixi,ie,ile)
           end do
        end do
     end do

    call integrate_moment (le_lo, gtmp, v0y0)    ! v0y0
!    call integrate_moment (le_lo, gle*aj0le, v0y0, vns)    ! v0y0

    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t1)
!    t1tot = t1tot + real(t1-t0)/tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

!     do ile = le_lo%llim_proc, le_lo%ulim_proc
!        ig = ig_idx(le_lo,ile)
!        it = it_idx(le_lo,ile)
!        ik = ik_idx(le_lo,ile)
!        is = is_idx(le_lo,ile)
!        gle(:,:,ile) = gle(:,:,ile) - v0y0(ig,it,ik,is)*bz0le(:,:,ile)
!     end do
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*v0y0(ile)*bz0le(:,:,ile)
    end do

    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t3)
!    t3tot = t3tot + real(t3-t1)/tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    vns(1,:,:) = vnmult(2)*delvnew(1,:,:)/tunits(1)
    do ixi = 1, nxi
       il = min(ixi, nxi+1-ixi)
       isgn = 1
       do ig = -ntgrid, ntgrid
       	  if (.not. forbid(ig,il)) &
               vpatmp(ig,ixi) = sign(isgn,nlambda-ixi)*sqrt(max(0.0,(1.0-al(il)*bmag(ig))))
       end do
    end do

    do is = 1, nspec
       do ie = 1, negrid
!          do ixi = 1, nxi
!             il = min(ixi, nxi+1-ixi)
!             isgn = 1
!             do ig = -ntgrid, ntgrid
!                if (.not. forbid(ig,il)) &
!                     vpadelnu(ig,ixi,ie,is) = sign(isgn, nlambda-ixi) &
!                     * vns(1,ie,is) * sqrt(max(0.0,(1.0-al(il)*bmag(ig))*e(ie,is) ))
!             end do
!          end do
          vpadelnu(:,:nxi,ie,is) = vns(1,ie,is) * sqrt(vpatmp*e(ie,is))
       end do
    end do

!    call system_clock (count=t4)
!    t4tot = t4tot + real(t4-t3)/tr

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       gtmp (:,:,ile) = vpadelnu(ig,:,:,is) * tunits(ik) * aj0le(:,:,ile) * gle(:,:,ile)
    end do

    call integrate_moment (le_lo, gtmp, v1y1)    ! v1y1

    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t6)
!    t6tot = t6tot + real(t6-t4)/tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

!     do ile = le_lo%llim_proc, le_lo%ulim_proc
!        it = it_idx(le_lo,ile)
!        ik = ik_idx(le_lo,ile)
!        ig = ig_idx(le_lo,ile)
!        is = is_idx(le_lo,ile)
!        gle(:,:,ile) = gle(:,:,ile) - bs0le(:,:,ile) * v1y1(ig,it,ik,is)
!     end do
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*bs0le(:,:,ile) * v1y1(ile)
    end do

    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t7)
!    t7tot = t7tot + real(t7-t6)/tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    do ik = 1, naky
       vns(ik,:,:) = delvnew(ik,:,:)*e
    end do

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       do ie=1, negrid
          do ixi = 1, nxi
             il = min(ixi, nxi+1-ixi)
!             gtmp(ixi,ie,ile) = delvnew(ik,ie,is) &
!                  * e(ie,is)*al(il)*aj1le(ixi,ie,ile) * gle(ixi,ie,ile)
             gtmp(ixi,ie,ile) = vns(ik,ie,is) &
                  * al(il)*aj1le(ixi,ie,ile) * gle(ixi,ie,ile)
          end do
       end do
    end do

    call integrate_moment (le_lo, gtmp, v2y2)    ! v2y2

    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t9)
!    t9tot = t9tot + real(t9-t7)/tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

!     do ile = le_lo%llim_proc, le_lo%ulim_proc
!        it = it_idx(le_lo,ile)
!        ik = ik_idx(le_lo,ile)
!        ig = ig_idx(le_lo,ile)
!        is = is_idx(le_lo,ile)
!        gle(:,:,ile) = gle(:,:,ile) - bw0le(:,:,ile) * v2y2(ig,it,ik,is)
!     end do
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*bw0le(:,:,ile) * v2y2(ile)
    end do

    ! TMP FOR TESTING -- MAB
!    call system_clock (count=t10)
!    t10tot = t10tot + real(t10-t9)/tr

    if (collision_model_switch == collision_model_ediffuse) &
         call scatter (g2le, gle, g)

!    write (*,*) 'conserve_diffuse', t1tot, t3tot, t4tot, t6tot, t7tot, t9tot, t10tot
!    write (*,*)

    deallocate (vpadelnu, vns, vpatmp, v0y0, v1y1, v2y2)

# else

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v1y1(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v2y2(-ntgrid:ntgrid, ntheta0, naky, nspec))

    allocate (vns(naky,negrid,nspec))
    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    vns = vnmult(2)*delvnew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0 f0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * g(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v0y0, all)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g(:,isgn,iglo) - ieqzip(it,ik)*v0y0(:,it,ik,is) &
               * bz0(:,isgn,iglo)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = (nus-nud) vpa J0 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*vpa(:,isgn,iglo)*aj0(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v1y1, all)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g1(:,isgn,iglo) - ieqzip(it,ik)*v1y1(:,it,ik,is) &
               * bs0(:,isgn,iglo)
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = (nus-nud) vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v2y2, all)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g(:,isgn,iglo) = g1(:,isgn,iglo) - ieqzip(it,ik)*v2y2(:,it,ik,is) &
               * bw0(:,isgn,iglo)
       end do
    end do

    deallocate (vns, v0y0, v1y1, v2y2)

# endif

    deallocate (gtmp)

  end subroutine conserve_diffuse

  subroutine g_adjust (g, phi, bpar, facphi, facbpar)
    use species, only: spec
    use theta_grid, only: ntgrid
    use le_grids, only: anon
    use dist_fn_arrays, only: vperp2, aj0, aj1
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar

    integer :: iglo, ig, ik, it, ie, is
    complex :: adj

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          adj = anon(ie,is)*2.0*vperp2(ig,iglo)*aj1(ig,iglo) &
                  *bpar(ig,it,ik)*facbpar &
               + spec(is)%z*anon(ie,is)*phi(ig,it,ik)*aj0(ig,iglo) &
                  /spec(is)%temp*facphi
          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do
  end subroutine g_adjust

  subroutine solfp_krook (g, g1)
    use species, only: spec, electron_species
    use theta_grid, only: ntgrid
    use le_grids, only: integrate, lintegrate, geint2g, gint2g
    use gs2_layouts, only: g_lo, gint_lo, geint_lo
    use gs2_layouts, only: ik_idx, il_idx, it_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1

    complex, dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc) &
         :: g0eint, g1eint
    complex, dimension (-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc) &
         :: g1int, g2int
    integer :: iglo, igint, igeint, ig, ik, it, il, is

    call prof_entering ("solfp_krook", "collisions")

    if (conserve_momentum) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
             g1(:,:,iglo) = 0.0
          else
             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
          end if
       end do
       call lintegrate (g1, g1int)
    end if

    if (conserve_number) call integrate (g, g0eint)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          g1(ig,1,iglo) = g(ig,1,iglo)/(1.0 + vnewfe(iglo))
          g1(ig,2,iglo) = g(ig,2,iglo)/(1.0 + vnewfe(iglo))
       end do
    end do

    if (conserve_number) then
       call integrate (g1, g1eint)
       do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
          g0eint(:,igeint) = (g0eint(:,igeint) - g1eint(:,igeint)) &
               *aintnorm(:,igeint)
       end do
       call geint2g (g0eint, g)
       g = g1 + g
    else
       g = g1
    end if

    if (conserve_momentum) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
             g1(:,:,iglo) = 0.0
          else
             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
          end if
       end do
       call lintegrate (g1, g2int)

       do igint = gint_lo%llim_proc, gint_lo%ulim_proc
          g1int(:,igint) = (g1int(:,igint) - g2int(:,igint))/g3int(:,igint)
       end do
       call gint2g (g1int, g1)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
          else if (spec(is)%type == electron_species) then
          else
             g(:,:,iglo) = g(:,:,iglo) + sq(:,il,:)*g1(:,:,iglo)
          end if
       end do
    end if

    call prof_leaving ("solfp_krook", "collisions")
  end subroutine solfp_krook

  subroutine solfp_lorentz (g, gc, gh, diagnostics, init)
    use species, only: spec, electron_species
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: nlambda, jend, lintegrate, ng2, al
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: ig_idx, ik_idx, il_idx, is_idx, it_idx, ie_idx
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use run_parameters, only: ieqzip
# ifdef USE_LE_LAYOUT
    use gs2_layouts, only: le_lo
    use le_grids, only: negrid, jend
# else
    use gs2_layouts, only: lz_lo
# endif

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gc, gh
    integer, optional, intent (in) :: diagnostics
    logical, optional, intent (in) :: init

# ifdef USE_LE_LAYOUT
    complex, dimension (:), allocatable :: gle0
    integer :: ile
# else
    complex, dimension (:,:), allocatable :: glz, glzc
    complex, dimension (:), allocatable :: glz0
    integer :: ilz
# endif
    complex, dimension (max(2*nlambda,2*ng2+1)) :: delta
    complex :: fac, gwfb
    integer :: iglo, ig, ik, il, is, je, it, ie, isgn, ixi, nxi

    nxi = max(2*nlambda-1, 2*ng2)
    call prof_entering ("solfp_lorentz", "collisions")

# ifdef USE_LE_LAYOUT

    if (collision_model_switch == collision_model_lorentz &
         .or. collision_model_switch == collision_model_lorentz_test) &
         call gather (g2le, g, gle)

    allocate (gle0(max(2*nlambda,2*ng2+1)))

    if (heating .and. present(diagnostics)) then
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ig = ig_idx(le_lo,ile)

          je = 2*jend(ig)          
          if (je == 0) then
             je = 2*ng2 
          end if

! when il=je-1 below, and we have trapped particles, gle is evaluated at gle(2*jend(ig),ie,ile).
! this seems like a bug, since there are only 2*jend(ig)-1 grid points and
! the value gle(2*jend(ig),ie,ile) corresponds to the value of g at xi = 0...this
! doesn't make any sense...MAB

          do ie = 1, negrid
             do il = 1, je-1
                fac = gle(il+1,ie,ile)-gle(il,ie,ile)
                glec(il,ie,ile) = conjg(fac)*fac*d1le(il,ie,ile)  ! d1le accounts for hC(h) entropy
             end do
          end do
       end do
       call scatter (g2le, glec, gc)

       if (hyper_colls) then
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             ig = ig_idx(le_lo,ile)

             je = 2*jend(ig)          
             if (je == 0) then
                je = 2*ng2 
             end if
             do ie = 1, negrid
                do il = 1, je-1
                   fac = gle(il+1,ie,ile)-gle(il,ie,ile)
                   glec(il,ie,ile) = conjg(fac)*fac*h1le(il,ie,ile)  ! h1le accounts for hH(h) entropy
                end do
             end do
          end do
          call scatter (g2le, glec, gh)
       end if
    end if

    ! solve for glz row by row
    do ile = le_lo%llim_proc, le_lo%ulim_proc

       ig = ig_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       is = is_idx(le_lo,ile)

       if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) cycle
       if (ieqzip(it_idx(le_lo,ile),ik_idx(le_lo,ile))==0) cycle

       je = jend(ig)

       if (je > ng2+1) then
          je = 2*je
       else
          je = 2*ng2+1
       end if

       do ie = 1, negrid
          ! deal with special case of wfb
          if (jend(ig) == ng2+1) then
             ! if wfb, remove vpa = 0 point (which has wgt of zero)
             gle0(:ng2) = gle(:ng2,ie,ile)
             gle0(ng2+1:je-1) = gle(ng2+2:je,ie,ile)
             ! save gwfb for reinsertion later
             gwfb = gle(ng2+1,ie,ile)
          else
             gle0 = gle(:,ie,ile)
          end if

          gle(:je-1,ie,ile) = gle0(:je-1)
          gle(je:,ie,ile) = 0.0

          ! right and left sweeps for tridiagonal solve:
          
          delta(1) = gle(1,ie,ile)
          do il = 1, je-1
             delta(il+1) = gle(il+1,ie,ile) - qle(il+1,ie,ile)*delta(il)
          end do
       
          gle(je,ie,ile) = delta(je)*betaale(je,ie,ile)
          do il = je-1, 1, -1
             gle(il,ie,ile) = (delta(il) - c1le(il,ie,ile)*gle(il+1,ie,ile))*betaale(il,ie,ile)
          end do
!       ! interpolate to obtain glz(vpa = 0) point for wfb
!       ! and insert this point into glz
          ! interpolation described above mysteriously causing numerical instability
          ! stabilized by using old (pre-collision) value of g for wfb
          if (jend(ig) == ng2+1) then
             gle0(ng2+2:je) = gle(ng2+1:je-1,ie,ile)
!          gle(ng2+1,ie,ile) = 0.5*(gle(ng2,ie,ile)+gle(ng2+1,ie,ile))
             gle(ng2+1,ie,ile) = gwfb
             gle(ng2+2:je,ie,ile) = gle0(ng2+2:je)
          end if

          if (jend(ig) /= 0) gle(2*jend(ig),ie,ile) = gle(jend(ig),ie,ile)

       end do
    end do

    if (present(init)) then
       call scatter (g2le, gle, g)
!    else if (collision_model_switch == collision_model_lorentz .and. .not.conserve_moments) then
    else if (.not.conserve_moments) then
       call scatter (g2le, gle, g)
    end if

    deallocate (gle0)
    
# else

    allocate (glz(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (glz0(max(2*nlambda,2*ng2+1)))
    glz = 0.0 ; glz0 = 0.0
    if (heating) then
       allocate (glzc(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       glzc = 0.0
    end if

!    call check_g ('beg', g)

    call gather (lorentz_map, g, glz)


    if (heating .and. present(diagnostics)) then
       do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
          ig = ig_idx(lz_lo,ilz)

          je = 2*jend(ig)          
          if (je == 0) then
             je = 2*ng2 
          end if

! when il=je-1 below, and we have trapped particles, glz is evaluated at glz(2*jend(ig),ilz).
! this seems like a bug, since there are only 2*jend(ig)-1 grid points and
! the value glz(2*jend(ig),ilz) corresponds to the value of g at xi = 0...this
! doesn't make any sense...MAB

          do il = 1, je-1
             fac = glz(il+1,ilz)-glz(il,ilz)
             glzc(il,ilz) = conjg(fac)*fac*d1(il,ilz)  ! d1 accounts for hC(h) entropy
          end do
       end do
       call scatter (lorentz_map, glzc, gc)

       if (hyper_colls) then
          do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
             ig = ig_idx(lz_lo,ilz)
             
             je = 2*jend(ig)          
             if (je == 0) then
                je = 2*ng2 
             end if
             
             do il = 1, je-1
                fac = glz(il+1,ilz)-glz(il,ilz)
                glzc(il,ilz) = conjg(fac)*fac*h1(il,ilz)  ! h1 accounts for hH(h) entropy
             end do
          end do
          call scatter (lorentz_map, glzc, gh)
       end if
    end if

!    call check_glz ('beg', glz)
    ! solve for glz row by row
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ie = ie_idx(lz_lo,ilz)
       ig = ig_idx(lz_lo,ilz)
       ik = ik_idx(lz_lo,ilz)
       is = is_idx(lz_lo,ilz)

       if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) cycle
       if (ieqzip(it_idx(lz_lo,ilz),ik_idx(lz_lo,ilz))==0) cycle

       je = jend(ig)

       if (je > ng2+1) then
          je = 2*je
       else
          je = 2*ng2+1
       end if

       ! deal with special case of wfb
       if (jend(ig) == ng2+1) then
          ! if wfb, remove vpa = 0 point (which has wgt of zero)
          glz0(:ng2) = glz(:ng2,ilz)
          glz0(ng2+1:je-1) = glz(ng2+2:je,ilz)
          ! save gwfb for reinsertion later
          gwfb = glz(ng2+1,ilz)
       else
          glz0 = glz(:,ilz)
       end if

       glz(:je-1,ilz) = glz0(:je-1)
       glz(je:,ilz) = 0.0

       ! right and left sweeps for tridiagonal solve:

       delta(1) = glz(1,ilz)
       do il = 1, je-1
          delta(il+1) = glz(il+1,ilz) - ql(il+1,ilz)*delta(il)
       end do
       
       glz(je,ilz) = delta(je)*betaa(je,ilz)
       do il = je-1, 1, -1
          glz(il,ilz) = (delta(il) - c1(il,ilz)*glz(il+1,ilz))*betaa(il,ilz)
       end do

!       ! interpolate to obtain glz(vpa = 0) point for wfb
!       ! and insert this point into glz
       ! interpolation described above mysteriously causing numerical instability
       ! stabilized by using old (pre-collision) value of g for wfb
       if (jend(ig) == ng2+1) then
          glz0(ng2+2:je) = glz(ng2+1:je-1,ilz)
!          glz(ng2+1,ilz) = 0.5*(glz(ng2,ilz)+glz(ng2+1,ilz))
          glz(ng2+1,ilz) = gwfb
          glz(ng2+2:je,ilz) = glz0(ng2+2:je)
       end if

!
! bug fixed 4.14.03
! was only a problem with collisions but with no trapped particles
! because of an array index out of bounds 
! Overall, this was a rare bug.
!

!       if (jend(ig) /= 0) glz(je,ilz) = glz(jend(ig),ilz)
       if (jend(ig) /= 0) glz(2*jend(ig),ilz) = glz(jend(ig),ilz)

    end do

!    call check_glz ('end', glz)

    call scatter (lorentz_map, glz, g)

!    call check_g ('end', g)

    deallocate (glz, glz0)
    if (heating) deallocate (glzc)

# endif

    call prof_leaving ("solfp_lorentz", "collisions")
  end subroutine solfp_lorentz

  subroutine solfp_ediffuse (g, init)
    use species, only: spec, nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use le_grids, only: negrid, integrate_moment, forbid, ng2
    use gs2_layouts, only: ig_idx, il_idx, is_idx, e_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use run_parameters, only: ieqzip
# ifdef USE_LE_LAYOUT
    use le_grids, only: nlambda, jend
    use gs2_layouts, only: le_lo
# endif

    ! TEMP FOR TESTING -- MAB
    use gs2_layouts, only: ik_idx, ie_idx, il_idx, it_idx
    use gs2_time, only: code_dt
    use le_grids, only: e

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    logical, intent (out), optional :: init

    integer :: ie, is, iglo, ig, isgn, it, ik, il
    complex, dimension (negrid) :: delta
# ifdef USE_LE_LAYOUT
    integer :: ile, nxi, ixi
# else
    complex, dimension (:,:), allocatable :: ged
    integer :: ielo
# endif

# ifdef USE_LE_LAYOUT

    call gather (g2le, g, gle)

    nxi = max(2*nlambda-1, 2*ng2)

    ! solve for ged row by row
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       if (spec(is)%vnewk < 2.0*epsilon(0.0)) cycle
       if (ieqzip(it_idx(le_lo,ile),ik_idx(le_lo,ile))==0) cycle
       ig = ig_idx(le_lo,ile)
       do ixi = 1, nxi
          il = min(ixi, nxi+1-ixi)
          if (forbid(ig,il)) cycle

          delta(1) = gle(ixi,1,ile)
          do ie = 1, negrid-1
             delta(ie+1) = gle(ixi,ie+1,ile) - eqle(ixi,ie+1,ile)*delta(ie)
          end do

          gle(ixi,negrid+1,ile) = 0.0
          do ie = negrid, 1, -1
             gle(ixi,ie,ile) = (delta(ie) - ec1le(ixi,ie,ile)*gle(ixi,ie+1,ile)) &
                  * ebetaale(ixi,ie,ile)
          end do
       end do
    end do

    if (present(init)) then
       call scatter (g2le, gle, g)
    else if (collision_model_switch == collision_model_ediffuse &
         .and. .not. conserve_moments) then
       call scatter (g2le, gle, g)
    end if

# else

    allocate (ged(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc))
    ged = 0.0

    call gather (ediffuse_map, g, ged)

    ! solve for ged row by row
    do ielo = e_lo%llim_proc, e_lo%ulim_proc
       is = is_idx(e_lo,ielo)
       ig = ig_idx(e_lo,ielo)
       il = il_idx(e_lo,ielo)

       if (spec(is)%vnewk < 2.0*epsilon(0.0) .or. forbid(ig,il)) cycle
       if (ieqzip(it_idx(e_lo,ielo),ik_idx(e_lo,ielo))==0) cycle

       delta(1) = ged(1,ielo)
       do ie = 1, negrid-1
          delta(ie+1) = ged(ie+1,ielo) - eql(ie+1,ielo)*delta(ie)
       end do
       
       ged(negrid+1,ielo) = 0.0
       do ie = negrid, 1, -1
          ged(ie,ielo) = (delta(ie) - ec1(ie,ielo)*ged(ie+1,ielo))*ebetaa(ie,ielo)
       end do

    end do

    call scatter (ediffuse_map, ged, g)

    deallocate (ged)

# endif

  end subroutine solfp_ediffuse

  subroutine check_g (str, g)

    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    
    character (3) :: str 
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    integer :: ig, iglo

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, 0
          if (g(ig,1,iglo) /= g(-ig,2,iglo)) then
             if (abs((g(ig,1,iglo)-g(-ig,2,iglo))/g(ig,1,iglo)) > 1.e-8) then
                write(*,*) str,ig,iglo,g(ig,1,iglo)-g(-ig,2,iglo)
             endif
          endif
       end do
    end do

  end subroutine check_g

  subroutine check_glz (str, glz)

    use gs2_layouts, only: lz_lo, ig_idx, ie_idx
    use le_grids, only: jend
    use theta_grid, only: ntgrid
    
    character (3) :: str 
    complex, dimension (:,lz_lo%llim_proc:), intent (in) :: glz
    integer :: il, ilz, ig, ie, je, ilzp, ilp

    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx (lz_lo, ilz)
       ie = ie_idx (lz_lo, ilz)
       je = jend(ig)
       ilzp = -ig + ntgrid + (ie-1)*(2*ntgrid+1)
       do il = 1, je
          ilp = 2*je+1-il
          if (glz(il,ilz) /= glz(ilp,ilzp)) then
             if (abs((glz(il,ilz)-glz(ilp,ilzp))/glz(il,ilz)) > 1.e-8) then
                write(*,*) str,il,ilz,glz(il,ilz),glz(ilp,ilzp)
             endif
          endif
       end do
    end do

  end subroutine check_glz

  subroutine init_vpdiff

    use le_grids, only: al, wl, nlambda, jend, ng2
    use theta_grid, only: ntgrid, bmag

    implicit none

    integer :: il, ig, je, te
    real :: slb0, slb1, slb2, slbl, slbr
    real, dimension (nlambda) :: xi
    real, dimension (2,nlambda) :: vptmp

    if (.not. allocated(vpdiff) .and. conservative) then

       allocate (vpdiff(-ntgrid:ntgrid,2,nlambda))
       
       vpdiff = 0.0
       
       do ig = -ntgrid, ntgrid
          je = jend(ig)
          if (je <= ng2+1) then
             te = ng2
          else
             te = je
          end if
          do il = 2, te-1
             slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
             slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1))) 
             
             slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
             slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
             
             vpdiff(ig,1,il) = (slbl**2 - slbr**2)/wl(ig,il)
          end do

          ! boundary at xi = 1
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(1))) 
          slb2 = sqrt(abs(1.0 - bmag(ig)*al(2))) 
          slbr = 0.5*(slb1 + slb2)
          vpdiff(ig,1,1) = (1.0 - slbr**2)/wl(ig,1)
          
          ! boundary at xi = 0
          il = te
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
          if (te == ng2) then
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
             slb2 = -slb1
          else
             slb1 = 0.0
             slb2 = -slb0
          end if
          
          slbl = (slb1 + slb0)*0.5
          slbr = (slb1 + slb2)*0.5
          vpdiff(ig,1,il) = (slbl**2 - slbr**2)/wl(ig,il)
          
          vpdiff(ig,2,:) = -vpdiff(ig,1,:)
          
       end do
       
    end if
    
  end subroutine init_vpdiff

  subroutine check_g2le

    use file_utils, only: error_unit
    use mp, only: finish_mp, iproc, nproc, proc0
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, negrid
    use gs2_layouts, only: g_lo, le_lo
    use gs2_layouts, only: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use redistribute, only: gather, scatter, report_map_property
    integer :: iglo, ile, ig, ik, it, il, ie, is, ierr
    integer :: ip
    complex, dimension (:,:,:), allocatable :: gtmp, letmp

    if (proc0) then
       ierr = error_unit()
    else
       ierr = 6
    end if

    ! report the map property
    if (proc0) print *, '=== g2le map property ==='
    call report_map_property (g2le)

    allocate (gtmp(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (letmp(nlambda*2+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))

    ! gather check
    gtmp = 0.0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       gtmp(0,1,iglo) = rule(ik,it,il,ie,is)
    end do
    call gather (g2le, gtmp, letmp)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ig = ig_idx(le_lo,ile)
! MAB>
! following line changed to account for extra elements in gle layout
! (i.e. il=2*nlambda+1 and ie=negrid+1)
!       if (ig /= 0 .and. all(letmp(:,:,ile)==0.0)) cycle
       if (ig /= 0 .and. all(letmp(:2*nlambda,:negrid,ile)==0.0)) cycle
! <MAB
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       is = is_idx(le_lo,ile)
       do il=1, nlambda
          do ie=1, negrid
             if (int(real(letmp(il,ie,ile))) /= rule(ik,it,il,ie,is)) &
                  write (ierr,'(a,8i6)') 'ERROR: gather by g2le broken!', iproc
          end do
       end do
    end do
    if (proc0) write (ierr,'(a)') 'g2le gather check done'

    ! scatter check
    gtmp = 0.0
    call scatter (g2le, letmp, gtmp)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       if (gtmp(0,1,iglo) /= rule(ik,it,il,ie,is)) &
            write (ierr,'(a,i6)') 'ERROR: scatter by g2le broken!', iproc
    end do
    if (proc0) write (ierr,'(a)') 'g2le scatter check done'

    deallocate (gtmp,letmp)

!!$    call finish_mp
!!$    stop

  contains

    function rule (ik, it, il, ie, is)
      integer, intent (in) :: ik, it, il, ie, is
      integer :: rule
      rule = ie + ik  ! make whatever you want
    end function rule

  end subroutine check_g2le

  subroutine reset_init
!
! forces recalculation of coefficients in collision operator
! when timestep changes.
!    
    initialized = .false.  

  end subroutine reset_init

  subroutine finish_collisions

    use dist_fn_arrays, only: c_rate

    implicit none

    vnmult = 0.0
    accelerated_x = .false. ; accelerated_v = .false.
    ediffinit = .false. ; lzinit = .false. ; leinit = .false.
    call reset_init

    if (allocated(c_rate)) deallocate (c_rate)
    if (allocated(z0)) deallocate (z0, w0, s0)
    if (allocated(bz0)) deallocate (bz0, bw0, bs0)
    if (allocated(vnew)) deallocate (vnew, vnew_s, vnew_D, vnew_E, delvnew)
    if (allocated(vnewh)) deallocate (vnewh)
    if (allocated(sq)) deallocate (sq)
    if (allocated(g3int)) deallocate (g3int)
    if (allocated(aintnorm)) deallocate (aintnorm)
    if (allocated(vnewfe)) deallocate (vnewfe)

# ifdef USE_LE_LAYOUT
    if (allocated(c1le)) then
       deallocate (c1le, betaale, qle)
       if (heating) deallocate (d1le, h1le)
    end if
    if (allocated(ec1le)) deallocate (ec1le, ebetaale, eqle)
    if (allocated(gle)) deallocate (gle)
    if (allocated(glec)) deallocate (glec)
# else
    if (allocated(c1)) then
       deallocate (c1, betaa, ql)
       if (heating) deallocate (d1, h1)
    end if
    if (allocated(ec1)) deallocate (ec1, ebetaa, eql)
# endif
    if (allocated(vpdiff)) deallocate (vpdiff)
    if (allocated(dtot)) deallocate (dtot, fdf, fdb)

  end subroutine finish_collisions

end module collisions

! OLD MOMENTUM CONSERVATION CODING
!    if (conserve_momentum) then
!       allocate (g1int(-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc))
!       do iglo = g_lo%llim_proc, g_lo%ulim_proc
!          ik = ik_idx(g_lo,iglo)
!          il = il_idx(g_lo,iglo)
!          is = is_idx(g_lo,iglo)
!          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
!             g1(:,:,iglo) = 0.0
!          else
!             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
!          end if
!       end do
!       call lintegrate (g1, g1int)
!    end if

! ...

!    if (conserve_momentum) then
!       do iglo = g_lo%llim_proc, g_lo%ulim_proc
!          ik = ik_idx(g_lo,iglo)
!          il = il_idx(g_lo,iglo)
!          is = is_idx(g_lo,iglo)
!          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
!             g1(:,:,iglo) = 0.0
!          else
!             g1(:,:,iglo) = g(:,:,iglo)*sq(:,il,:)
!          end if
!       end do
!       allocate (g2int(-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_alloc))
!       call lintegrate (g1, g2int)
!
!       do igint = gint_lo%llim_proc, gint_lo%ulim_proc
!          g1int(:,igint) = (g1int(:,igint) - g2int(:,igint))/g3int(:,igint)
!       end do
!       call gint2g (g1int, g1)
!       deallocate (g1int, g2int)
!
!       do iglo = g_lo%llim_proc, g_lo%ulim_proc
!          ik = ik_idx(g_lo,iglo)
!          il = il_idx(g_lo,iglo)
!          is = is_idx(g_lo,iglo)
!          if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) then
!          else if (spec(is)%type == electron_species) then
!          else
!             g(:,:,iglo) = g(:,:,iglo) + sq(:,il,:)*g1(:,:,iglo)
!          end if
!       end do
!    end if

