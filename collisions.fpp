# include "define.inc"

!>  Routines for implementing the model collision operator
!!defined by Barnes, Abel et al. 2009. The collision operator causes 
!! physically motivated smoothing of structure in velocity space which 
!!is necessary to prevent buildup of structure at fine scales in 
!!velocity space, while conserving energy and momentum.

module collisions
  
  use redistribute, only: redist_type

  implicit none

  private 

  public :: init_collisions, finish_collisions
  public :: read_parameters, wnml_collisions, check_collisions
  public :: solfp1
  public :: reset_init
  public :: dtot, fdf, fdb, vnmult, vnfac
  public :: ncheck, vnslow, vary_vnew
  public :: etol, ewindow, etola, ewindowa
  public :: init_lorentz, init_ediffuse
  public :: init_lorentz_conserve, init_diffuse_conserve
  public :: init_lorentz_error, collision_model_switch
  public :: colls, hyper_colls, heating, adjust
  public :: use_le_layout, collision_model_ediffuse
  public :: collision_model_lorentz, collision_model_none
  public :: collision_model_lorentz_test, collision_model_full
  
  interface solfp1
     module procedure solfp1_le_layout
     module procedure solfp1_standard_layout
  end interface

  interface solfp_lorentz
     module procedure solfp_lorentz_le_layout
     module procedure solfp_lorentz_standard_layout
  end interface

  interface conserve_lorentz
     module procedure conserve_lorentz_le_layout
     module procedure conserve_lorentz_standard_layout
  end interface

  interface conserve_diffuse
     module procedure conserve_diffuse_le_layout
     module procedure conserve_diffuse_standard_layout
  end interface 

  ! knobs
  logical :: use_le_layout = .false.
  logical :: const_v, conserve_moments
  logical :: conservative, resistivity
  integer :: collision_model_switch
  integer :: lorentz_switch, ediff_switch
  logical :: adjust
  logical :: heating
  logical :: hyper_colls
  logical :: ei_coll_only
  logical :: test
  logical :: special_wfb_lorentz
  integer, parameter :: collision_model_lorentz = 1      ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_lorentz_test = 5 ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_full = 6
  integer, parameter :: collision_model_ediffuse = 7

  integer, parameter :: lorentz_scheme_default = 1
  integer, parameter :: lorentz_scheme_old = 2

  integer, parameter :: ediff_scheme_default = 1
  integer, parameter :: ediff_scheme_old = 2

  real, dimension (2) :: vnmult = 0.0
  integer :: ncheck
  logical :: vary_vnew
  real :: vnfac, vnslow
  real :: etol, ewindow, etola, ewindowa

  real, dimension (:,:,:), allocatable :: dtot
  ! (-ntgrid:ntgrid,nlambda,max(ng2,nlambda-ng2)) lagrange coefficients for derivative error estimate

  real, dimension (:,:), allocatable :: fdf, fdb
  ! (-ntgrid:ntgrid,nlambda) finite difference coefficients for derivative error estimate

  real, dimension (:,:,:), allocatable :: vnew, vnew_s, vnew_D, vnew_E, delvnew
  ! (naky,negrid,nspec) replicated

  real, dimension (:,:,:), allocatable :: vpdiff
  ! (-ntgrid:ntgrid,2,nlambda) replicated

  complex, dimension (:,:,:), allocatable :: glec
  ! (2*nlambda+1,negrid+1,le_lo%llim_proc:le_lo%ulim_alloc)

  ! only for hyper-diffusive collisions
  real, dimension (:,:,:,:), allocatable :: vnewh
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec) replicated

  ! only for momentum conservation due to Lorentz operator (8.06)
  real, dimension(:,:,:), allocatable :: s0, w0, z0
  ! The following (between the start and finish comments) are only used for LE layout
  ! start
  real, dimension (:,:,:), allocatable :: s0le, w0le, z0le
  real, dimension (:,:,:), allocatable :: aj0le, aj1le
  ! finish

  ! needed for momentum and energy conservation due to energy diffusion (3.08)
  real, dimension(:,:,:), allocatable :: bs0, bw0, bz0

  ! The following (between the start and finish comments) are only used for LE layout
  ! start
  real, dimension (:,:), allocatable :: vpatmp
  real, dimension (:,:,:), allocatable :: bs0le, bw0le, bz0le
  ! finish

  real :: cfac

  ! The following (between the start and finish comments) are only used for LE layout
  ! start
  ! only for lorentz
  real, dimension (:,:,:), allocatable :: c1le, betaale, qle, d1le, h1le
  ! only for energy diffusion
  real, dimension (:,:,:), allocatable :: ec1le, ebetaale, eqle
  ! finish
  ! The following (between the start and finish comments) are only used for none LE layout
  ! start
  ! only for lorentz
  real, dimension (:,:), allocatable :: c1, betaa, ql, d1, h1
  ! only for energy diffusion
  real, dimension (:,:), allocatable :: ec1, ebetaa, eql
  ! finish

  logical :: drag = .false.
  logical :: heating_flag = .false.
  logical :: colls = .true.

  logical :: hypermult
  logical :: initialized = .false.
  logical :: ediffinit = .false.
  logical :: lzinit = .false., leinit = .false.
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.
  logical :: exist

contains

  subroutine check_collisions(report_unit)
    implicit none
    integer, intent(in) :: report_unit
    select case (collision_model_switch)
    case (collision_model_lorentz,collision_model_lorentz_test)
       write (report_unit, fmt="('A Lorentz collision operator has been selected.')")
       if (cfac > 0) write (report_unit, fmt="('This has both terms of the Lorentz collision operator: cfac=',e12.4)") cfac
       if (cfac == 0) write (report_unit, fmt="('This is only a partial Lorentz collision operator (cfac=0.0)')")
       if (const_v) write (report_unit, fmt="('This is an energy independent Lorentz collision operator (const_v=true)')")  
!          if (hypercoll) call init_hyper_lorentz
    case (collision_model_full)
       write (report_unit, fmt="('Full GS2 collision operator has been selected.')")
    end select
  end subroutine check_collisions

  subroutine wnml_collisions(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "collisions_knobs"
    select case (collision_model_switch)
    case (collision_model_lorentz)
       write (unit, fmt="(' collision_model = ',a)") '"lorentz"'
       if (hypermult) write (unit, fmt="(' hypermult = ',L1)") hypermult
    case (collision_model_lorentz_test)
       write (unit, fmt="(' collision_model = ',a)") '"lorentz-test"'
    case (collision_model_none)
       write (unit, fmt="(' collision_model = ',a)") '"collisionless"'
    end select
    write (unit, fmt="(' cfac = ',f6.3)") cfac
    write (unit, fmt="(' heating = ',L1)") heating
    write (unit, fmt="(' /')")
  end subroutine wnml_collisions

  subroutine init_collisions
    use species, only: init_species, nspec, spec
    use theta_grid, only: init_theta_grid
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
    call init_dist_fn_layouts (naky, ntheta0, nlambda, negrid, nspec)
    call read_parameters

    call init_arrays

  end subroutine init_collisions

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    use run_parameters, only: beta, fapar
    use species, only: nspec
    implicit none
    type (text_option), dimension (6), parameter :: modelopts = &
         (/ text_option('default', collision_model_full), &
            text_option('lorentz', collision_model_lorentz), &
            text_option('ediffuse', collision_model_ediffuse), &
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
         use_le_layout, special_wfb_lorentz
    integer :: ierr, in_file

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
       conservative = .true.
       resistivity = .true.
       const_v = .false.
       heating = .false.
       test = .false.
       ei_coll_only = .false.
       special_wfb_lorentz=.true.
       !USE_LE_LAYOUT default set at module level
       !Note: If we want to move default here then we need to note
       !      that if LOWFLOW then we want default to be .true.
       in_file = input_unit_exist ("collisions_knobs", exist)
!       if (exist) read (unit=input_unit("collisions_knobs"), nml=collisions_knobs)
       if (exist) read (unit=in_file, nml=collisions_knobs)

       ierr = error_unit()
       call get_option_value &
            (collision_model, modelopts, collision_model_switch, &
            ierr, "collision_model in collisions_knobs",.true.)

       call get_option_value &
            (lorentz_scheme, schemeopts, lorentz_switch, &
            ierr, "lorentz_scheme in collisions_knobs",.true.)

       call get_option_value &
            (ediff_scheme, eschemeopts, ediff_switch, &
            ierr, "ediff_scheme in collisions_knobs",.true.)
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
    call broadcast (use_le_layout)
    call broadcast (special_wfb_lorentz)

    drag = resistivity .and. (beta > epsilon(0.0)) .and. (nspec > 1) .and. (fapar.gt.0)
  end subroutine read_parameters

  subroutine init_arrays
    use species, only: nspec
    use le_grids, only: negrid, init_map
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: c_rate

    implicit none

    real, dimension (negrid) :: hee
    logical :: use_lz_layout = .false.
    logical :: use_e_layout = .false.
! lowflow terms include higher-order corrections to GK equation
! such as parallel nonlinearity that require derivatives in v-space.
! most efficient way to take these derivatives is to go from g_lo to le_lo,
! i.e., bring all energies and lambdas onto each processor
!<DD>Note the user can still disable use_le_layout in the input file
!    this just changes the default.
# ifdef LOWFLOW
    use_le_layout = .true.
# endif

    !<DD>This (potentially large) array should only be needed if 
    !write_hrate (gs2_diagnostics) is true.
    if (.not. allocated(c_rate)) then
       allocate (c_rate(-ntgrid:ntgrid, ntheta0, naky, nspec, 3))
       c_rate = 0.
    end if

    !<DD>Should probably put this first
    if (collision_model_switch == collision_model_none) then
       colls = .false.
       return
    end if

    if( .not. use_le_layout ) then
       select case (collision_model_switch)
       case (collision_model_full)
          use_lz_layout = .true.
          use_e_layout = .true.
       case (collision_model_lorentz,collision_model_lorentz_test)
          use_lz_layout = .true.
       case (collision_model_ediffuse)
          use_e_layout = .true.
       end select
    end if

    call init_map (use_lz_layout, use_e_layout, use_le_layout, test)

    call init_vnew (hee)
    if (all(abs(vnew(:,1,:)) <= 2.0*epsilon(0.0))) then
       collision_model_switch = collision_model_none
       colls = .false.
       return
    end if

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
    end select

  end subroutine init_arrays

  subroutine init_lorentz_conserve

! Precompute three quantities needed for momentum and energy conservation:
! z0, w0, s0 (z0le, w0le, s0le if use_le_layout chosen in the input file defined)
    
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec, electron_species
    use kt_grids, only: kperp2, naky, ntheta0
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy, al, integrate_moment, negrid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: tunits
    use le_grids, only: g2le, nxi
    use gs2_layouts, only: le_lo
    use redistribute, only: gather, scatter

    implicit none
    
    complex, dimension (1,1,1) :: dum1 = 0., dum2 = 0.
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: duinv, dtmp
    real, dimension (:,:,:,:), allocatable :: vns
    integer :: ie, il, ik, is, isgn, iglo,  it
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
    complex, dimension(:,:,:), allocatable :: s0tmp, w0tmp, z0tmp
    logical :: all = .true.

    if(use_le_layout) then
       allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end if
    
! TO DO: 
! tunits not included anywhere yet

    allocate(s0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(w0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(z0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (dtmp(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (vns(naky,negrid,nspec,3))
    
    duinv=0.0 ; dtmp=0.0; !This initialisation is needed in case kwork_filter is true anywhere

    vns(:,:,:,1) = vnmult(1)*vnew_D
    vns(:,:,:,2) = vnmult(1)*vnew_s
    vns(:,:,:,3) = 0.0

    if (drag) then
       do is = 1, nspec
          if (spec(is)%type /= electron_species) cycle
          do ik = 1, naky
             vns(ik,:,is,3) = vnmult(1)*spec(is)%vnewk*tunits(ik)/energy**1.5
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
             z0tmp(:,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpdiff(:,isgn,il) &
                  * sqrt(energy(ie))*aj0(:,iglo)
          else
             z0tmp(:,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpa(:,isgn,iglo)*aj0(:,iglo)
          end if
       end do
    end do

    if(use_le_layout) then    
       call gather (g2le, z0tmp, ctmp)
       call solfp_lorentz (ctmp)
       call scatter (g2le, ctmp, z0tmp)   ! z0 is redefined below
    else
       call solfp_lorentz (z0tmp,dum1,dum2)   ! z0 is redefined below
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*z0tmp(:,isgn,iglo)
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
          z0tmp(:,isgn,iglo) = z0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
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
                  * vpdiff(:,isgn,il)*sqrt(energy(ie))
!             gtmp(:,isgn,iglo)  = vns(ik,ie,is,2)*energy(ie)
          end do
       end do

       all = .true.
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

       all = .true.
       call integrate_moment (gtmp, duinv, all)  ! not 1/du yet

!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(1)*spec(is)%vnewk*sqrt(2./pi)
!       end do
    end if

   where (abs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                        !duinv=0 iff vnew=0 so ok to keep duinv=0.
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
             s0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*vpdiff(:,isgn,il)*sqrt(energy(ie)) &
                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
          else
             s0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
          end if
       end do
    end do

    if(use_le_layout) then    
       call gather (g2le, s0tmp, ctmp)
       call solfp_lorentz (ctmp)
       call scatter (g2le, ctmp, s0tmp)   ! s0
    else
       call solfp_lorentz (s0tmp,dum1,dum2)   ! s0
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0s0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*s0tmp(:,isgn,iglo)
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
          s0tmp(:,isgn,iglo) = s0tmp(:,isgn,iglo) - dtmp(:,it,ik,is) &
               * z0tmp(:,isgn,iglo)
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
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*sqrt(energy(ie))*vpdiff(:,isgn,il) &
                  * aj0(:,iglo)*s0tmp(:,isgn,iglo)
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * s0tmp(:,isgn,iglo)
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
          s0tmp(:,isgn,iglo) = s0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
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
             w0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
                  * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                  / bmag
          else
             w0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
                  * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                  / bmag
          end if
       end do
    end do

    if(use_le_layout) then    
       call gather (g2le, w0tmp, ctmp)
       call solfp_lorentz (ctmp)
       call scatter (g2le, ctmp, w0tmp)
    else
       call solfp_lorentz (w0tmp,dum1,dum2)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0w0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*w0tmp(:,isgn,iglo)
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
          w0tmp(:,isgn,iglo) = w0tmp(:,isgn,iglo) - z0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
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
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*sqrt(energy(ie))*vpdiff(:,isgn,il) &
                  * aj0(:,iglo)*w0tmp(:,isgn,iglo)
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * w0tmp(:,isgn,iglo)
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
          w0tmp(:,isgn,iglo) = w0tmp(:,isgn,iglo) - s0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
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
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
               * w0tmp(:,isgn,iglo)
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
          w0tmp(:,isgn,iglo) = w0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

    deallocate (gtmp, duinv, dtmp, vns)


    if(use_le_layout) then    
       allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
       
       ! first set s0le, w0le & z0le
       z_big = 0.0
       z_big = cmplx(real(s0tmp), real(w0tmp))

       ! get rid of z0, s0, w0 now that we've converted to z0le, s0le, w0le
       if (allocated(s0tmp)) deallocate(s0tmp)
       if (allocated(w0tmp)) deallocate(w0tmp)

      call gather (g2le, z_big, ctmp)

       if (.not. allocated(s0le)) then
          allocate (s0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (w0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       end if
       
       s0le = real(ctmp)
       w0le = aimag(ctmp)
       
       ! set z0le
       call gather (g2le, z0tmp, ctmp)
       if (allocated(z0tmp)) deallocate(z0tmp)
       if (.not. allocated(z0le)) allocate (z0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
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
    else
       !Only need the real components (imaginary part should be zero, just
       !use complex arrays to allow reuse of existing integrate routines etc.) 
       if (.not. allocated(s0)) then
          allocate (s0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          s0=real(s0tmp)
          deallocate(s0tmp)
       endif

       if (.not. allocated(w0)) then
          allocate (w0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          w0=real(w0tmp)
          deallocate(w0tmp)
       endif

       if (.not. allocated(z0)) then
          allocate (z0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          z0=real(z0tmp)
          deallocate(z0tmp)
       endif
    end if
    
  end subroutine init_lorentz_conserve

  subroutine init_diffuse_conserve

! Precompute three quantities needed for momentum and energy conservation:
! bz0, bw0, bs0
    
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0, kperp2
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy, al, integrate_moment, negrid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, vpa
    use le_grids, only: g2le, nxi
    use gs2_layouts, only: le_lo
    use redistribute, only: gather, scatter

    implicit none

    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: duinv, dtmp
    real, dimension (:,:,:,:), allocatable :: vns
    integer :: ie, il, ik, is, isgn, iglo, it
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
    complex, dimension (:,:,:), allocatable :: bs0tmp, bw0tmp, bz0tmp    
    logical :: all

    if(use_le_layout) then
      allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end  if

! TO DO: 
! tunits not included anywhere yet


    allocate(bs0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(bw0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(bz0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (dtmp(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (vns(naky,negrid,nspec,2))
    duinv=0.0 ; dtmp=0.0    
    vns(:,:,:,1) = vnmult(2)*delvnew
    vns(:,:,:,2) = vnmult(2)*vnew_s

    ! first obtain 1/du
!    if (conservative) then
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          gtmp(:,isgn,iglo) = energy(ie)*vnmult(2)*vnew_E(ik,ie,is)
       end do
    end do
       
    all = .true.
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
!    else
!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(2)*spec(is)%vnewk*sqrt(2./pi)
!       end do
!    end if

   where (abs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                      !  duinv=0 iff vnew=0 so ok to keep duinv=0.
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
          bz0tmp(:,isgn,iglo) = -code_dt*vnmult(2)*vnew_E(ik,ie,is) &
               * aj0(:,iglo)*duinv(:,it,ik,is)
       end do
    end do

    if(use_le_layout) then    
       call gather (g2le, bz0tmp, ctmp)
       call solfp_ediffuse_le_layout (ctmp)
       call scatter (g2le, ctmp, bz0tmp)   ! bz0 is redefined below
    else
       call solfp_ediffuse_standard_layout (bz0tmp)   ! bz0 is redefined below
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0 f_0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * bz0tmp(:,isgn,iglo)
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
          bz0tmp(:,isgn,iglo) = bz0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
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
!             gtmp(:,isgn,iglo)  = vns(ik,ie,is,2)*energy(ie)

       end do
    end do
       
    all = .true.
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
!    else
!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(2)*spec(is)%vnewk*sqrt(2./pi)
!       end do
!    end if

   where (abs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                      !  duinv=0 iff vnew=0 so ok to keep duinv=0.
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
          bs0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
               * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
!          else
!             bs0(:,isgn,iglo) = -3.0*vns(ik,ie,is,2)*vpa(:,isgn,iglo) &
!                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
!          end if
       end do
    end do

    if(use_le_layout) then
       call gather (g2le, bs0tmp, ctmp)
       call solfp_ediffuse_le_layout (ctmp)
       call scatter (g2le, ctmp, bs0tmp)   ! bs0
    else
       call solfp_ediffuse_standard_layout (bs0tmp)    ! s0
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0s0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * bs0tmp(:,isgn,iglo)
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
          bs0tmp(:,isgn,iglo) = bs0tmp(:,isgn,iglo) - dtmp(:,it,ik,is) &
               * bz0tmp(:,isgn,iglo)
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
               * bs0tmp(:,isgn,iglo)
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
          bs0tmp(:,isgn,iglo) = bs0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
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
          bw0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
               * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
               / bmag
       end do
    end do

    if(use_le_layout) then    
       call gather (g2le, bw0tmp, ctmp)
       call solfp_ediffuse_le_layout (ctmp)
       call scatter (g2le, ctmp, bw0tmp)
    else
       call solfp_ediffuse_standard_layout (bw0tmp)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0w0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * bw0tmp(:,isgn,iglo)
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
          bw0tmp(:,isgn,iglo) = bw0tmp(:,isgn,iglo) - bz0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
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
               * bw0tmp(:,isgn,iglo)
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
          bw0tmp(:,isgn,iglo) = bw0tmp(:,isgn,iglo) - bs0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
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
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
               * bw0tmp(:,isgn,iglo)
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
          bw0tmp(:,isgn,iglo) = bw0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do

    deallocate (gtmp, duinv, dtmp, vns)

    if(use_le_layout) then
       allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))    
       z_big=cmplx(real(bs0tmp),real(bw0tmp))
       deallocate (bs0tmp)
       deallocate (bw0tmp)
       call gather (g2le, z_big, ctmp)
       deallocate(z_big)
       if (.not. allocated(bs0le)) allocate (bs0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       if (.not. allocated(bw0le)) allocate (bw0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       bs0le = real(ctmp)
       bw0le = aimag(ctmp)

       call gather (g2le, bz0tmp, ctmp)
       deallocate (bz0tmp)
       if (.not. allocated(bz0le)) allocate (bz0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       bz0le = real(ctmp)

       if (.not. allocated(aj0le)) then
          allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))    
          ! next set aj0le & aj1l
          z_big(:,1,:) = cmplx(aj0,aj1)
          z_big(:,2,:) = z_big(:,1,:)
          call gather (g2le, z_big, ctmp)
          deallocate (z_big)
          allocate (aj0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (aj1le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
          aj0le = real(ctmp)
          aj1le = aimag(ctmp)
       end if
       
       deallocate (ctmp)
    else
       !Only need the real components (imaginary part should be zero, just
       !use complex arrays to allow reuse of existing integrate routines etc.)
       if (.not. allocated(bs0)) then
          allocate (bs0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          bs0=real(bs0tmp)
          deallocate(bs0tmp)
       endif

       if (.not. allocated(bw0)) then
          allocate (bw0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          bw0=real(bw0tmp)
          deallocate(bw0tmp)
       endif

       if (.not. allocated(bz0)) then
          allocate (bz0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          bz0=real(bz0tmp)
          deallocate(bz0tmp)
       endif
    end if
 end subroutine init_diffuse_conserve

 subroutine init_vnew (hee)
   use species, only: nspec, spec, electron_species, has_electron_species
   use le_grids, only: negrid, energy, w
   use kt_grids, only: naky, ntheta0, kperp2
   use theta_grid, only: ntgrid
   use run_parameters, only: zeff, tunits
    use constants, only: pi
    use spfunc, only: erf => erf_ext
    real, dimension (:), intent (out) :: hee
    real,dimension (negrid)::heevth, hsg, hsgvth
    integer :: ik, ie, is, it, ig
    real :: k4max
    real :: vl, vr, dv2l, dv2r
!    real :: erf ! this is needed for PGI: RN


    do ie = 1, negrid
       hee(ie) = exp(-energy(ie))/sqrt(pi*energy(ie)) &
            + (1.0 - 0.5/energy(ie))*erf(sqrt(energy(ie)))

!>MAB
! hsg is the G of Hirshman and Sigmar
! added to allow for momentum conservation with energy diffusion
       hsg(ie) = hsg_func(sqrt(energy(ie)))
!<MAB
    end do

!heevth is hee but using only thermal velocity (energy independent)

    do ie = 1, negrid
       heevth(ie) = exp(-1.0)/sqrt(pi) + 0.5*erf(1.0)

!>MAB
! hsg is the G of Helander and Sigmar
! added to allow for momentum conservation with energy diffusion
       hsgvth(ie) = hsg_func(1.0)
!<MAB
    end do                                                                  

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
                        *(zeff+heevth(ie))*0.5*tunits(ik)                   
                   vnew_s(ik,ie,is) = spec(is)%vnewk &
                        *hsgvth(ie)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = spec(is)%vnewk &
                        *heevth(ie)*tunits(ik)                   
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                           - 2.0*vnew_D(ik,ie,is)
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                        *(zeff + hee(ie))*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(energy(ie)) &
                        *hsg(ie)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                        *hee(ie)*tunits(ik)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = energy(ie)*(vnew_s(ik,ie,is)*(2.0-0.5/energy(ie)) &
                           - 2.0*vnew_D(ik,ie,is))
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                end if
                if (ei_coll_only) then
                   vnew(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
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
                        *heevth(ie)*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk &
                        *hsgvth(ie)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                           - 2.0*vnew_D(ik,ie,is)
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                        *hee(ie)*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(energy(ie)) &
                        *hsg(ie)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = energy(ie)*(vnew_s(ik,ie,is)*(2.0-0.5/energy(ie)) &
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
             vr = 0.5*(sqrt(energy(ie+1)) + sqrt(energy(ie)))
             vl = 0.5*(sqrt(energy(ie  )) + sqrt(energy(ie-1)))
             dv2r = (energy(ie+1) - energy(ie)) / (sqrt(energy(ie+1)) - sqrt(energy(ie)))
             dv2l = (energy(ie) - energy(ie-1)) / (sqrt(energy(ie)) - sqrt(energy(ie-1)))
          
             vnew_E(:,ie,is) = spec(is)%vnewk*tunits*(vl*exp(-vl**2)*dv2l*hsg_func(vl) &
                  - vr*exp(-vr**2)*dv2r*hsg_func(vr)) / (sqrt(pi)*w(ie))
             delvnew(:,ie,is) = spec(is)%vnewk*tunits*(vl*exp(-vl**2)*hsg_func(vl) &
                     - vr*exp(-vr**2)*hsg_func(vr)) / (sqrt(pi*energy(ie))*w(ie))
          end do

          ! boundary at v = 0
          vr = 0.5*(sqrt(energy(2)) + sqrt(energy(1)))
          dv2r = (energy(2) - energy(1)) / (sqrt(energy(2)) - sqrt(energy(1)))

          vnew_E(:,1,is) = -spec(is)%vnewk*tunits*vr*exp(-vr**2)*hsg_func(vr)*dv2r &
               / (sqrt(pi)*w(1))
          delvnew(:,1,is) = -spec(is)%vnewk*tunits*vr*exp(-vr**2)*hsg_func(vr) &
               / (sqrt(pi*energy(1))*w(1))

          ! boundary at v -> infinity
          vl = 0.5*(sqrt(energy(negrid)) + sqrt(energy(negrid-1)))
          dv2l = (energy(negrid) - energy(negrid-1)) / (sqrt(energy(negrid)) - sqrt(energy(negrid-1)))

          vnew_E(:,negrid,is) = spec(is)%vnewk*tunits*vl*exp(-vl**2)*hsg_func(vl)*dv2l &
               / (sqrt(pi)*w(negrid))
          delvnew(:,negrid,is) = spec(is)%vnewk*tunits*vl*exp(-vl**2)*hsg_func(vl) &
               / (sqrt(pi*energy(negrid))*w(negrid))

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

  subroutine init_ediffuse (vnmult_target)
    use le_grids, only: negrid, nxi
    use le_grids, only: forbid, ixi_to_il
    use egrid, only: zeroes, x0
    use gs2_layouts, only: le_lo, e_lo, il_idx
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx

    implicit none
    
    real, intent (in), optional :: vnmult_target

    integer :: ie, is, ik, il, ig, it
    real, dimension (:), allocatable :: aa, bb, cc, xe
    integer :: ile, ixi
    integer :: ielo

    allocate (aa(negrid), bb(negrid), cc(negrid))
    allocate (xe(negrid))

    ! want to use x variables instead of e because we want conservative form
    ! for the x-integration
    xe(1:negrid-1) = zeroes
    xe(negrid) = x0

    if (present(vnmult_target)) then
       vnmult(2) = max (vnmult_target, 1.0)
    end if


    if(use_le_layout) then    

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
            il = ixi_to_il(ig,ixi)
            if (forbid(ig,il)) cycle
            call get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, xe)
            ec1le(ixi,:,ile) = cc
            ebetaale(ixi,1,ile) = 1.0 / bb(1)
            do ie = 1, negrid-1
               eqle(ixi,ie+1,ile) = aa(ie+1) * ebetaale(ixi,ie,ile)
               ebetaale(ixi,ie+1,ile) = 1.0 / ( bb(ie+1) - eqle(ixi,ie+1,ile) &
                    * ec1le(ixi,ie,ile) )
            end do
         end do
      end do

    else

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
         call get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, xe)
         ec1(:,ielo) = cc

! fill in the arrays for the tridiagonal
         ebetaa(1,ielo) = 1.0/bb(1)
         do ie = 1, negrid-1
            eql(ie+1,ielo) = aa(ie+1)*ebetaa(ie,ielo)
            ebetaa(ie+1,ielo) = 1.0/(bb(ie+1)-eql(ie+1,ielo)*ec1(ie,ielo))
         end do

      end do

    end if

    deallocate(aa, bb, cc, xe)

  end subroutine init_ediffuse

  subroutine get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, xe)

    use species, only: spec
    use run_parameters, only: tunits
    use theta_grid, only: bmag
    use le_grids, only: al, negrid, w
    use spfunc, only: erf => erf_ext
    use constants, only: pi
    use kt_grids, only: kperp2
    use gs2_time, only: code_dt

    implicit none

    integer, intent (in) :: ig, ik, it, il, is
    real, dimension (:), intent (in) :: xe
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
         
          capgr = 2.0*exp(-xer**2)*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/(xer*sqrt(pi))
          capgl = 2.0*exp(-xel**2)*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/(xel*sqrt(pi))
          
          ee = 0.125*(1.-slb1**2)*vnew_s(ik,ie,is) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(ie) = -0.25*vn*code_dt*capgr/(w(ie)*(xe2 - xe1))
          aa(ie) = -0.25*vn*code_dt*capgl/(w(ie)*(xe1 - xe0))
          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*code_dt
          
             ! coefficients for entropy heating calculation
!             if (heating) then
!                dd(ie) =vnc*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!                hh(ie) =vnh*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!             end if
       end do

       ! boundary at v = 0
       xe1 = xe(1)
       xe2 = xe(2)
       
       xer = (xe1 + xe2)*0.5
       
       capgr = 2.0*exp(-xer**2)*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer/sqrt(pi)
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,1,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -0.25*vn*code_dt*capgr/(w(1)*(xe2 - xe1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*code_dt
       
!          if (heating) then
!             dd(1) =vnc*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!             hh(1) =vnh*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!          end if

       ! boundary at v = infinity
       xe0 = xe(negrid-1)
       xe1 = xe(negrid)
       
       xel = (xe1 + xe0)*0.5
       
       capgl = 2.0*exp(-xel**2)*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel/sqrt(pi)

       ee = 0.125*(1.-slb1**2)*vnew_s(ik,negrid,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(negrid) = 0.0
       aa(negrid) = -0.25*vn*code_dt*capgl/(w(negrid)*(xe1 - xe0))
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
          
          capgr = 0.5*exp(xe1**2-xer**2)/xe1**2*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer
          capgl = 0.5*exp(xe1**2-xel**2)/xe1**2*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel
          
          ee = 0.125*(1.-slb1**2)*vnew_s(ik,ie,is) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(ie) = -vn*code_dt*capgr/((xer-xel)*(xe2 - xe1))
          aa(ie) = -vn*code_dt*capgl/((xer-xel)*(xe1 - xe0))
          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*code_dt
          
          ! coefficients for entropy heating calculation
!             if (heating) then
!                dd(ie) =vnc*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!                hh(ie) =vnh*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!             end if
       end do

       ! boundary at xe = 0
       xe1 = xe(1)
       xe2 = xe(2)
       
       xer = (xe1 + xe2)*0.5
       
       capgr = 0.5*exp(xe1**2-xer**2)/xe1**2*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,1,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -vn*code_dt*capgr/(xer*(xe2 - xe1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*code_dt
       
!          if (heating) then
!             dd(1) =vnc*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!             hh(1) =vnh*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!          end if

       ! boundary at xe = 1
       xe0 = xe(negrid-1)
       xe1 = xe(negrid)
       
       xel = (xe1 + xe0)*0.5
       
       capgl = 0.5*exp(xe1**2-xel**2)/xe1**2*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel
       
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

  subroutine init_lorentz (vnmult_target)
    use le_grids, only: negrid, jend, ng2, nxi
    use gs2_layouts, only: le_lo
    use gs2_layouts, only: lz_lo, g_lo
    use gs2_layouts, only: ig_idx, ik_idx, ie_idx, is_idx, it_idx
    use redistribute, only: gather
    use le_grids, only: g2le
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0, aj1

    implicit none

    real, intent (in), optional :: vnmult_target

    integer :: ig, il, it, ik, ie, is, je, te2
    real, dimension (:), allocatable :: aa, bb, cc
    real, dimension (:), allocatable :: dd, hh
    integer :: ile
    integer :: ilz
    complex, dimension (:,:,:), allocatable :: ctmp, z_big

    allocate (aa(nxi+1), bb(nxi+1), cc(nxi+1), dd(nxi+1), hh(nxi+1))

    call init_vpdiff

    if(use_le_layout) then

       !Make sure aj0le is available if drag term included
       if(drag)then
          if(.not.allocated(aj0le))then
             allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
             allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))

             ! next set aj0le & aj1l
             z_big(:,1,:) = cmplx(aj0,aj1)
             z_big(:,2,:) = z_big(:,1,:)

             call gather (g2le, z_big, ctmp)
             deallocate (z_big)

             allocate (aj0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
             allocate (aj1le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))

             aj0le = real(ctmp)
             aj1le = aimag(ctmp)

             deallocate (ctmp)
          endif
       endif

      if (heating .and. .not. allocated(glec)) then
         allocate (glec(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
         glec = 0.0
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

    else

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

    end if

    deallocate (aa, bb, cc, dd, hh)

  end subroutine init_lorentz

  subroutine get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)

    use species, only: spec
    use le_grids, only: al, energy, xi, ng2
    use le_grids, only: wl, jend, al
    use gs2_time, only: code_dt
    use kt_grids, only: kperp2
    use theta_grid, only: bmag
    use run_parameters, only: tunits

    implicit none

    real, dimension (:), intent (out) :: aa, bb, cc, dd, hh
    integer, intent (in) :: ig, ik, it, ie, is

    integer :: il, je, te, te2, teh
    real :: slb0, slb1, slb2, slbl, slbr, vn, vnh, vnc, ee
    real deltaxi

    je = jend(ig)
!
!CMR, 17/2/2014:
!         te, te2, teh: indices in xi, which runs from +1 -> -1.
!         te   :  index of minimum xi value >= 0.
!         te2  :  total #xi values = index of minimum xi value (= -1) 
!         teh  :  index of minimum xi value > 0.  
!                 teh = te if no bouncing particle at this location
!              OR teh = te-1 if there is a bouncing particle
!
    if (special_wfb_lorentz) then
!CMRDDGC, 17/2/2014:
!   This clause is appropriate for Lorentz collisons with 
!         SPECIAL (unphysical) treatment of wfb at its bounce point
       if (je <= ng2+1) then
          te = ng2
          te2 = 2*ng2
          teh = ng2
       else
          te = je
          te2 = 2*je-1
          teh = je-1
       end if
    else
!CMRDDGC, 17/2/2014:
!   This clause is appropriate for Lorentz collisons with 
!         STANDARD treatment of wfb at its bounce point
       if (je <= ng2) then
          te = ng2
          te2 = 2*ng2
          teh = ng2
       else
          te = je
          te2 = 2*je-1
          teh = je-1
       end if
    endif 
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
          
          ee = 0.5*energy(ie)*(1+slb1**2) &
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
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
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
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
!
!CMR, 6/3/2014:  
! STANDARD treatment of pitch angle scattering must resolve T-P boundary.
! NEED special_wfb= .false. to resolve T-P boundary at wfb bounce point 
!     (special_wfb= .true. AVOIDS TP boundary at wfb bounce point)
!
! Original code (pre-r2766) used eq.(42) Barnes et al, Phys Plasmas 16, 072107 
! (2009), with pitch angle weights to enhance accuracy in derivatives.
! NB THIS FAILS at wfb bounce point, giving aa=cc=infinite, 
!    because weight wl(ig,il)=0 at wfb bounce point.
!    UNPHYSICAL as d/dxi ((1-xi^2)g) IS NOT resolved numerically for wl=0.
! MUST accept limitations of the grid resolution and USE FINITE coefficients! 
! FIX here by setting a FINITE width of the trapped region at wfb B-P 
!              deltaxi=xi(ig,ng2)-xi(ig,ng2+2)
! ASIDE: NB    deltaxi=wl is actually 2*spacing in xi !!!
!              which explains upfront factor 2 in definition of aa, cc
       deltaxi=wl(ig,il)
       if (.not. special_wfb_lorentz .and. deltaxi .eq. 0. .and. il .eq. ng2+1) then 
          deltaxi=xi(ig,ng2)-xi(ig,ng2+2)
       endif
       cc(il) = 2.0*vn*code_dt*(1.0 - slbr**2)/(deltaxi*(slb2 - slb1))
       aa(il) = 2.0*vn*code_dt*(1.0 - slbl**2)/(deltaxi*(slb1 - slb0))
       bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
       
       if (heating) then
          dd(il) =vnc*(-2.0*(1.0-slbr**2)/(deltaxi*(slb2-slb1)) + ee)
          hh(il) =vnh*(-2.0*(1.0-slbr**2)/(deltaxi*(slb2-slb1)) + ee)
       end if
!CMRend

    case (lorentz_scheme_old)
       
       do il = 2, te-1
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1))) 
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il))) 
          slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1))) 
          
          slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
          slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
          
          ee = 0.5*energy(ie)*(1+slb1**2) &
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
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
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
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
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
    cc(te+1:te2) = aa(teh:1:-1)
       
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

  subroutine solfp1_standard_layout (g, g1, gc1, gc2, diagnostics, gtoc, ctog)

    use gs2_layouts, only: g_lo, it_idx, ik_idx, ie_idx, is_idx
    use theta_grid, only: ntgrid
    use run_parameters, only: beta
    use gs2_time, only: code_dt
    use le_grids, only: energy
    use species, only: spec, electron_species
    use dist_fn_arrays, only: vpa, aj0
    use fields_arrays, only: aparnew
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter, kperp2
    ! TMP FOR TESTING -- MAB
!    use mp, only: proc0

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1, gc1, gc2
    integer, optional, intent (in) :: diagnostics
    logical, optional, intent (in) :: gtoc, ctog
    integer :: ig, it, ik, ie, is, iglo
!CMR, 12/9/2013: 
!CMR   New logical optional input parameters gtoc, ctog used to set
!CMR   flags (g_to_c and c_to_g) to control whether redistributes required
!CMR   to map g_lo to collision_lo, and collision_lo to g_lo.
!CMR   All redistributes are performed by default.
!CMR  

    logical :: g_to_c, c_to_g

    if (present(gtoc)) then 
       g_to_c=gtoc 
    else 
       g_to_c=.true.
    endif

    if (present(ctog)) then 
       c_to_g=ctog 
    else 
       c_to_g=.true.
    endif


    ! TMP FOR TESTING -- MAB
!    integer :: t0, t1, t2, t3, t4, t5, tr
!    real :: t1tot = 0., t2tot = 0., t3tot = 0., t4tot = 0., t5tot = 0.

    heating_flag = heating .and. present(diagnostics)

    select case (collision_model_switch)
    case (collision_model_full)

       ! TMP FOR TESTING -- MAB
!       if (proc0) call system_clock (count=t0, count_rate=tr)

       call solfp_ediffuse_standard_layout (g)

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t1)
!          t1tot = t1tot + real(t1-t0)/tr
!       end if

       if (conserve_moments) call conserve_diffuse (g, g1)

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t2)
!          t2tot = t2tot + real(t2-t1)/tr
!       end if

       if (drag) then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             is = is_idx(g_lo,iglo)
             if (spec(is)%type /= electron_species) cycle
             it = it_idx(g_lo,iglo)
             ik = ik_idx(g_lo,iglo)
             if(kwork_filter(it,ik)) cycle
             ie = ie_idx(g_lo,iglo)
             do ig = -ntgrid, ntgrid
                g(ig,:,iglo) = g(ig,:,iglo) + ieqzip(it,ik) * &
                     vnmult(1)*spec(is)%vnewk*code_dt &
                     * vpa(ig,:,iglo)*kperp2(ig,it,ik)*aparnew(ig,it,ik)*aj0(ig,iglo) &
                     / (beta*spec(is)%stm*energy(ie)**1.5)
                ! probably need 1/(spec(is_ion)%z*spec(is_ion)%dens) above
             end do
          end do
       end if

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t3)
!          t3tot = t3tot + real(t3-t2)/tr
!       end if

       if (heating_flag) then
          call solfp_lorentz (g, gc1, gc2, diagnostics)
       else
          call solfp_lorentz (g, gc1, gc2)
       end if

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t4)
!          t4tot = t4tot + real(t4-t3)/tr
!       end if

       if (conserve_moments) call conserve_lorentz (g, g1)

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t5)
!          t5tot = t5tot + real(t5-t4)/tr
!          write (*,'(a10,5e14.5)') 'solfp1: ', t1tot, t2tot, t3tot, t4tot, t5tot
!       end if

    case (collision_model_lorentz,collision_model_lorentz_test)

       if (drag) then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             is = is_idx(g_lo,iglo)
             if (spec(is)%type /= electron_species) cycle
             it = it_idx(g_lo,iglo)
             ik = ik_idx(g_lo,iglo)
             if(kwork_filter(it,ik))cycle
             ie = ie_idx(g_lo,iglo)
             do ig = -ntgrid, ntgrid
                g(ig,:,iglo) = g(ig,:,iglo) + ieqzip(it,ik) * &
                     vnmult(1)*spec(is)%vnewk*code_dt &
                     * vpa(ig,:,iglo)*kperp2(ig,it,ik)*aparnew(ig,it,ik)*aj0(ig,iglo) &
                     / (beta*spec(is)%stm*energy(ie)**1.5)
             end do
          end do
       end if
       
       if (heating_flag) then
          call solfp_lorentz (g, gc1, gc2, diagnostics)
       else
          call solfp_lorentz (g, gc1, gc2)
       end if

       if (conserve_moments) call conserve_lorentz (g, g1)

    case (collision_model_ediffuse)

       call solfp_ediffuse_standard_layout (g)
       if (conserve_moments) call conserve_diffuse (g, g1)

    end select
   
  end subroutine solfp1_standard_layout

  subroutine solfp1_le_layout (gle, diagnostics)

    use gs2_layouts, only: le_lo, it_idx, ik_idx, ig_idx, is_idx
    use theta_grid, only: bmag
    use run_parameters, only: beta, ieqzip
    use gs2_time, only: code_dt
    use le_grids, only: energy, negrid, nxi, ixi_to_il, ixi_to_isgn, sgn, al
    use species, only: spec, electron_species
    use fields_arrays, only: aparnew
    use kt_grids, only: kwork_filter, kperp2
    ! TMP FOR TESTING -- MAB
!    use mp, only: proc0

    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    integer, optional, intent (in) :: diagnostics

    integer :: ig, it, ik, il, ie, is, ile, ixi, isgn

    ! TMP FOR TESTING -- MAB
!    integer :: t0, t1, t2, t3, t4, t5, tr
!    real :: t1tot = 0., t2tot = 0., t3tot = 0., t4tot = 0., t5tot = 0.

    heating_flag = heating .and. present(diagnostics)

    select case (collision_model_switch)
    case (collision_model_full)

       ! TMP FOR TESTING -- MAB
!       if (proc0) call system_clock (count=t0, count_rate=tr)

       call solfp_ediffuse_le_layout (gle)

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t1)
!          t1tot = t1tot + real(t1-t0)/tr
!       end if

       if (conserve_moments) call conserve_diffuse (gle)

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t2)
!          t2tot = t2tot + real(t2-t1)/tr
!       end if

       if (drag) then
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             is = is_idx(le_lo,ile)
             if (spec(is)%type /= electron_species) cycle
             it = it_idx(le_lo,ile)
             ik = ik_idx(le_lo,ile)
             if(kwork_filter(it,ik)) cycle
             ig = ig_idx(le_lo,ile)
             do ie = 1, negrid
                do ixi = 1, nxi
                   il = ixi_to_il(ig,ixi)
                   isgn = ixi_to_isgn(ig,ixi)
                   gle(ixi,ie,ile) = gle(ixi,ie,ile) + ieqzip(it,ik) * &
                        vnmult(1)*spec(is)%vnewk*code_dt &
                        * kperp2(ig,it,ik)*aparnew(ig,it,ik)*aj0le(ixi,ie,ile) &
                        / (beta*spec(is)%stm*energy(ie)) &
                        * sgn(isgn)*sqrt(MAX(1.0-al(il)*bmag(ig),0.0))
                   ! probably need 1/(spec(is_ion)%z*spec(is_ion)%dens) above
                end do
             end do
          end do
       end if

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t3)
!          t3tot = t3tot + real(t3-t2)/tr
!       end if

       if (heating_flag) then
          call solfp_lorentz (gle, diagnostics)
       else
          call solfp_lorentz (gle)
       end if

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t4)
!          t4tot = t4tot + real(t4-t3)/tr
!       end if

       if (conserve_moments) call conserve_lorentz (gle)

       ! TMP FOR TESTING -- MAB
!       if (proc0) then
!          call system_clock (count=t5)
!          t5tot = t5tot + real(t5-t4)/tr
!          write (*,'(a10,5e14.5)') 'solfp1: ', t1tot, t2tot, t3tot, t4tot, t5tot
!       end if

    case (collision_model_lorentz,collision_model_lorentz_test)

       if (drag) then
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             is = is_idx(le_lo,ile)
             if (spec(is)%type /= electron_species) cycle
             it = it_idx(le_lo,ile)
             ik = ik_idx(le_lo,ile)
             if(kwork_filter(it,ik)) cycle
             ig = ig_idx(le_lo,ile)
             do ie = 1, negrid
                do ixi = 1, nxi
                   il = ixi_to_il(ig,ixi)
                   isgn = ixi_to_isgn(ig,ixi)
                   gle(ixi,ie,ile) = gle(ixi,ie,ile) + ieqzip(it,ik) * &
                        vnmult(1)*spec(is)%vnewk*code_dt &
                        * kperp2(ig,it,ik)*aparnew(ig,it,ik)*aj0le(ixi,ie,ile) &
                        / (beta*spec(is)%stm*energy(ie)) &
                        * sgn(isgn)*sqrt(MAX(1.0-al(il)*bmag(ig),0.0))
                   ! probably need 1/(spec(is_ion)%z*spec(is_ion)%dens) above
                end do
             end do
          end do
       end if

       if (heating_flag) then
          call solfp_lorentz (gle, diagnostics)
       else
          call solfp_lorentz (gle)
       end if

       if (conserve_moments) call conserve_lorentz (gle)

    case (collision_model_ediffuse)

       call solfp_ediffuse_le_layout (gle)

       if (conserve_moments) call conserve_diffuse (gle)

    end select
   
  end subroutine solfp1_le_layout


  subroutine conserve_lorentz_standard_layout (g, g1)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: energy, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y1, v2y2

    integer :: isgn, iglo, ik, ie, il, is, it
    logical :: all = .true. 

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v1y1(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v2y2(-ntgrid:ntgrid, ntheta0, naky, nspec))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (vns(naky,negrid,nspec))
    vns = vnmult(1)*vnew_D

    if (drag) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          if(kwork_filter(it,ik))cycle
          do isgn = 1, 2
             ! v0 = vpa J0 f0, y0 = g
             gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*g(:,isgn,iglo)
          end do
       end do
       
       call integrate_moment (gtmp, v0y0, all)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          if(kwork_filter(it,ik))cycle
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g1(:,isgn,iglo) = g(:,isgn,iglo) - ieqzip(it,ik)*v0y0(:,it,ik,is) &
                  * z0(:,isgn,iglo)
          end do
       end do

    else

       g1 = g

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nud vpa J0 f0, y1 = g1
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is)*sqrt(energy(ie))*vpdiff(:,isgn,il) &
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
       if(kwork_filter(it,ik))cycle
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
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = nud vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*energy(ie)*al(il)*aj1(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v2y2, all)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g(:,isgn,iglo) = g1(:,isgn,iglo) - ieqzip(it,ik)*v2y2(:,it,ik,is) &
               * w0(:,isgn,iglo)
       end do
    end do

    deallocate (vns, v0y0, v1y1, v2y2, gtmp)

  end subroutine conserve_lorentz_standard_layout


  subroutine conserve_lorentz_le_layout (gle)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: energy, al, integrate_moment, negrid
    use le_grids, only: ixi_to_il, ixi_to_isgn
    use run_parameters, only: ieqzip
    use le_grids, only: sgn, nxi
    use gs2_layouts, only: le_lo, ig_idx
    use redistribute, only: scatter
    use run_parameters, only: tunits
    use theta_grid, only: bmag
    use kt_grids, only: kwork_filter
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:,:), allocatable :: vpanud
    complex, dimension (:), allocatable :: v0y0, v1y1, v2y2

    integer :: ig, isgn, ik, ie, il, is, it
    integer :: ile, ixi

    allocate (v0y0(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v1y1(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v2y2(le_lo%llim_proc:le_lo%ulim_alloc))

    ! Let's work on gle directly instead of g for the moment
    allocate (gtmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc)) ; gtmp = 0.0
    allocate (vpanud(-ntgrid:ntgrid, nxi+1, negrid+1, nspec)) ; vpanud = 0.0

!    if (resistivity) then
    if (drag) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

       ! initially define vpanud to be vpa
       do is = 1, nspec
          do ie = 1, negrid
             do ixi = 1, nxi
                do ig = -ntgrid, ntgrid
                   il = ixi_to_il(ig,ixi)
                   isgn = ixi_to_isgn(ig,ixi)
                   vpanud(ig,ixi,ie,is) = sgn(isgn)*sqrt(MAX(1.0-al(il)*bmag(ig),0.0)*energy(ie))
                end do
             end do
          end do
       end do
       
       ! v0 = vpa J0 f0, y0 = gle
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          is = is_idx(le_lo,ile)
          ig = ig_idx(le_lo,ile)
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          if(kwork_filter(it,ik)) cycle
          gtmp(:,:,ile) = vpanud(ig,:,:,is) * aj0le(:,:,ile) * gle(:,:,ile)
       end do
       call integrate_moment (le_lo, gtmp, v0y0)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

       do ile = le_lo%llim_proc, le_lo%ulim_proc
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          if(kwork_filter(it,ik)) cycle
          gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)* z0le(:,:,ile) * v0y0(ile)
       end do

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    ! redefine vpanud to be vpa * nu_D
    if (conservative) then
       do is = 1, nspec
          do ie = 1, negrid
             do ixi = 1, nxi
                do ig = -ntgrid, ntgrid
                   il = ixi_to_il(ig,ixi)
                   isgn = ixi_to_isgn(ig,ixi)
                   vpanud(ig,ixi,ie,is) = vpdiff(ig,isgn,il) * sqrt(energy(ie)) * vnmult(1)*vnew_D(1,ie,is)/tunits(1)
                end do
             end do
          end do
       end do
    else
       do is = 1, nspec
          do ie = 1, negrid
             do ixi = 1, nxi
                do ig = -ntgrid, ntgrid
                   il = ixi_to_il(ig,ixi)
                   isgn = ixi_to_isgn(ig,ixi)
                   vpanud(ig,ixi,ie,is) = sgn(isgn)*sqrt(MAX(1.0-al(il)*bmag(ig),0.0)*energy(ie))*vnmult(1)*vnew_D(1,ie,is)/tunits(1)
                end do
             end do
          end do
       end do
    end if

    ! v1 = nud vpa J0 f0, y1 = gle
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       is = is_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       gtmp(:,:,ile) = vpanud(ig,:,:,is) * tunits(ik) * aj0le(:,:,ile) &
            * gle(:,:,ile)
    end do

    call integrate_moment (le_lo, gtmp, v1y1)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle

       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*s0le(:,:,ile) * v1y1(ile)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       is = is_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       do ixi=1, nxi
          il = ixi_to_il(ig,ixi)
          ! aj1vp2 = 2 * J1(arg)/arg * vperp^2
          gtmp(ixi,:negrid,ile) = vnmult(1) * vnew_D(ik,:negrid,is) * aj1le(ixi,:negrid,ile) &
               * energy(:) * al(il) * gle(ixi,:negrid,ile)
       end do
    end do

    call integrate_moment (le_lo, gtmp, v2y2)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*w0le(:,:,ile) * v2y2(ile)
    end do

    deallocate (gtmp, vpanud, v0y0, v1y1, v2y2)

  end subroutine conserve_lorentz_le_layout


  subroutine conserve_diffuse_standard_layout (g, g1)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: energy, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns

    integer :: isgn, iglo, ik, ie, il, is, it 
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y1, v2y2    
    logical :: all = .true.

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v1y1(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v2y2(-ntgrid:ntgrid, ntheta0, naky, nspec))

    allocate (vns(naky,negrid,nspec))
    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    vns = vnmult(2)*delvnew

    !This is needed to to ensure the it,ik values we don't set aren't included
    !in the integral (can also be enforced in integrate_moment routine)
    if(any(kwork_filter)) gtmp=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
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
       if(kwork_filter(it,ik))cycle
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
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
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
       if(kwork_filter(it,ik))cycle
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
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = (nus-nud) vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*energy(ie)*al(il)*aj1(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do

    call integrate_moment (gtmp, v2y2, all)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g(:,isgn,iglo) = g1(:,isgn,iglo) - ieqzip(it,ik)*v2y2(:,it,ik,is) &
               * bw0(:,isgn,iglo)
       end do
    end do

    deallocate (vns, v0y0, v1y1, v2y2, gtmp)

  end subroutine conserve_diffuse_standard_layout

  subroutine conserve_diffuse_le_layout (gle)

    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: energy, al, integrate_moment, negrid
    use le_grids, only: ixi_to_il, ixi_to_isgn
    use run_parameters, only: ieqzip
    use le_grids, only: forbid, sgn, speed, nxi
    use gs2_layouts, only: le_lo, ig_idx, ik_idx, it_idx
    use redistribute, only: scatter
    use theta_grid, only: bmag
    use run_parameters, only: tunits
    use kt_grids, only: kwork_filter
    ! TMP FOR TESTING -- MAB
!    use mp, only: proc0

    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns

    integer :: ig, isgn, ik, ie, il, is, it
    integer :: ile, ixi
    real, dimension (:,:,:,:), allocatable :: vpadelnu
    complex, dimension (:), allocatable :: v0y0, v1y1, v2y2

    ! TMP FOR TESTING -- MAB
!    integer :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, tr
!    real :: t1tot=0., t2tot = 0., t3tot = 0., t4tot = 0., t5tot = 0., t6tot = 0., t7tot = 0., t8tot = 0., t9tot = 0., t10tot=0.

    allocate (v0y0(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v1y1(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v2y2(le_lo%llim_proc:le_lo%ulim_alloc))

    allocate (gtmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (vpadelnu(-ntgrid:ntgrid, nxi+1, negrid+1, nspec)) ; vpadelnu = 0.0
    allocate (vns(naky,negrid,nspec))

    !This is needed to to ensure the it,ik values we don't set aren't included
    !in the integral (can also be enforced in integrate_moment routine)
    !if(any(kwork_filter)) gtmp=0.
    gtmp=0.

    if (.not. allocated(vpatmp)) then
       allocate(vpatmp(-ntgrid:ntgrid,nxi)) ; vpatmp = 0.0
       do ixi = 1, nxi
          do ig = -ntgrid, ntgrid
             il = ixi_to_il(ig,ixi)
             isgn = ixi_to_isgn(ig,ixi)
             if (.not. forbid(ig,il)) &
                  vpatmp(ig,ixi) = sgn(isgn)*sqrt(max(0.0,(1.0-al(il)*bmag(ig))))
          end do
       end do
    end if
    
    vns = vnmult(2)*vnew_E
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0
    
    ! TMP FOR TESTING -- MAB
!    if (proc0) call system_clock (count=t0, count_rate=tr)

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       do ie=1, negrid
          do ixi = 1, nxi
             gtmp(ixi,ie,ile) = vns(ik,ie,is)*aj0le(ixi,ie,ile)*gle(ixi,ie,ile)
          end do
       end do
    end do
    
    call integrate_moment (le_lo, gtmp, v0y0)    ! v0y0
!    call integrate_moment (le_lo, gle*aj0le, v0y0, vns)    ! v0y0

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t1)
!       t1tot = t1tot + real(t1-t0)/tr
!    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle

       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*v0y0(ile)*bz0le(:,:,ile)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    vns(1,:,:) = vnmult(2)*delvnew(1,:,:)/tunits(1)

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t2)
!       t2tot = t2tot + real(t2-t1)/tr
!    end if

    do is = 1, nspec
       do ie = 1, negrid
          vpadelnu(:,:nxi,ie,is) = vns(1,ie,is) * vpatmp * speed(ie)
       end do
    end do

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t3)
!       t3tot = t3tot + real(t3-t2)/tr
!    end if

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle

       gtmp (:,:,ile) = vpadelnu(ig,:,:,is) * tunits(ik) * aj0le(:,:,ile) * gle(:,:,ile)
    end do

    call integrate_moment (le_lo, gtmp, v1y1)    ! v1y1

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t4)
!       t4tot = t4tot + real(t4-t3)/tr
!    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle

       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*bs0le(:,:,ile) * v1y1(ile)
    end do

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t5)
!       t5tot = t5tot + real(t5-t4)/tr
!    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    do is = 1, nspec
       do ik = 1, naky
          vns(ik,:,is) = delvnew(ik,:,is)*energy 
       end do
    end do

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle

       ig = ig_idx(le_lo,ile)
       do ie=1, negrid
          do ixi = 1, nxi
             il = ixi_to_il(ig,ixi)
             gtmp(ixi,ie,ile) = vns(ik,ie,is) &
                  * al(il)*aj1le(ixi,ie,ile) * gle(ixi,ie,ile)
          end do
       end do
    end do

    call integrate_moment (le_lo, gtmp, v2y2)    ! v2y2

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t6)
!       t6tot = t6tot + real(t6-t5)/tr
!    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle

       gle(:,:,ile) = gle(:,:,ile) - ieqzip(it,ik)*bw0le(:,:,ile) * v2y2(ile)
    end do

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       call system_clock (count=t7)
!       t7tot = t7tot + real(t7-t6)/tr
!    end if

    ! TMP FOR TESTING -- MAB
!    if (proc0) then
!       write (*,'(a20,7e14.5)') 'conserve_diffuse: ', t1tot, t2tot, t3tot, t4tot, t5tot, t6tot, t7tot
!    end if

    deallocate (vpadelnu, vns, v0y0, v1y1, v2y2, gtmp)

  end subroutine conserve_diffuse_le_layout

  subroutine solfp_lorentz_standard_layout (g, gc, gh, diagnostics)

    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, jend, ng2, lambda_map, nxi
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: ig_idx, ik_idx, il_idx, is_idx, it_idx, ie_idx
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    use gs2_layouts, only: lz_lo

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gc, gh
    integer, optional, intent (in) :: diagnostics

    complex, dimension (:,:), allocatable :: glz, glzc
    integer :: ilz
    complex, dimension (nxi+1) :: delta
    complex :: fac, gwfb
    integer :: ig, ik, il, is, it, je, ie
    integer :: nxi_scatt

    call prof_entering ("solfp_lorentz", "collisions")

    allocate (glz(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
    glz = 0.0
    if (heating) then
       allocate (glzc(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       glzc = 0.0
    end if

    call gather (lambda_map, g, glz)

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
       call scatter (lambda_map, glzc, gc)

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
          call scatter (lambda_map, glzc, gh)
       end if
    end if

    ! solve for glz row by row
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc

       ik = ik_idx(lz_lo,ilz)
       it = it_idx(lz_lo,ilz)
       if(kwork_filter(it,ik))cycle
       if (ieqzip(it,ik)==0) cycle
       is = is_idx(lz_lo,ilz)
       if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) cycle
       ie = ie_idx(lz_lo,ilz)
       ig = ig_idx(lz_lo,ilz)

!CMRDDGC, 10/2/2014: 
! Fixes for wfb treatment below, use same je definition in ALL cases
!   je  = #physical xi values + 1 
!       NB +1 above WITH TRAPPED is duplicate xi=vpar=0 point with isign=2
!          +1 above WITHOUT TRAPPED is entirely unphysical extra grid point
!  je-1 = #physical xi values removing unphysical/duplicate extra point
       je=max(2*jend(ig),2*ng2+1)
       nxi_scatt=je-1
       glz(je,ilz)=0.0  ! zero unphysical/duplicate extra point
       if (jend(ig) == ng2+1 .and.special_wfb_lorentz) then
!CMRDDGC:  special_wfb_lorentz = t  => unphysical handling of wfb at bounce pt: 
!          remove wfb from collisions, reinsert later
!
! first save gwfb for reinsertion later
          gwfb = glz(ng2+1,ilz) 
! then remove vpa = 0 point, weight 0: (CMR confused by this comment!)  
          glz(ng2+1:je-2,ilz) = glz(ng2+2:je-1,ilz)
          nxi_scatt=nxi_scatt-1
       endif
!CMRDDGCend

       ! right and left sweeps for tridiagonal solve:

       delta(1) = glz(1,ilz)
       do il = 1, nxi_scatt
          delta(il+1) = glz(il+1,ilz) - ql(il+1,ilz)*delta(il)
       end do
       
!CMR, 5/9/14:
! Correction to next line fixes a bug introduced at r2699 which
! affected Lorentz scattering operator at wfb bounce point (theta=pi)
! with default setting of special_wfb_lorentz=t        
       glz(nxi_scatt+1,ilz) = delta(nxi_scatt+1)*betaa(nxi_scatt+1,ilz)
       do il = nxi_scatt, 1, -1
          glz(il,ilz) = (delta(il) - c1(il,ilz)*glz(il+1,ilz))*betaa(il,ilz)
       end do

!       ! interpolate to obtain glz(vpa = 0) point for wfb
!       ! and insert this point into glz
       ! interpolation described above mysteriously causing numerical instability
       ! stabilized by using old (pre-collision) value of g for wfb
       if (jend(ig) == ng2+1.and.special_wfb_lorentz) then
          glz(ng2+2:je-1,ilz) = glz(ng2+1:je-2,ilz)
          glz(ng2+1,ilz) = gwfb
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

    call scatter (lambda_map, glz, g)

    deallocate (glz)
    if (heating) deallocate (glzc)

    call prof_leaving ("solfp_lorentz", "collisions")

  end subroutine solfp_lorentz_standard_layout


  subroutine solfp_lorentz_le_layout (gle, diagnostics)

    use le_grids, only: jend, ng2, negrid, nxi, integrate_moment
    use gs2_layouts, only: ig_idx, ik_idx, is_idx, it_idx
    use prof, only: prof_entering, prof_leaving
    use run_parameters, only: ieqzip
    use gs2_layouts, only: le_lo
    use dist_fn_arrays, only: c_rate
    use kt_grids, only: kwork_filter
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    integer, optional, intent (in) :: diagnostics

    complex, dimension (:), allocatable :: tmp
    integer :: ile
    complex, dimension (nxi+1) :: delta
    complex :: fac, gwfb
    integer :: ig, ik, il, is, je, it, ie
    integer :: nxi_scatt

    call prof_entering ("solfp_lorentz", "collisions")

    allocate (tmp(le_lo%llim_proc:le_lo%ulim_alloc)) ; tmp = 0.0

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
       call integrate_moment (le_lo, glec, tmp)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ig = ig_idx(le_lo,ile)
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          is = is_idx(le_lo,ile)
          c_rate(ig,it,ik,is,1) = tmp(ile)
       end do

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
          call integrate_moment (le_lo, glec, tmp)
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             ig = ig_idx(le_lo,ile)
             it = it_idx(le_lo,ile)
             ik = ik_idx(le_lo,ile)
             is = is_idx(le_lo,ile)
             c_rate(ig,it,ik,is,2) = tmp(ile)
          end do
       end if
    end if

    ! solve for gle row by row
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if (kwork_filter(it,ik)) cycle
       if (ieqzip(it,ik)==0) cycle
       is = is_idx(le_lo,ile)
       if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) cycle
       ig = ig_idx(le_lo,ile)

!CMRDDGC, 10/2/1014: 
! Fixes for wfb treatment below, use same je definition in ALL cases
!   je  = #physical xi values at location, includes duplicate point at vpar=0
!  je-1 = #physical xi values removing duplicate vpar=0 point
       je=max(2*jend(ig),2*ng2+1)
       nxi_scatt=je-1
       if (jend(ig) == ng2+1 .and.special_wfb_lorentz) then
          nxi_scatt=nxi_scatt-1
       endif
! These fixes address previous comment by CMR: 
! "surely wfb's bounce point should be handled like any other trapped bp?"
!CMR, 1/11/2013:
! Numerical instability referred to below more likely arises from wfb failing
! to satisfy trapping condition after invert_rhs.
!
       do ie = 1, negrid
          gle(je,ie,ile)=0.0d0  ! zero redundant duplicate xi, isign=2 for vpar=0!
          if (jend(ig) == ng2+1 .and.special_wfb_lorentz) then
!CMRDDGC:  special_wfb_lorentz = t  => unphysical handling of wfb at bounce pt: 
!          remove wfb from collisions, reinsert later
!
! first save gwfb for reinsertion later
             gwfb = gle(ng2+1,ie,ile)
! then remove vpa = 0 point, weight 0: (CMR confused by this comment!)  
             gle(ng2+1:je-2,ie,ile) = gle(ng2+2:je-1,ie,ile)
          endif
!CMRDDGCend

          ! right and left sweeps for tridiagonal solve:
          
          delta(1) = gle(1,ie,ile)
          do il = 1, nxi_scatt
             delta(il+1) = gle(il+1,ie,ile) - qle(il+1,ie,ile)*delta(il)
          end do
!CMR, 5/9/14:
! Correction to next line fixes a bug introduced at r2699 which
! affected Lorentz scattering operator at wfb bounce point (theta=pi)
! with default setting of special_wfb_lorentz=t        
          gle(nxi_scatt+1,ie,ile) = delta(nxi_scatt+1)*betaale(nxi_scatt+1,ie,ile)
          do il = nxi_scatt, 1, -1
             gle(il,ie,ile) = (delta(il) - c1le(il,ie,ile)*gle(il+1,ie,ile))*betaale(il,ie,ile)
          end do
!       ! interpolate to obtain glz(vpa = 0) point for wfb
!       ! and insert this point into glz
          ! interpolation described above mysteriously causing numerical instability
          ! stabilized by using old (pre-collision) value of g for wfb
          if (jend(ig) == ng2+1.and.special_wfb_lorentz) then
             gle(ng2+2:je-1,ie,ile)=gle(ng2+1:je-2,ie,ile)
             gle(ng2+1,ie,ile)=gwfb
!          gle(ng2+1,ie,ile) = 0.5*(gle(ng2,ie,ile)+gle(ng2+1,ie,ile))
          end if

!CMR, 1/11/2013:
! next line ensures bounce condition is satisfied after lorentz collisions
! this is right thing to do, but it would mask any prior bug in trapping condition!

          if (jend(ig) /= 0) gle(2*jend(ig),ie,ile) = gle(jend(ig),ie,ile)

       end do
    end do

    deallocate (tmp)
    
    call prof_leaving ("solfp_lorentz", "collisions")

  end subroutine solfp_lorentz_le_layout

  ! energy diffusion subroutine used with energy layout (not le_layout)
  ! this is always the case when initializing the conserving terms,
  ! otherwise is the case if use_le_layout is no specified in the input file.
  subroutine solfp_ediffuse_standard_layout (g)

    use species, only: spec
    use theta_grid, only: ntgrid
    use le_grids, only: negrid, forbid
    use le_grids, only: energy_map
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, il_idx, is_idx, e_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g

    integer :: ie, is, ig, il, it, ik
    complex, dimension (negrid) :: delta
    complex, dimension (:,:), allocatable :: ged
    integer :: ielo

    allocate (ged(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc)) ; ged = 0.0

    call gather (energy_map, g, ged)

    ! solve for ged row by row
    do ielo = e_lo%llim_proc, e_lo%ulim_proc

       it = it_idx(e_lo,ielo)       
       ik = ik_idx(e_lo,ielo)
       if(kwork_filter(it,ik))cycle
       if (ieqzip(it,ik)==0) cycle
       is = is_idx(e_lo,ielo)
       if (spec(is)%vnewk < 2.0*epsilon(0.0)) cycle
       ig = ig_idx(e_lo,ielo)
       il = il_idx(e_lo,ielo)
       if (forbid(ig,il)) cycle

       delta(1) = ged(1,ielo)
       do ie = 1, negrid-1
          delta(ie+1) = ged(ie+1,ielo) - eql(ie+1,ielo)*delta(ie)
       end do
       
       ged(negrid+1,ielo) = 0.0
       do ie = negrid, 1, -1
          ged(ie,ielo) = (delta(ie) - ec1(ie,ielo)*ged(ie+1,ielo))*ebetaa(ie,ielo)
       end do

    end do

    call scatter (energy_map, ged, g)

    deallocate (ged)

  end subroutine solfp_ediffuse_standard_layout


  subroutine solfp_ediffuse_le_layout (gle)

    use species, only: spec
    use le_grids, only: nxi, negrid, forbid, ixi_to_il
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx, le_lo
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle

    integer :: ie, is, ig, il
    complex, dimension (negrid) :: delta
    integer :: ile, ixi, ik, it

    ! solve for gle row by row
    do ile = le_lo%llim_proc, le_lo%ulim_proc

       is = is_idx(le_lo,ile)
       if (spec(is)%vnewk < 2.0*epsilon(0.0)) cycle
       it=it_idx(le_lo,ile)
       ik=ik_idx(le_lo,ile)
       if (kwork_filter(it,ik)) cycle
       if (ieqzip(it,ik)==0) cycle
       ig = ig_idx(le_lo,ile)

       do ixi = 1, nxi
          il = ixi_to_il(ig,ixi)
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

  end subroutine solfp_ediffuse_le_layout


  subroutine check_g (str, g)

    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    implicit none    
    character (3), intent(in) :: str 
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
    implicit none
    character (3), intent(in) :: str 
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

    if(use_le_layout) then
      if (allocated(c1le)) then
         deallocate (c1le, betaale, qle)
         if (heating) deallocate (d1le, h1le)
      end if
      if (allocated(ec1le)) deallocate (ec1le, ebetaale, eqle)
      if (allocated(glec)) deallocate (glec)
      if (allocated(vpatmp)) deallocate (vpatmp)
      if (allocated(s0le)) deallocate(s0le)
      if (allocated(w0le)) deallocate(w0le)
      if (allocated(z0le)) deallocate(z0le)
      if (allocated(aj0le)) deallocate(aj0le)
      if (allocated(aj1le)) deallocate(aj1le)
      if (allocated(bs0le)) deallocate(bs0le)
      if (allocated(bw0le)) deallocate(bw0le)
      if (allocated(bz0le)) deallocate(bz0le)
    else
      if (allocated(c1)) then
         deallocate (c1, betaa, ql)
         if (heating) deallocate (d1, h1)
      end if
      if (allocated(ec1)) deallocate (ec1, ebetaa, eql)
    end if
    if (allocated(vpdiff)) deallocate (vpdiff)
    if (allocated(dtot)) deallocate (dtot, fdf, fdb)

  end subroutine finish_collisions

end module collisions

