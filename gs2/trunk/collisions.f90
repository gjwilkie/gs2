module collisions
  
  use redistribute, only: redist_type

  implicit none

  public :: init_collisions
  public :: solfp1
  public :: reset_init
  public :: dtot, fdf, fdb, lorentz_map, vnmult, vnfac
  public :: ncheck, vnslow, vary_vnew
  public :: etol, ewindow, etola, ewindowa
  public :: init_lorentz, init_ediffuse
  public :: init_lorentz_conserve, init_diffuse_conserve
  public :: init_lorentz_error, collision_model_switch
  public :: ntot_diff, upar_diff, uperp_diff, ttot_diff

  private

  ! knobs
  real :: vncoef, absom
  integer :: ivnew
  logical :: conserve_number, conserve_momentum, const_v, conserve_moments
  logical :: test_mom_conserve, new_coll_scheme
  integer :: collision_model_switch
  logical :: use_shmem, adjust
  logical :: heating
  logical :: hyper_colls
  logical :: diffuse_energy

  integer, parameter :: collision_model_lorentz = 1      ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_krook = 2
  integer, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_krook_test = 4
  integer, parameter :: collision_model_lorentz_test = 5 ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_full = 6
  integer, parameter :: collision_model_ediffuse = 7

  real, dimension (2), save :: vnmult = 0.0
  integer :: ncheck
  logical :: vary_vnew
  real :: vnfac, vnslow
  real :: etol, ewindow, etola, ewindowa

  real, dimension (:,:,:), allocatable :: dtot
  ! (-ntgrid:ntgrid,nlambda,max(ng2,nlambda-ng2)) lagrange coefficients for derivative error estimate

  real, dimension (:,:), allocatable :: fdf, fdb
  ! (-ntgrid,ntgrid,nlambda) finite difference coefficients for derivative error estimate

  real, dimension (:,:,:), allocatable :: vnew, vnew_s, vnew_D, vnew_E
  ! (naky,negrid,nspec) replicated

  ! only for hyper-diffusive collisions
  real, dimension (:,:,:,:), allocatable :: vnewh
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec) replicated

  ! TEMPORARY FOR TESTING -- MAB
  complex, dimension (:,:,:,:), allocatable :: ntot_diff, upar_diff, uperp_diff, ttot_diff
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec) replicated

  ! only for krook
  real, dimension (:), allocatable :: vnewfe
  ! (-*- g_layout -*-)

  ! only for krook
  real, dimension (:,:), allocatable :: aintnorm
  ! (-ntgrid:ntgrid, -*- geint_layout -*-)

  ! only for momentum conservation due to Lorentz operator (8.06)
  complex, dimension(:,:,:), allocatable :: z0, z1

  ! needed for momentum and energy conservation due to energy diffusion (3.08)
  complex, dimension(:,:,:), allocatable :: bs0, bw0, bz0
  complex, dimension(:,:,:,:), allocatable :: v2w1

  ! only for original parallel mom conservation (not used nowadays)
  real, dimension (:,:,:), allocatable :: sq
  ! (-ntgrid:ntgrid,nlambda,2) replicated

  ! only for energy diffusion
  real, dimension (:,:), allocatable :: ec1, ebetaa, eql
  real, dimension (:,:), allocatable :: era1, erb1, erc1

! only around for testing -- MAB
!  real, dimension (:), allocatable :: ea1, eb1

  ! only for lorentz
  real :: cfac
  real, dimension (:,:), allocatable :: c1, betaa, ql, d1, h1
  real, dimension (:,:), allocatable :: ra1, rb1, rc1

! not necessary to make these variables global (only used in solfp_lorentz)
! MAB
!  complex, dimension (:,:), allocatable :: glz, glzc
  ! (max(2*nlambda,2*ng2+1), -*- lz_layout -*-)
  ! (-ntgrid:ntgrid,naky,max(2*nlambda,2*ng2+1),negrid,nspec)

  ! momentum conservation
  complex, dimension (:,:), allocatable :: g3int
  ! (-ntgrid:ntgrid, -*- gint_layout -*-)

  type (redist_type), save :: lorentz_map
  type (redist_type), save :: ediffuse_map

  logical :: hypermult
  logical :: initialized = .false.
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
    call init_arrays
  end subroutine init_collisions

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options
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
    character(20) :: collision_model
    namelist /collisions_knobs/ collision_model, vncoef, absom, ivnew, &
         conserve_number, conserve_momentum, use_shmem, heating, &
         adjust, const_v, cfac, hypermult, diffuse_energy, vnfac, &
         etol, ewindow, ncheck, vnslow, vary_vnew, etola, ewindowa, &
         test_mom_conserve, conserve_moments, new_coll_scheme
    integer :: ierr, in_file
    logical :: exist

    if (proc0) then
       hypermult = .false.
       cfac = 1.   ! DEFAULT CHANGED TO INCLUDE CLASSICAL DIFFUSION: APRIL 18, 2006
       vnfac = 1.1
       vnslow = 0.9
       etol = 2.e-2
       ewindow = 1.e-2
       etola = 2.e-2
       ewindowa = 1.e-2
       ncheck = 100
       adjust = .true.
       collision_model = 'default'
       vncoef = 0.6
       absom = 0.5
       ivnew = 0
       conserve_number = .true.
       conserve_momentum = .true.  ! DEFAULT CHANGED TO REFLECT IMPROVED MOMENTUM CONSERVATION, 8/06
       conserve_moments = .false.
       new_coll_scheme = .true.
       test_mom_conserve = .false.
       diffuse_energy = .false.
       vary_vnew = .false.
       const_v = .false.
       heating = .false.
       in_file = input_unit_exist ("collisions_knobs", exist)
       if (exist) read (unit=input_unit("collisions_knobs"), nml=collisions_knobs)

       ierr = error_unit()
       call get_option_value &
            (collision_model, modelopts, collision_model_switch, &
            ierr, "collision_model in collisions_knobs")
    end if

    call broadcast (hypermult)
    call broadcast (cfac)
    call broadcast (vnfac)
    call broadcast (vnslow)
    call broadcast (vary_vnew)
    call broadcast (etol)
    call broadcast (ewindow)
    call broadcast (etola)
    call broadcast (ewindowa)
    call broadcast (ncheck)
    call broadcast (vncoef)
    call broadcast (absom)
    call broadcast (ivnew)
    call broadcast (conserve_number)
    call broadcast (conserve_momentum)
    call broadcast (conserve_moments)
    call broadcast (new_coll_scheme)
    call broadcast (test_mom_conserve)
    call broadcast (diffuse_energy)
    call broadcast (const_v)
    call broadcast (collision_model_switch)
    call broadcast (heating)
    call broadcast (adjust)

  end subroutine read_parameters

  subroutine init_arrays
    use species, only: nspec
    use le_grids, only: negrid
    use gs2_layouts, only: g_lo
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: c_rate
    implicit none
    real, dimension (negrid,nspec) :: hee
    logical :: first_time = .true.

    if (first_time) then
       allocate (c_rate(-ntgrid:ntgrid, ntheta0, naky, nspec, 3))
       c_rate = 0.
       first_time = .false.
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

!
! Precompute two quantities needed for momentum conservation:
! z0, z1
    
    use mp, only: proc0
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: e, al, integrate_moment, negrid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, kperp2, vpa

    logical, save :: first = .true.
    complex, dimension (1,1,1) :: dum1 = 0., dum2 = 0.
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: v0z0, v0z1, v1z0, v1z1, duinv
    real, dimension (:,:,:), allocatable :: vns
!    real :: vnm1, vnm2
    integer :: ie, il, ik, is, ig, isgn, iglo, all, it

!    vnm1 = vnmult(1)
!    vnm2 = vnmult(2)

! TO DO: 
! tunits not included anywhere yet

    if (first) then
       allocate (z0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (z1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       first = .false.
    end if

! First, get du and then 1/du == duinv

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))       
    allocate (vns(naky,negrid,nspec))

!    vns = (vnm1-vnm2)*vnew_D + vnm2*vnew_s
    vns = vnmult(1)*vnew_D

!
! du == int (E nu_s f_0);  du = du(z, kx, ky, s)
! duinv = 1/du
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             gtmp(ig,isgn,iglo) = e(ie,is)*vns(ik,ie,is)
          end do
       end do
    end do

    all = 1
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet

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
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! recall that aj1 == J1 (arg) / arg, arg = sqrt(kperp2*e*al*(T*m/B/q**2))
! u0 = -3 nu_s V_perp dt a f_0 / du
! where a = kperp2 * (T m / q**2) / b(theta)**2
!>MAB
! note that u0 here is defined slightly differently from u0 in write-up.
! in particular, the u0 from write-up has only factor of sqrt(a)...
! the additional factor of sqrt(a) appearing here is taken from v0
!<MAB
             z0(ig,isgn,iglo) = - 3.*vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * code_dt * spec(is)%smz**2 * kperp2(ig,it,ik) * duinv(ig,it,ik,is) &
                  / bmag(ig)
          end do
       end do
    end do

    call solfp_lorentz (z0, dum1, dum2)   ! z0 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get z1 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa 
! V_parallel == v_parallel J0
! u1 = -3 nu V_parallel dt f_0 / du
!
! No factor of sqrt(T/m) here on purpose (see derivation) 
!
             z1(ig,isgn,iglo) = - 3.*vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * code_dt * duinv(ig,it,ik,is)
          end do
       end do
    end do

    deallocate (duinv)  ! Done with this variable

    call solfp_lorentz (z1, dum1, dum2)    ! z1 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    allocate (v0z0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! v0 = nu V_perp
! no sqrt(kperp2 smz**2 / bmag) because already absorbed into z0 (u0 from write-up)
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * z0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v0z0, all)    ! v0z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1z0

    allocate (v1z0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa
! V_parallel == v_parallel J0
! v1 = nu V_parallel f_0
!
! No factor of sqrt(T/m) here on purpose (see derivation) 
!
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * z0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v1z0, all)    ! v1z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z1

    allocate (v0z1(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! v0 = nu V_perp
! no sqrt(kperp2 smz**2 / bmag) in v0 because it was absorbed into u0 -- MAB 
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * z1(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v0z1, all)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now redefine z1 == z1 - z0 [v0 . z1]/(1+v0.z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             z1(ig,isgn,iglo) = z1(ig,isgn,iglo) - z0(ig,isgn,iglo)*v0z1(ig,it,ik,is) &
                  / (1.+v0z0(ig,it,ik,is))
          end do
       end do
    end do

    deallocate (v0z1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1z1

    allocate (v1z1(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa
! V_parallel == v_parallel J0
! v1 = nu V_parallel f_0
!
! No factor of sqrt(T/m) here on purpose (see derivation) 
!
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * z1(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v1z1, all)    ! redefined below

    deallocate (gtmp, vns)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now redefine z1 == z1/(1 + v1z1)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             z1(ig,isgn,iglo) = z1(ig,isgn,iglo) / (1.+v1z1(ig,it,ik,is))
          end do
       end do
    end do

    deallocate (v1z1)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now redefine z0 == (z1 * v1z0 - z0) / (1 + v0z0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             z0(ig,isgn,iglo) = (z1(ig,isgn,iglo) * v1z0(ig,it,ik,is) - z0(ig,isgn,iglo)) &
                  / (1.+v0z0(ig,it,ik,is))
          end do
       end do
    end do

    deallocate (v0z0, v1z0)
    
  end subroutine init_lorentz_conserve

  subroutine init_diffuse_conserve

!
! Precompute four quantities needed for momentum and energy conservation:
! z0, w0, s0, v2w1
    
    use mp, only: proc0
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: e, al, integrate_moment, negrid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, kperp2, vpa

    logical, save :: first = .true.
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: v0z0, duinv
    real, dimension (:,:,:), allocatable :: vns
!    real :: vnm1, vnm2
    integer :: ie, il, ik, is, ig, isgn, iglo, all, it

!    vnm1 = vnmult(1)
!    vnm2 = vnmult(2)

! TO DO: 
! tunits not included anywhere yet

    if (first) then
       allocate (bz0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (bw0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (bs0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (v2w1(-ntgrid:ntgrid,ntheta0,naky,nspec))
       first = .false.
    end if

! First, get du and then 1/du == duinv

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))       
    allocate (vns(naky,negrid,nspec))

!    vns = (vnm1-vnm2)*vnew_D + vnm2*vnew_s
    vns = vnmult(2)*(vnew_s - vnew_D)

!
! du == int (E nu_s f_0);  du = du(z, kx, ky, s)
! duinv = 1/du
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             gtmp(ig,isgn,iglo) = e(ie,is)*vns(ik,ie,is)
          end do
       end do
    end do

    all = 1
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet

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
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! recall that aj1 == J1 (arg) / arg, arg = sqrt(kperp2*e*al*(T*m/B/q**2))
! u0 = -3 nu_s V_perp dt a f_0 / du
! where a = kperp2 * (T m / q**2) / b(theta)**2
!>MAB
! note that u0 here is defined slightly differently from u0 in write-up.
! in particular, the u0 from write-up has only factor of sqrt(a)...
! the additional factor of sqrt(a) appearing here is taken from v0
!<MAB
             bs0(ig,isgn,iglo) = - 3.*vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * code_dt * spec(is)%smz**2 * kperp2(ig,it,ik) * duinv(ig,it,ik,is) &
                  / bmag(ig)
          end do
       end do
    end do

    call solfp_ediffuse (bs0)   ! s0 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get w0 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa 
! V_parallel == v_parallel J0
! u1 = -3 nu V_parallel dt f_0 / du
!
             bw0(ig,isgn,iglo) = - 3.*vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * code_dt * duinv(ig,it,ik,is)
          end do
       end do
    end do

    call solfp_ediffuse (bw0)    ! w0 is redefined below

! redefine du = dq (for energy-conserving term)
! du == int (E^2 nu_E f_0);  du = du(z, kx, ky, s)
! duinv = 1/du
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             gtmp(ig,isgn,iglo)  = e(ie,is)**2*vnmult(2)*vnew_E(ik,ie,is)
          end do
       end do
    end do

    all = 1
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet

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
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! u2 = -nu_E v^2 dt J0 f_0 / du
             bz0(ig,isgn,iglo) = - vnmult(2)*vnew_E(ik,ie,is)*e(ie,is)*aj0(ig,iglo) &
                  * code_dt * duinv(ig,it,ik,is)
          end do
       end do
    end do

    deallocate (duinv)  ! Done with this variable

    call solfp_ediffuse (bz0)    ! z0 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2w1 (will redefine later...really getting v0s0 right now)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! v0 = nu V_perp
! no sqrt(kperp2 smz**2 / bmag) because already absorbed into bs0 (u0 from write-up)
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * bs0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v2w1, all)    ! v0s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 / (1 + v0s0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             bs0(ig,isgn,iglo) = bs0(ig,isgn,iglo) / (1.0 + v2w1(ig,it,ik,is))
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine v2w1 to actually be v0w0 (instead of v0s0)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! v0 = nu V_perp
! no sqrt(kperp2 smz**2 / bmag) because already absorbed into s0 (u0 from write-up)
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * bw0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v2w1, all)    ! v0w0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w0 - v0w0 * s0 (this is w1 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             bw0(ig,isgn,iglo) = bw0(ig,isgn,iglo) - bs0(ig,isgn,iglo)*v2w1(ig,it,ik,is)
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)*bmag(ig)
! v0 = nu V_perp
! no sqrt(kperp2 smz**2 / bmag) because already absorbed into s0 (u0 from write-up)
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * bz0(ig,isgn,iglo)
          end do
       end do
    end do

    allocate (v0z0(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (gtmp, v0z0, all)    ! v0z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine z0 = z0 - v0z0 * s0 (this is z1 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             bz0(ig,isgn,iglo) = bz0(ig,isgn,iglo) - bs0(ig,isgn,iglo)*v0z0(ig,it,ik,is)
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine v2w1 = v1 . w1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa
! V_parallel == v_parallel J0
! v1 = nu V_parallel f_0
!
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * bw0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v2w1, all)    ! v1w1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w1 / (1 + v1 . w1)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             bw0(ig,isgn,iglo) = bw0(ig,isgn,iglo) / (1.0 + v2w1(ig,it,ik,is))
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine v0z0 = v1 . z1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa
! V_parallel == v_parallel J0
! v1 = nu V_parallel f_0
!
             gtmp(ig,isgn,iglo) = vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * bz0(ig,isgn,iglo)
          end do
       end do
    end do

    deallocate (vns)

    call integrate_moment (gtmp, v0z0, all)    ! v1z1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine z0 = z1 - v1z1 * w1 (this is z2 from MAB notes)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             bz0(ig,isgn,iglo) = bz0(ig,isgn,iglo) - bw0(ig,isgn,iglo)*v0z0(ig,it,ik,is)
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine v2w1 = v2 . w1

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v2 = nu_E*e*J0 
             gtmp(ig,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*e(ie,is)*aj0(ig,iglo) &
                  * bw0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v2w1, all)   ! v2w1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine v0z0 = v2 . z2

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v2 = nu_E*e*J0 
             gtmp(ig,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*e(ie,is)*aj0(ig,iglo) &
                  * bz0(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v0z0, all)   ! v2z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now redefine z0 = z2 / (1 + v2z2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             bz0(ig,isgn,iglo) = bz0(ig,isgn,iglo) / (1.0 + v0z0(ig,it,ik,is))
          end do
       end do
    end do

    deallocate (v0z0)
    
  end subroutine init_diffuse_conserve

  subroutine init_vnew (hee)
    use species, only: nspec, spec, electron_species, has_electron_species
    use le_grids, only: negrid, e
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use run_parameters, only: zeff, tunits
    use dist_fn_arrays, only: kperp2
    use constants
    real, dimension (:,:), intent (out) :: hee
    real,dimension (negrid,nspec)::heevth, hsg, hsgvth
    integer :: ik, ie, is, it, ig
    real :: v, k4max

    do is = 1, nspec
       do ie = 1, negrid
          v = sqrt(e(ie,is))
          hee(ie,is) = 1.0/sqrt(pi)/v*exp(-e(ie,is)) &  ! hee is vnew_D from MAB notes
               + (1.0 - 0.5/e(ie,is)) &
               *(1.0 - 1.0/(1.0          + v &          ! this line on is erf(v)
               *(0.0705230784 + v &
               *(0.0422820123 + v &
               *(0.0092705272 + v &
               *(0.0001520143 + v &
               *(0.0002765672 + v &
               *(0.0000430638)))))))**16)

!>MAB
! hsg is the G of Hirshman and Sigmar
! added to allow for momentum conservation with energy diffusion
          hsg(ie,is) = -1.0/sqrt(pi)/v*exp(-e(ie,is)) &
               + 0.5/e(ie,is) &
               *(1.0 - 1.0/(1.0 + v &  ! this line on is erf(v)
               *(0.0705230784 + v &
               *(0.0422820123 + v &
               *(0.0092705272 + v &
               *(0.0001520143 + v &
               *(0.0002765672 + v &
               *(0.0000430638)))))))**16)
!<MAB
       end do
    end do

!heevth is hee but using only thermal velocity (energy independent)

    do is = 1, nspec
       do ie = 1, negrid
          v = 1           
          heevth(ie,is) = 1.0/sqrt(pi)/v*exp(-v**2) &
               + (1.0 - 0.5/v**2) &
               *(1.0 - 1.0/(1.0          + v &
               *(0.0705230784 + v &
               *(0.0422820123 + v &
               *(0.0092705272 + v &
               *(0.0001520143 + v &
               *(0.0002765672 + v &
               *(0.0000430638)))))))**16)

!>MAB
! hsg is the G of Helander and Sigmar
! added to allow for momentum conservation with energy diffusion
          hsgvth(ie,is) = -1.0/sqrt(pi)/v*exp(-v**2) &
               + 0.5/v**2 &
               *(1.0 - 1.0/(1.0 + v &
               *(0.0705230784 + v &
               *(0.0422820123 + v &
               *(0.0092705272 + v &
               *(0.0001520143 + v &
               *(0.0002765672 + v &
               *(0.0000430638)))))))**16)
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
                   vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                        - 2.0*vnew_D(ik,ie,is)
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *(zeff + hee(ie,is))*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(e(ie,is)) &
                        *hsg(ie,is)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *hee(ie,is)*tunits(ik)
                   vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*(2.0-0.5/e(ie,is)) &
                        - 2.0*vnew_D(ik,ie,is)
                end if
             end do
          end do
       else
          do ie = 1, negrid
             do ik = 1, naky
                if (const_v) then
                   vnew(ik,ie,is) = spec(is)%vnewk &
                        *heevth(ie,is)*0.5*tunits(ik)
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *hee(ie,is)*0.5*tunits(ik)
                end if
                vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(e(ie,is)) &
                     *hsg(ie,is)*4.0*tunits(ik)
                vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*(2.0-0.5/e(ie,is)) &
                     - 2.0*vnew_D(ik,ie,is)
             end do
          end do
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
    use mp, only: proc0
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, ng2, e, ecut, integrate_moment, al, w
    use egrid, only: zeroes, x0, energy
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: init_ediffuse_layouts
    use gs2_layouts, only: e_lo, g_lo, ie_idx, ik_idx, is_idx
    use gs2_layouts, only: ig_idx, it_idx, il_idx
    use file_utils, only: open_output_file, close_output_file
    use dist_fn_arrays, only: kperp2

    implicit none
    
    real, intent (in), optional :: vnmult_target

    integer :: ie, is, iglo, ik, ielo, il, ig, it
    real, dimension (:), allocatable :: aa, bb, cc, xe, ba, rba, eec
    real :: vn, xe0, xe1, xe2, xer, xel, er, el, fac
    real :: dela, delb, delc, delfac
    real :: capgl, capgr, slb1, ee, eea, eeb
    real :: delp, delm, del, capg, mu, eta, wxe, capgp

    logical :: first_time = .true.

    if (first_time) then
       vnmult(2) = max(1.0, vnmult(2))
       first_time = .false.
    end if

    call init_ediffuse_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)
    call init_ediffuse_redistribute

! no point in making these global variables and allocating here (only used in solfp_ediffuse)
! MAB
!    if (.not.allocated(ged)) then
!       allocate (ged(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc))
!       allocate (ged0(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc))
!       ged = 0.0 ; ged0 = 0.0
!    end if

    allocate (aa(negrid), bb(negrid), cc(negrid))
    allocate (ba(negrid), rba(negrid))
    allocate (xe(negrid), eec(negrid))

! want to use x variables instead of e because we want conservative form
! for the x-integration
    xe(1:negrid-1) = zeroes
    xe(negrid) = x0

    if (.not.allocated(ec1)) then
       allocate (ec1   (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (era1   (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (erb1   (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (erc1   (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (ebetaa(negrid,e_lo%llim_proc:e_lo%ulim_alloc))
       allocate (eql   (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
    endif

    if (present(vnmult_target)) then
       vnmult(2) = max (vnmult_target, 1.0)
    end if

    ec1 = 0.0 ; ebetaa = 0.0 ; eql = 0.0
    era1 = 0.0 ; erb1 = 0.0 ; erc1 = 0.0

    do ielo = e_lo%llim_proc, e_lo%ulim_proc
       is = is_idx(e_lo, ielo)
       ik = ik_idx(e_lo, ielo)
       il = il_idx(e_lo, ielo)
       ig = ig_idx(e_lo, ielo)
       it = it_idx(e_lo, ielo)

       vn = vnmult(2)*spec(is)%vnewk*tunits(ik)

       slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))     ! xi_j

       ee =cfac*0.25*(1-slb1**2)*kperp2(ig,it,ik)*vn*code_dt &
            / (bmag(ig)*spec(is)%zstm)**2

       ! boundary at xe = 0
       wxe = 4.0*w(1,is)
       xe1 = xe(1)
       xe2 = xe(2)

!       ee = 0.25/e(1,is)**1.5*(1-slb1**2)*xe1 &
!            / (bmag(ig)*spec(is)%zstm)**2 &
!            * kperp2(ig,it,ik)*cfac

       delp = xe(2) - wxe
       delm = xe(1) - wxe
       del  = xe(2) - xe(1)

       capg = 8.0*xe1*sqrt(e(1,is))*exp(-2.0*e(1,is))/pi
       capgp = 8.0*exp(-e(1,is))*(sqrt(e(1,is))*exp(-e(1,is)) &
            + 0.25*sqrt(pi)*xe1*(1.-4.*e(1,is))/e(1,is))/pi

       ! TEMP FOR TESTING -- MAB
!       capg = 1.0 ; capgp = 0.0 ; vn = 1.0

       mu = 0.5*(delp + delm)
       eta = (1.0 - mu*capgp/capg)*del

       eec(1) = -ee*wxe*delm*mu/(eta*energy(wxe,ecut)**1.5)
       eeb = -ee*wxe*delp*mu/(eta*energy(wxe,ecut)**1.5)

       erc1(1,ielo) = -delm*mu/eta
       cc(1) = -(code_dt*vn*capg + delm*mu)/eta
       
       rba(1) = -delp*mu/eta
       erb1(1,ielo) = 1.0 - 0.25*rba(1)/w(1,is)
       ba(1) = -(code_dt*vn*capg + delp*mu)/eta
       bb(1) = 1.0 - 0.25*(ba(1) + eeb)/w(1,is) + ee*xe1/e(1,is)**1.5
       
       era1(1,ielo) = 0.0
       aa(1) = 0.0

! non-conservative scheme
!       xer = (xe2 + xe1)*0.5
!       
!       er = energy(xer,ecut)
!
!       fac = -8.0*vn*code_dt/pi
!
!       ee = 0.25/e(1,is)**1.5*(1-slb1**2)*xe1 &
!            / (bmag(ig)*spec(is)%zstm)**2 &
!            * kperp2(ig,it,ik)*cfac
!
!       aa(1) = 0.0
!       cc(1) = fac*sqrt(er)*exp(-2.0*er)/(xe2 - xe1)
!       bb(1) = 1.0 - cc(1) + ee*vn*code_dt

       do ie = 2, negrid-1

          eea = ee*wxe*delp*mu/(eta*energy(wxe,ecut)**1.5)

          wxe = 4.0*sum(w(:ie,is))
          xe1 = xe(ie)
          xe2 = xe(ie+1)

          delp = xe2 - wxe
          delm = xe1 - wxe
          del  = xe2 - xe1

          capg = 8.0*xe1*sqrt(e(ie,is))*exp(-2.0*e(ie,is))/pi
          capgp = 8.0*exp(-e(ie,is))*(sqrt(e(ie,is))*exp(-e(ie,is)) &
               + 0.25*sqrt(pi)*xe1*(1.-4.*e(ie,is))/e(ie,is))/pi

          ! TEMP FOR TESTING -- MAB
!          capg = 1.0 ; capgp = 0.0 ; vn = 1.0

          mu = 0.5*(delp + delm)
          eta = (1.0 - mu*capgp/capg)*del

          eeb = -ee*wxe*delp*mu/(eta*energy(wxe,ecut)**1.5) - eec(ie-1)
          eec(ie) = -ee*wxe*delm*mu/(eta*energy(wxe,ecut)**1.5)

          erc1(ie,ielo) = -delm*mu/eta
          cc(ie) = -(code_dt*vn*capg + delm*mu)/eta

          rba(ie) = -delp*mu/eta
          erb1(ie,ielo) = 1.0 - 0.25*(rba(ie) + erc1(ie-1,ielo))/w(ie,is)
          ba(ie) = -(code_dt*vn*capg + delp*mu)/eta
          bb(ie) = 1.0 - 0.25*(cc(ie-1) + ba(ie) + eeb)/w(ie,is) + ee*xe1/e(ie,is)**1.5

          era1(ie,ielo) = 0.25*rba(ie-1)/w(ie,is)
          aa(ie) = 0.25*(ba(ie-1) + eea)/w(ie,is)

! old non-conservative scheme
!          xe0 = xe(ie-1)
!          xe1 = xe(ie)
!          xe2 = xe(ie+1)
!
!          capgl = 0.5*(xe0*sqrt(e(ie-1,is))*exp(-2.0*e(ie-1,is)) &
!               +xe1*sqrt(e(ie,is))*exp(-2.0*e(ie,is)))
!          capgr = 0.5*(xe2*sqrt(e(ie+1,is))*exp(-2.0*e(ie+1,is)) &
!               +xe1*sqrt(e(ie,is))*exp(-2.0*e(ie,is)))
!
!          dela = xe1 - xe0
!          delb = xe2 - xe0
!          delc = xe2 - xe1
!          delfac = delc/dela - dela/delc
!       
!          fac = -16.0*vn*code_dt/pi
!
!          ee = 0.25/e(ie,is)**1.5*(1-slb1**2)*xe1 &
!               / (bmag(ig)*spec(is)%zstm)**2 &
!               * kperp2(ig,it,ik)*cfac
!
!          aa(ie) = fac*delc/(delb*dela)*(capgl/dela &
!               - xe1*sqrt(e(ie,is))*exp(-2.0*e(ie,is))*delfac/delb)
!          cc(ie) = fac*dela/(delb*delc)*(capgr/delc &
!               + xe1*sqrt(e(ie,is))*exp(-2.0*e(ie,is))*delfac/delb)
!          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*vn*code_dt

       end do

       ! boundary at xe = 1

       ie = negrid

       eea = ee*wxe*delp*mu/(eta*energy(wxe,ecut)**1.5)
       
       xe1 = xe(ie)
 
       eeb = -eec(ie-1)
       eec(ie) = 0.0

       erc1(ie,ielo) = 0.0
       cc(ie) = 0.0
          
       erb1(ie,ielo) = 1.0 - 0.25*erc1(ie-1,ielo)/w(ie,is)
       bb(ie) = 1.0 - 0.25*(cc(ie-1) + eeb)/w(ie,is) + ee*xe1/e(ie,is)**1.5
       
       era1(ie,ielo) = 0.25*rba(ie-1)/w(ie,is)
       aa(ie) = 0.25*(ba(ie-1) + eea)/w(ie,is)

!       xe0 = xe(negrid-1)
!       xe1 = xe(negrid)
!       xe2 = 1.0
!
!       xel = (xe1 + xe0)*0.5
!
!       el = energy(xel,ecut)
!
!       fac = -8.0*vn*code_dt/pi/(xe2 - xel)
!
!       ee = 0.25/e(negrid,is)**1.5*(1-slb1**2)*xe1 &
!            / (bmag(ig)*spec(is)%zstm)**2 &
!            * kperp2(ig,it,ik)*cfac
!
!       aa(negrid) = fac*sqrt(el)*xel*exp(-2.0*el)/(xe1-xe0)
!       cc(negrid) = 0.0
!       bb(negrid) = 1.0 - aa(negrid) + ee*vn*code_dt

! fill in the arrays for the tridiagonal
       erc1(:,ielo) = 0.25*erc1(:,ielo)/w(:,is)
       ec1(:,ielo) = 0.25*(cc + eec)/w(:,is)
       ebetaa(1,ielo) = 1.0/bb(1)
       do ie = 1, negrid-1
          eql(ie+1,ielo) = aa(ie+1)*ebetaa(ie,ielo)
          ebetaa(ie+1,ielo) = 1.0/(bb(ie+1)-eql(ie+1,ielo)*ec1(ie,ielo))
       end do

!       do ie = 1, negrid
!          write (*,*) aa(ie), bb(ie), ec1(ie,ielo), era1(ie,ielo), erb1(ie,ielo), erc1(ie,ielo), ielo, ie
!       end do

    end do

    deallocate(aa, bb, cc, ba, rba, xe, eec)

  end subroutine init_ediffuse

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
    logical :: done = .false.

    if (done) return

    call init_ediffuse_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
! ik = ik_idx, etc. should be moved outside of isign and ig loops here and later
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             il = il_idx(g_lo,iglo)
             is = is_idx(g_lo,iglo)
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
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             il = il_idx(g_lo,iglo)
             is = is_idx(g_lo,iglo)
             ie = ie_idx(g_lo,iglo)
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

    done = .true.

  end subroutine init_ediffuse_redistribute

  subroutine init_lorentz (vnmult_target)
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, al, jend, ng2,e, wl
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: kperp2
    use gs2_layouts, only: init_lorentz_layouts
    use gs2_layouts, only: lz_lo
    use gs2_layouts, only: ig_idx, ik_idx, ie_idx, is_idx, it_idx
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0

    implicit none

    real, intent (in), optional :: vnmult_target

    integer :: ig, il, ilz, it, ik, ie, is, je, te, te2, teh
    real, dimension (nlambda+1) :: aa, bb, cc
    real, dimension (nlambda+1) :: ra, rba, rb, rc, ba
    real, dimension (max(2*nlambda,2*ng2+1)) :: a1, b1
    real, dimension (nlambda+1) :: eec
    real, dimension (:), allocatable :: dd, hh
    real :: slb0, slb1, slb2, slbl, slbr, vn, ee, vnh, vnc
    real :: wslb, eta, mu, capg, delp, delm, del, wll
    real :: eea, eeb
    logical :: first_time = .true.

    if (heating) allocate (dd(nlambda+1), hh(nlambda+1))

    if (first_time) then
       vnmult(1) = max(1.0, vnmult(1))
       first_time = .false.
    end if

    call init_lorentz_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    call init_lorentz_redistribute

! no point in making these global variables and allocating here (only used in solfp_lorentz)
! MAB
!    if (.not.allocated(glz)) then
!       allocate (glz(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
!       glz = 0.0
!       if (heating) then
!          allocate (glzc(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
!          glzc = 0.0
!       end if
!    end if

    if (.not.allocated(c1)) then
       allocate (c1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (ra1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (rb1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (rc1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (betaa(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (ql   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
! changed to reduce memory storage -- MAB
!       allocate (d1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
!       allocate (h1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       if (heating) then
          allocate (d1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
          allocate (h1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
          d1 = 0.0 ; h1 = 0.0
       end if
    endif

    c1 = 0.0 ; betaa = 0.0 ; ql = 0.0
    ra1 = 0.0; rb1 = 0.0; rc1 = 0.0

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
          te = ng2
          te2 = 2*ng2
          teh = ng2
       else
          te = je 
          te2 = 2*je - 1 
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

       ! TEMP FOR TESTING -- MAB
!       vn = 1.0

       if (new_coll_scheme) then

          ee = cfac*0.5*e(ie,is)*kperp2(ig,it,ik)*vn*code_dt &
               / (bmag(ig)*spec(is)%zstm)**2

          il = 1
          wslb = 1.0 - 0.5*sum(wl(ig,:il))
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
          slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1)))

          delp = wslb - slb2
          delm = wslb - slb1
          del  = slb1 - slb2

          capg = 1.0 - wslb**2

          mu = 0.5*(delp + delm)
          eta = (1.0 + wslb*(delp + delm)/capg)*del

          eec(il) = -ee*(1.0 + wslb**2)*delm*mu/eta
!          eeb = -ee*((1.0 + wslb**2)*delp*mu/eta + (1+slb1**2))
          eeb = -ee*(1.0 + wslb**2)*delp*mu/eta

          rc(il) = -delm*mu/eta
          cc(il) = -(code_dt*vn*capg + delm*mu)/eta
       
          rba(il) = -delp*mu/eta
          rb(il) = 1.0 - 2.0*rba(il)/wl(ig,il)
          ba(il) = -(code_dt*vn*capg + delp*mu)/eta
          bb(il) = 1.0 - 2.0*(ba(il) + eeb)/wl(ig,il) + ee*(1.+slb1**2)
       
          ra(il) = 0.0
          aa(il) = 0.0

          ! coefficients for entropy heating calculation
          if (heating) then
             dd(il) =vnc*(2.0*(capg + delm*mu)/(eta*wl(ig,il)) + ee)
             hh(il) =vnh*(2.0*(capg + delm*mu)/(eta*wl(ig,il)) + ee)
          end if

          do il = 2, te-1
             
             eea = ee*(1.0 + wslb**2)*delp*mu/eta

             wslb = 1.0 - 0.5*sum(wl(ig,:il))
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
             slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1)))
             
             delp = wslb - slb2
             delm = wslb - slb1
             del  = slb1 - slb2

             capg = 1.0 - wslb**2

             mu = 0.5*(delp + delm)
             eta = (1.0 + wslb*(delp + delm)/capg)*del

!             eeb = -ee*((1.0 + wslb**2)*delp*mu/eta + (1.0 + slb1**2)) - eec(il-1)
             eeb = -ee*(1.0 + wslb**2)*delp*mu/eta - eec(il-1)
             eec(il) = -ee*(1.0 + wslb**2)*delm*mu/eta

             rc(il) = -delm*mu/eta
             cc(il) = -(code_dt*vn*capg + delm*mu)/eta

             rba(il) = -delp*mu/eta
             rb(il) = 1.0 - 2.0*(rba(il) + rc(il-1))/wl(ig,il)
             ba(il) = -(code_dt*vn*capg + delp*mu)/eta
             bb(il) = 1.0 - 2.0*(cc(il-1) + ba(il) + eeb)/wl(ig,il) + ee*(1.+slb1**2)

             ra(il) = 2.0*rba(il-1)/wl(ig,il)
             aa(il) = 2.0*(ba(il-1) + eea)/wl(ig,il)

             ! coefficients for entropy heating calculation
             if (heating) then
                dd(il) =vnc*(2.0*(capg + delm*mu)/(eta*wl(ig,il)) + ee)
                hh(il) =vnh*(2.0*(capg + delm*mu)/(eta*wl(ig,il)) + ee)
             end if

          end do

          ! boundary at xi = 0
          il = te

          eea = ee*(1.0 + wslb**2)*delp*mu/eta
 
          if (te == ng2) then
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
             slb2 = -slb1
             wll = wl(ig,il)
          else
             slb1 = 0.0
             slb2 = -sqrt(abs(1.0 - bmag(ig)*al(il-1)))
             wll = 2.0*wl(ig,il)
          end if
          wslb = 1.0 - 0.5*(sum(wl(ig,:il-1))+wll)

!          ee = 0.5*e(ie,is)*(1+slb1**2) / (bmag(ig)*spec(is)%zstm)**2 &
!               * kperp2(ig,it,ik)*cfac
             
          delp = wslb - slb2
          delm = wslb - slb1
          del  = slb1 - slb2
             
          capg = 1.0 - wslb**2
             
          mu = 0.5*(delp + delm)
          eta = (1.0 + wslb*(delp + delm)/capg)*del
          
!          eeb = -ee*((1.0 + wslb**2)*delp*mu/eta - (1.0 + slb1**2)) + eec(il-1)
          eeb = -ee*(1.0 + wslb**2)*delp*mu/eta - eec(il-1)
          eec(il) = -ee*(1.0 + wslb**2)*delm*mu/eta

          rc(il) = -delm*mu/eta
          cc(il) = -(code_dt*vn*capg + delm*mu)/eta
          
          rba(il) = -delp*mu/eta
          rb(il) = 1.0 - 2.0*(rba(il) + rc(il-1))/wll
          ba(il) = -(code_dt*vn*capg + delp*mu)/eta
          bb(il) = 1.0 - 2.0*(cc(il-1) + ba(il) + eeb)/wll + ee*(1.+slb1**2)
       
          ra(il) = 2.0*rba(il-1)/wll
          aa(il) = 2.0*(ba(il-1) + eea)/wll

          ! coefficients for entropy heating calculation
          if (heating) then
             dd(il) =vnc*(2.0*(capg + delm*mu)/(eta*wll) + ee)
             hh(il) =vnh*(2.0*(capg + delm*mu)/(eta*wll) + ee)
          end if

          cc(:te-1) = 2.0*(cc(:te-1) + eec(:te-1))/wl(ig,:te-1)
          rc(:te-1) = 2.0*rc(:te-1)/wl(ig,:te-1)

          cc(te) = 2.0*(cc(te) + eec(te))/wll
          rc(te) = 2.0*rc(te)/wll

       else
 
! old finite difference scheme

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

          ra = 0.0 ; rb = 1.0 ; rc = 0.0

       end if

! start to fill in the arrays for the tridiagonal
       a1(:te) = aa(:te)
       b1(:te) = bb(:te)
       c1(:te,ilz) = cc(:te)
       
       ra1(:te,ilz) = ra(:te)
       rb1(:te,ilz) = rb(:te)
       rc1(:te,ilz) = rc(:te)
       
       if (heating) then
          d1(:te,ilz) = dd(:te)
          h1(:te,ilz) = hh(:te)
       end if

! assuming symmetry in xi, fill in the rest of the arrays.
       a1(te+1:te2) = cc(teh:1:-1)
       b1(te+1:te2) = bb(teh:1:-1)
       c1(te+1:te2,ilz) =aa(teh:1:-1)
       
       ra1(te+1:te2,ilz) = rc(teh:1:-1)
       rb1(te+1:te2,ilz) = rb(teh:1:-1)
       rc1(te+1:te2,ilz) = ra(teh:1:-1)
       
!       do il = 1, te2
!          write (*,*) ig, il, sqrt(max(0.0,1.0-al(il)*bmag(ig))), 0.5*wl(ig,il), 1.0-0.5*sum(wl(ig,:il))
!       end do
!       write (*,*)

!       do il = 1, te2
!          write (*,*) ig, il, a1(il), b1(il), c1(il,ilz), ra1(il,ilz), rb1(il,ilz), rc1(il,ilz)
!       end do
!       write (*,*)

       if (heating) then
          d1(te+1:te2,ilz) = dd(teh:1:-1)
          h1(te+1:te2,ilz) = hh(teh:1:-1)
          d1(te2+1:,ilz) = 0.0
          h1(te2+1:,ilz) = 0.0
       end if

       betaa(1,ilz) = 1.0/b1(1)
       do il = 1, te2-1
          ql(il+1,ilz) = a1(il+1)*betaa(il,ilz)
          betaa(il+1,ilz) = 1.0/(b1(il+1)-ql(il+1,ilz)*c1(il,ilz))
       end do

       ql(1,ilz) = 0.0
       ql(te2+1:,ilz) = 0.0
       c1(te2+1:,ilz) = 0.0
       betaa(te2+1:,ilz) = 0.0
       
    end do

    if (heating) deallocate (dd, hh)

  end subroutine init_lorentz

  subroutine init_lorentz_error

    use mp, only: proc0
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
    use gs2_layouts, only: gidx2lzidx 
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension(3) :: from_low, from_high
    integer, dimension(2) :: to_high
    integer :: to_low
    integer :: ig, isign, iglo, il, ilz
    integer :: n, ip
    logical :: done = .false.

    if (done) return

    call init_lorentz_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
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
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
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

    done = .true.

  end subroutine init_lorentz_redistribute

  subroutine solfp1 (g, gold, g1, phi, bpar, phinew, bparnew, diagnostics)
    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use run_parameters, only: tunits
    use run_parameters, only: fphi, fbpar
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    use le_grids, only: e, integrate_moment
    use species, only: nspec
    use dist_fn_arrays, only: c_rate
    use mp, only: proc0
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: g, gold, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar, phinew, bparnew
    complex, dimension (:,:,:), allocatable :: gc1, gc2, gc3
    integer, optional, intent (in) :: diagnostics

    integer :: ik, ie, is, iglo

    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot, upar, uperp, ttot
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky,nspec) :: ntot2, upar2, uperp2, ttot2

    logical :: first = .true.

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
    
    if (adjust) call g_adjust (g, phi, bpar, fphi, fbpar)
    if (adjust) call g_adjust (gold, phi, bpar, fphi, fbpar)

    call conserve_test (g, phinew, ntot, upar, uperp, ttot)

    select case (collision_model_switch)
    case (collision_model_full)

       if (heating .and. present(diagnostics)) then
          gc3 = g
          call solfp_lorentz (g, gc1, gc2, diagnostics)
       else
          call solfp_lorentz (g, gc1, gc2)
       end if
       if (conserve_moments) call conserve_lorentz (g, g1)

       call solfp_ediffuse (g)
       if (conserve_moments) call conserve_diffuse (g, g1)

    case (collision_model_lorentz,collision_model_lorentz_test)

       if (heating .and. present(diagnostics)) then
          gc3 = g
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

    call conserve_test (g, phinew, ntot2, upar2, uperp2, ttot2)

    if (.not. allocated(ntot_diff)) then
       allocate (ntot_diff(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (upar_diff(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (uperp_diff(-ntgrid:ntgrid,ntheta0,naky,nspec))
       allocate (ttot_diff(-ntgrid:ntgrid,ntheta0,naky,nspec))
       ntot_diff = 0.0 ; upar_diff = 0.0 ; uperp_diff = 0.0 ; ttot_diff = 0.0
    end if

    ntot_diff = ntot_diff + (ntot2 - ntot)
    upar_diff = upar_diff + (upar2 - upar)
    uperp_diff = uperp_diff + (uperp2 - uperp)
    ttot_diff = ttot_diff + (ttot2 - ttot)

    if (adjust) call g_adjust (g, phi, bpar, -fphi, -fbpar)
    if (adjust) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
    
  end subroutine solfp1

  subroutine conserve_lorentz (g, g1)

    use mp, only: proc0
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: e, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y0
    real, dimension (:,:,:), allocatable :: vns

    integer :: ig, isgn, iglo, ik, ie, il, is, it, all = 1
!    real :: vnm1, vnm2

    allocate (vns(naky,negrid,nspec))

!    vnm1 = vnmult(1)
!    vnm2 = vnmult(2)

!    vns = (vnm1 - vnm2)*vnew_D + vnm2*vnew_s
    vns = vnmult(1)*vnew_D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)
! v0 = nu_s V_perp
! old form had vnew_s instead of vns -- MAB (1.23.08)
! no factor of sqrt( kperp2 smz**2 / bmag ) in v0 because it was absorbed into u0 -- MAB
             g1(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * g(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (g1, v0y0, all)    ! v0y0
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y0

    allocate (v1y0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa
! V_parallel == v_parallel J0
! v1 = nu V_parallel f_0
!
! No factor of sqrt(T/m) here on purpose (see derivation) 
!
             g1(ig,isgn,iglo) = vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * g(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (g1, v1y0, all)    ! v1y0

    deallocate (vns)

!    if (proc0) then
!       write (*,*) v1y0
!    end if

! Conserve momentum:

!    write (*,*) sum(z0)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             g(ig,isgn,iglo) = g(ig,isgn,iglo) + z0(ig,isgn,iglo)*v0y0(ig,it,ik,is) &
                  - z1(ig,isgn,iglo) * v1y0(ig,it,ik,is)
          end do
       end do
    end do

    deallocate (v0y0, v1y0)

  end subroutine conserve_lorentz

  subroutine conserve_diffuse (g, g1)

    use mp, only: proc0
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: e, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y0, v2y0
    real, dimension (:,:,:), allocatable :: vns

    integer :: ig, isgn, iglo, ik, ie, il, is, it, all = 1
!    real :: vnm1, vnm2

    allocate (vns(naky,negrid,nspec))

!    vnm1 = vnmult(1)
!    vnm2 = vnmult(2)

!    vns = (vnm1 - vnm2)*vnew_D + vnm2*vnew_s
    vns = vnmult(2)*(vnew_s - vnew_D)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)
! v0 = nu_s V_perp
! no factor of sqrt( kperp2 smz**2 / bmag ) in v0 because it was absorbed into u0 -- MAB
             g1(ig,isgn,iglo) = vns(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * g(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (g1, v0y0, all)    ! v0y0
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y0 and v2y0 (actually v1 and v2 . (y0 - v0y0 * s0)...i.e. v1y1 and v2y1)

    allocate (v1y0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa
! V_parallel == v_parallel J0
! v1 = nu V_parallel f_0
!
             g1(ig,isgn,iglo) = vns(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * (g(ig,isgn,iglo) - bs0(ig,isgn,iglo)*v0y0(ig,it,ik,is))
          end do
       end do
    end do

    call integrate_moment (g1, v1y0, all)    ! v1y1

    allocate (v2y0(-ntgrid:ntgrid, ntheta0, naky, nspec))         

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v2 = nu_E e J0
             g1(ig,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*e(ie,is)*aj0(ig,iglo) &
                  * (g(ig,isgn,iglo) - bs0(ig,isgn,iglo)*v0y0(ig,it,ik,is))
          end do
       end do
    end do

    call integrate_moment (g1, v2y0, all)    ! v2y1

    deallocate (vns)

! Conserve momentum and energy:

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             g(ig,isgn,iglo) = g(ig,isgn,iglo) - bs0(ig,isgn,iglo)*v0y0(ig,it,ik,is) &
                  - bz0(ig,isgn,iglo)*v2y0(ig,it,ik,is) + (bz0(ig,isgn,iglo)*v2w1(ig,it,ik,is) &
                  - bw0(ig,isgn,iglo))*v1y0(ig,it,ik,is)
          end do
       end do
    end do

    deallocate (v0y0, v1y0, v2y0)

  end subroutine conserve_diffuse

!>MAB
  subroutine conserve_test (g, phi, ntot, upar, uperp, ttot)

    use mp, only: proc0
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ie_idx, il_idx, is_idx
    use gs2_layouts, only: ik_idx, it_idx
    use dist_fn_arrays, only: aj0, aj1, vpa
    use le_grids, only: integrate_moment, e, al, anon
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0

    implicit none

    integer :: ig, iglo, ie, il, is, ik, it, all = 1
    complex, dimension (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc) :: g0

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: ntot, upar, uperp, ttot

    g0 = 0.0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          g0(ig,:,iglo) =  anon(ie,is)*phi(ig,it,ik)*spec(is)%zt*spec(is)%dens
!          g0(ig,:,iglo) =  anon(ie,is)*aj0(ig,iglo)*phi(ig,it,ik)*spec(is)%zt*spec(is)%dens
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig=-ntgrid, ntgrid
          g0(ig,:,iglo) = aj0(ig,iglo)*g(ig,:,iglo) - g0(ig,:,iglo)
!          g0(ig,:,iglo) = g(ig,:,iglo) + g0(ig,:,iglo)
       end do
    end do
!    call integrate_moment (g0, ntot)
    call integrate_moment (g, ntot)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid
          g0(ig,:,iglo) = vpa(ig,:,iglo)*aj0(ig,iglo)*g(ig,:,iglo)
!          g0(ig,:,iglo) = vpa(ig,:,iglo)*g(ig,:,iglo)
       end do
    end do
    call integrate_moment (g0, upar, all)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       il = il_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          g0(ig,:,iglo) = e(ie,is)*al(il)*aj1(ig,iglo)*g(ig,:,iglo)
!          g0(ig,:,iglo) = e(ie,is)*al(il)*g0(ig,:,iglo)
       end do
    end do
    call integrate_moment (g0, uperp, all)

   do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo, iglo)
       il = il_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)
       do ig = -ntgrid, ntgrid
          g0(ig,:,iglo) = e(ie,is)*aj0(ig,iglo)*g(ig,:,iglo) &
               - 0.5*anon(ie,is)*phi(ig,it,ik)*spec(is)%zt*spec(is)%dens
!          g0(ig,:,iglo) = e(ie,is)*al(il)*g0(ig,:,iglo)
       end do
    end do
    call integrate_moment (g0, ttot, all)

  end subroutine conserve_test
!<MAB

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

  subroutine solfp_lorentz (g, gc, gh, diagnostics)
    use species, only: spec, electron_species
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: nlambda, jend, lintegrate, ng2, al
    use gs2_layouts, only: g_lo, gint_lo, lz_lo
    use gs2_layouts, only: ig_idx, ik_idx, il_idx, is_idx, it_idx, ie_idx
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gc, gh
    integer, optional, intent (in) :: diagnostics

    complex, dimension (:,:), allocatable :: glz, glzc
    complex, dimension (:), allocatable :: glz0
    complex, dimension (max(2*nlambda,2*ng2+1)) :: delta
    complex :: fac
    integer :: iglo, igint, ilz, ig, ik, il, is, je, it, ie

    logical :: first = .true.
    
    call prof_entering ("solfp_lorentz", "collisions")

    allocate (glz(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (glz0(max(2*nlambda,2*ng2+1)))
    glz = 0.0 ; glz0 = 0.0
    if (heating) then
       allocate (glzc(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       glzc = 0.0
    end if

!    call check_g ('beg', g)

    ! TEMPORARY FOR TESTING -- MAB
!    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!       ie = ie_idx(g_lo,iglo)
!       il = il_idx(g_lo,iglo)
!       do ig = -ntgrid, ntgrid
!          g(ig,1,iglo) = il
!          g(ig,2,iglo) = -il
!       end do
!    end do

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

       je = jend(ig)

       if (je > ng2+1) then
          je = 2*je
       else
          je = 2*ng2+1
       end if

!       je = 2*jend(ig)

!       if (je <= ng2) then
!          je = 2*ng2+1
!       end if

       ! deal with special case of wfb
       if (jend(ig) == ng2+1) then
          ! if wfb, remove vpa = 0 point (which has wgt of zero)
          glz0(:ng2) = glz(:ng2,ilz)
          glz0(ng2+1:je-1) = glz(ng2+2:je,ilz)
       else
          glz0 = glz(:,ilz)
       end if

!       do il = 1, je
!          if (ilz == 1) write (*,*) ie, il, ig, real(glz0(il)), 'glz(il,ilz=1)', sqrt(max(0.0,1.0-al(il)*bmag(ig)))
!       end do

       glz(1,ilz) = rb1(1,ilz)*glz0(1) + rc1(1,ilz)*glz0(2)
       glz(je-1,ilz) = ra1(je-1,ilz)*glz0(je-2) + rb1(je-1,ilz)*glz0(je-1)
       do il = 2, je-2
          glz(il,ilz) = ra1(il,ilz)*glz0(il-1) + rb1(il,ilz)*glz0(il) + rc1(il,ilz)*glz0(il+1)
       end do

       glz(je:,ilz) = 0.0

!       do il = 1, je
!          if (ilz == 1) write (*,*) ie, il, ig, real(glz(il,ilz)), 'glz0(il,ilz=1)', sqrt(max(0.0,1.0-al(il)*bmag(ig)))
!       end do

       ! right and left sweeps for tridiagonal solve:

       delta(1) = glz(1,ilz)
       do il = 1, je-1
          delta(il+1) = glz(il+1,ilz) - ql(il+1,ilz)*delta(il)
       end do
       
       glz(je,ilz) = delta(je)*betaa(je,ilz)
       do il = je-1, 1, -1
          glz(il,ilz) = (delta(il) - c1(il,ilz)*glz(il+1,ilz))*betaa(il,ilz)
       end do

       ! interpolate to obtain glz(vpa = 0) point for wfb
       ! and insert this point into glz
       if (jend(ig) == ng2+1) then
          glz0(ng2+2:je) = glz(ng2+1:je-1,ilz)
          glz(ng2+1,ilz) = 0.5*(glz(ng2,ilz)+glz(ng2+1,ilz))
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

!       do il = 1, je
!          if (ilz == 1) write (*,*) ie, il, ig, real(glz(il,ilz)), 'final glz', sqrt(max(0.0,1.0-al(il)*bmag(ig)))
!       end do

    end do

!    call check_glz ('end', glz)
    call scatter (lorentz_map, glz, g)

!    call check_g ('end', g)

    deallocate (glz, glz0)
    if (heating) deallocate (glzc)

    call prof_leaving ("solfp_lorentz", "collisions")
  end subroutine solfp_lorentz

  subroutine solfp_ediffuse (g)
    use mp, only: proc0
    use species, only: spec, nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use le_grids, only: negrid, integrate_moment!, e
    use gs2_layouts, only: is_idx, e_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use file_utils, only: open_output_file, close_output_file
    use egrid, only: zeroes, x0

    use gs2_layouts, only: ik_idx, it_idx, is_idx, ie_idx, il_idx, ig_idx

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g

    complex, dimension (:,:), allocatable :: ged
    complex, dimension (:), allocatable :: ged0
    complex, dimension (negrid) :: delta
    integer :: ielo, ie, is, iglo, ig, isgn, it, ik, il

    allocate (ged(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc))
    allocate (ged0(negrid+1))
    ged = 0.0 ; ged0 = 0.0

    ! TEMP FOR TESTING -- MAB
!    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!       ie = ie_idx(g_lo,iglo)
!       g(:,:,iglo) = ie
!    end do

    call gather (ediffuse_map, g, ged)

    ! solve for ged row by row
    do ielo = e_lo%llim_proc, e_lo%ulim_proc
       is = is_idx(e_lo,ielo)
       ! TEMP -- MAB
!       ik = ik_idx(e_lo,ielo)
!       ig = ig_idx(e_lo,ielo)
!       il = il_idx(e_lo,ielo)
!       it = it_idx(e_lo,ielo)

       if (spec(is)%vnewk < 2.0*epsilon(0.0)) cycle

!       do ie = 1, negrid
!          write (*,*) ie, ged(ie,ielo), 'pre', ielo
!       end do

       ged0(1) = erb1(1,ielo)*ged(1,ielo) + erc1(1,ielo)*ged(2,ielo)
       do ie = 2, negrid-1
          ged0(ie) = era1(ie,ielo)*ged(ie-1,ielo) + erb1(ie,ielo)*ged(ie,ielo) + erc1(ie,ielo)*ged(ie+1,ielo)
       end do
       ged0(negrid) = era1(negrid,ielo)*ged(negrid-1,ielo) + erb1(negrid,ielo)*ged(negrid,ielo)

!       do ie = 1, negrid
!          write (*,*) ie, ged0(ie), 'intermediate', ielo
!       end do

       ged(:,ielo) = ged0

       delta(1) = ged(1,ielo)
! should be ie = 1, negrid? and pad ged, eql, etc. with extra element?
       do ie = 1, negrid-1
          delta(ie+1) = ged(ie+1,ielo) - eql(ie+1,ielo)*delta(ie)
       end do
       
! also following should change to mimic lorentz?
       ged(negrid+1,ielo) = 0.0
       do ie = negrid, 1, -1
          ged(ie,ielo) = (delta(ie) - ec1(ie,ielo)*ged(ie+1,ielo))*ebetaa(ie,ielo)
       end do

!       do ie = 1, negrid
!          write (*,*) ie, ged(ie,ielo), 'after', ielo
!       end do

    end do

    call scatter (ediffuse_map, ged, g)

    deallocate (ged, ged0)

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

  subroutine reset_init
!
! forces recalculation of coefficients in collision operator
! when timestep changes.
!    
    initialized = .false.  

  end subroutine reset_init

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

