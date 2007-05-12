module collisions

  use redistribute, only: redist_type

  implicit none

  public :: init_collisions
  public :: solfp1
  public :: reset_init
  public :: dtot, fdf, fdb, lorentz_map, vnmult, vnfac
  public :: ewindow, etol, ncheck, vnslow, vary_vnew
  public :: init_lorentz, init_lz_mom_conserve, init_escatter
  public :: init_lorentz_error, collision_model_switch

  private

  ! knobs
  real :: vncoef, absom
  integer :: ivnew
  logical :: conserve_number, conserve_momentum, const_v
  integer :: collision_model_switch
  logical :: use_shmem, adjust
  logical :: heating
  logical :: hyper_colls
  logical :: scatter_energy

  integer, parameter :: collision_model_lorentz = 1      ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_krook = 2
  integer, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_krook_test = 4
  integer, parameter :: collision_model_lorentz_test = 5 ! if this changes, check gs2_diagnostics

  real, dimension (2), save :: vnmult
  integer :: ncheck
  logical :: vary_vnew
  real :: vnfac, etol, ewindow, vnslow

  real, dimension (:,:,:), allocatable :: dtot
  ! (-ntgrid:ntgrid,nlambda,max(ng2,nlambda-ng2)) lagrange coefficients for derivative error estimate

  real, dimension (:,:), allocatable :: fdf, fdb
  ! (-ntgrid,ntgrid,nlambda) finite difference coefficients for derivative error estimate

  real, dimension (:,:,:), allocatable :: vnew, vnew_ss
  ! (naky,negrid,nspec) replicated

  ! only for hyper-diffusive collisions
  real, dimension (:,:,:,:), allocatable :: vnewh
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec) replicated

  ! only for krook
  real, dimension (:), allocatable :: vnewfe
  ! (-*- g_layout -*-)

  ! only for krook
  real, dimension (:,:), allocatable :: aintnorm
  ! (-ntgrid:ntgrid, -*- geint_layout -*-)

  ! only for "new" momentum conservation (8.06)
  complex, dimension(:,:,:), allocatable :: z0, z1

  ! only for original parallel mom conservation (not used nowadays)
  real, dimension (:,:,:), allocatable :: sq
  ! (-ntgrid:ntgrid,nlambda,2) replicated

  ! only for energy scattering
  real, dimension (:), allocatable :: ec1, ebetaa, eql
  complex, dimension (:,:), allocatable :: gesc

  ! only for lorentz
  real :: cfac
  real, dimension (:,:), allocatable :: c1, betaa, ql, d1, h1
  complex, dimension (:,:), allocatable :: glz, glzc
  ! (max(2*nlambda,2*ng2+1), -*- lz_layout -*-)
  ! (-ntgrid:ntgrid,naky,max(2*nlambda,2*ng2+1),negrid,nspec)

  ! momentum conservation
  complex, dimension (:,:), allocatable :: g3int
  ! (-ntgrid:ntgrid, -*- gint_layout -*-)

  type (redist_type), save :: lorentz_map
  type (redist_type), save :: escatter_map

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
    type (text_option), dimension (7), parameter :: modelopts = &
         (/ text_option('default', collision_model_lorentz), &
            text_option('lorentz', collision_model_lorentz), &
            text_option('krook', collision_model_krook), &
            text_option('krook-test', collision_model_krook_test), &
            text_option('lorentz-test', collision_model_lorentz_test), &
            text_option('none', collision_model_none), &
            text_option('collisionless', collision_model_none) /)
    character(20) :: collision_model
    namelist /collisions_knobs/ collision_model, vncoef, absom, ivnew, &
         conserve_number, conserve_momentum, use_shmem, heating, &
         adjust, const_v, cfac, hypermult, scatter_energy, vnfac, &
         etol, ewindow, ncheck, vnslow, vary_vnew
    integer :: ierr, in_file
    logical :: exist

    if (proc0) then
       hypermult = .false.
       cfac = 1.   ! DEFAULT CHANGED TO INCLUDE CLASSICAL DIFFUSION: APRIL 18, 2006
       vnfac = 1.
       vnslow = 0.9
       etol = 2.e-2
       ewindow = 1.e-2
       ncheck = 100
       adjust = .true.
       collision_model = 'default'
       vncoef = 0.6
       absom = 0.5
       ivnew = 0
       conserve_number = .true.
       conserve_momentum = .true.  ! DEFAULT CHANGED TO REFLECT IMPROVED MOMENTUM CONSERVATION, 8/06
       scatter_energy = .false.
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
    call broadcast (ncheck)
    call broadcast (vncoef)
    call broadcast (absom)
    call broadcast (ivnew)
    call broadcast (conserve_number)
    call broadcast (conserve_momentum)
    call broadcast (scatter_energy)
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
    case (collision_model_lorentz,collision_model_lorentz_test)
       call init_lorentz
       if (scatter_energy) call init_escatter
       if (conserve_momentum) call init_lz_mom_conserve
    case (collision_model_krook,collision_model_krook_test)
       call init_krook (hee)
    end select

  end subroutine init_arrays

  subroutine init_lz_mom_conserve

!
! Precompute two quantities needed for momentum conservation:
! z0, z1
    
    use mp, only: proc0
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use le_grids, only: e, al, integrate_moment
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, kperp2, vpa

    logical, save :: first = .true.
    complex, dimension (1,1,1) :: dum1 = 0., dum2 = 0.
    complex, dimension (:,:,:), allocatable :: gtmp
    complex, dimension (:,:,:,:), allocatable :: v0z0, v0z1, v1z0, v1z1, Enuinv
    real :: vnm
    integer :: ie, il, ik, is, ig, isgn, iglo, all, it

    vnm = vnmult(1)

! TO DO: 
! tunits not included anywhere yet

    if (first) then
       allocate (z0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (z1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       first = .false.
    end if

! First, get Enu and then 1/Enu == Enuinv

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (Enuinv(-ntgrid:ntgrid, ntheta0, naky, nspec))       

!
! Enu == int (E nu f_0);  Enu = Enu(z, kx, ky, s)
! Enuinv = 1/Enu
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             gtmp(ig,isgn,iglo) = e(ie,is)*vnm*vnew_ss(ik,ie,is)
          end do
       end do
    end do

    all = 1
    call integrate_moment (gtmp, Enuinv, all)  ! not 1/Enu yet

    where (cabs(Enuinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
                                   ! Enuinv=0 iff vnew=0 so ok to keep Enuinv=0.
       Enuinv = 1./Enuinv  ! now it is 1/Enu
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
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)
! u0 = -3 nu V_perp dt a f_0 / Enu
! where a = kperp2 * (T m / q**2)  ! missing factor of 1/B(theta)**2 ???
             z0(ig,isgn,iglo) = - 3.*vnm*vnew_ss(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
                  * code_dt * spec(is)%smz**2 * kperp2(ig,it,ik) * Enuinv(ig,it,ik,is)
          end do
       end do
    end do

    call solfp_lorentz (z0, dum1, dum2)   ! z0 is redefined below

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get z1 (first form)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
! v_parallel == vpa 
! V_parallel == v_parallel J0
! u1 = -3 nu V_parallel dt f_0 / Enu
!
! No factor of sqrt(T/m) here on purpose (see derivation) 
!
             z1(ig,isgn,iglo) = - 3.*vnm*vnew_ss(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * code_dt * Enuinv(ig,it,ik,is)
          end do
       end do
    end do

    deallocate (Enuinv)  ! Done with this variable

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
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)
! v0 = nu V_perp
             gtmp(ig,isgn,iglo) = vnm*vnew_ss(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
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
             gtmp(ig,isgn,iglo) = vnm*vnew_ss(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
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
! V_perp == e(ie,is)*al(il)*aj1(ig,iglo)
! v0 = nu V_perp
             gtmp(ig,isgn,iglo) = vnm*vnew_ss(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
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
             gtmp(ig,isgn,iglo) = vnm*vnew_ss(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * z1(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (gtmp, v1z1, all)    ! redefined below

    deallocate (gtmp)
    
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
    
  end subroutine init_lz_mom_conserve

  subroutine init_vnew (hee)
    use species, only: nspec, spec, electron_species, has_electron_species
    use le_grids, only: negrid, e
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use run_parameters, only: zeff, tunits
    use dist_fn_arrays, only: kperp2
    use constants
    real, dimension (:,:), intent (out) :: hee
    real,dimension (negrid,nspec)::heevth
    integer :: ik, ie, is, it, ig
    real :: v, k4max

    do is = 1, nspec
       do ie = 1, negrid
          v = sqrt(e(ie,is))
          hee(ie,is) = 1.0/sqrt(pi)/v*exp(-e(ie,is)) &
               + (1.0 - 0.5/e(ie,is)) &
               *(1.0 - 1.0/(1.0          + v &
               *(0.0705230784 + v &
               *(0.0422820123 + v &
               *(0.0092705272 + v &
               *(0.0001520143 + v &
               *(0.0002765672 + v &
               *(0.0000430638)))))))**16)
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
       end do
    end do                                                                  

!    do is = 1, nspec
!       if (spec(is) % nustar < 0.)
!
!    end do

    if(.not.allocated(vnew)) then
       allocate (vnew(naky,negrid,nspec))
       allocate (vnew_ss(naky,negrid,nspec))
    end if
    if(.not.allocated(vnewh)) allocate (vnewh(-ntgrid:ntgrid,ntheta0,naky,nspec))

    do is = 1, nspec
       if (spec(is)%type == electron_species) then
          do ie = 1, negrid
             do ik = 1, naky
                if (const_v) then
                   vnew(ik,ie,is) = spec(is)%vnewk &
                        *(zeff+heevth(ie,is))*0.5*tunits(ik)                   
                   vnew_ss(ik,ie,is) = spec(is)%vnewk &
                        *heevth(ie,is)*0.5*tunits(ik)                   
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *(zeff + hee(ie,is))*0.5*tunits(ik)
                   vnew_ss(ik,ie,is) = spec(is)%vnewk/e(ie,is)**1.5 &
                        *hee(ie,is)*0.5*tunits(ik)
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
                vnew_ss(ik,ie,is) = vnew(ik,ie,is)
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

  subroutine init_escatter (vnmult_target)
    use constants, only: pi
    use mp, only: proc0
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, ng2, e, ecut
    use egrid, only: zeroes, x0, energy
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: init_escatter_layouts
    use gs2_layouts, only: e_lo
    use file_utils, only: open_output_file, close_output_file

    implicit none
    
    real, intent (in), optional :: vnmult_target

    integer :: ie, is
    real, dimension (:), allocatable :: aa, bb, cc, xe, cder
    real :: vn, xe0, xe1, xe2, xer, xel, er, el, fac

    integer, save :: tmp_unit
    logical :: first_time = .true.

    if (first_time) then
       vnmult(2) = max(1.0, vnmult(2))
       first_time = .false.
    end if

    call init_escatter_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)
    call init_escatter_redistribute

    if (.not.allocated(gesc)) then
       allocate (gesc(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc))
       gesc = 0.0
    end if

    allocate (aa(negrid), bb(negrid), cc(negrid), cder(negrid))
    allocate (xe(negrid))
    xe(:negrid-1) = zeroes
    xe(negrid) = x0

    if (.not.allocated(ec1)) then
       allocate (ec1   (negrid))
       allocate (ebetaa(negrid))
       allocate (eql   (negrid))
    endif

    if (present(vnmult_target)) then
       vnmult(2) = max (vnmult_target, 1.0)
    end if

    ec1 = 0.0 ; ebetaa = 0.0 ; eql = 0.0
    is = 1

! should include multiplication by tunits(ik) in general
    vn = vnmult(2)*spec(is)%vnewk

    do ie = 2, negrid-1

       xe0 = xe(ie-1)
       xe1 = xe(ie)
       xe2 = xe(ie+1)

       xel = (xe1 + xe0)*0.5
       xer = (xe2 + xe1)*0.5

       el = energy(xel,ecut)
       er = energy(xer,ecut)

       fac = -8.0*vn*code_dt/pi/(xer - xel)

       aa(ie) = fac*xel*sqrt(el)*exp(-2.0*el)/(xe1 - xe0)
       cc(ie) = fac*xer*sqrt(er)*exp(-2.0*er)/(xe2 - xe1)
       bb(ie) = 1.0 - (aa(ie) + cc(ie))

    end do

! boundary at xe = 0
    xe0 = 0.0
    xe1 = xe(1)
    xe2 = xe(2)

    xer = (xe2 + xe1)*0.5

    er = energy(xer,ecut)

    fac = -8.0*vn*code_dt/pi/xer

    aa(1) = 0.0
    cc(1) = fac*xer*sqrt(er)*exp(-2.0*er)/(xe2 - xe1)
    bb(1) = 1.0 - cc(1)

! boundary at xe = 1

    xe0 = xe(negrid-1)
    xe1 = xe(negrid)
    xe2 = 1.0

    xel = (xe1 + xe0)*0.5

    el = energy(xel,ecut)

    fac = 8.0*vn*code_dt/pi/(xe2 - xel)

!    aa(negrid) = fac*xel*sqrt(el)*exp(-2.0*el)/(xe1 - xe0)
! derivative set to zero at xe(negrid)=x0 because we assume
! d h(xe) / d xe is zero between x0 and 1 in our integration scheme 
    aa(negrid) = 0.0
    cc(negrid) = 0.0
    bb(negrid) = 1.0 - aa(negrid)

! fill in the arrays for the tridiagonal
    ec1 = cc
    ebetaa(1) = 1.0/bb(1)
    do ie = 1, negrid-1
       eql(ie+1) = aa(ie+1)*ebetaa(ie)
       ebetaa(ie+1) = 1.0/(bb(ie+1)-eql(ie+1)*ec1(ie))
    end do

    deallocate(aa, bb, cc, xe, cder)

  end subroutine init_escatter

  subroutine init_escatter_redistribute
    use mp, only: nproc
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, ng2
    use gs2_layouts, only: init_escatter_layouts
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

    call init_escatter_layouts &
         (ntgrid, naky, ntheta0, nlambda, nspec)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
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

    call init_redist (escatter_map, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

    done = .true.

  end subroutine init_escatter_redistribute

  subroutine init_lorentz (vnmult_target)
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, al, jend, ng2,e
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: kperp2
    use gs2_layouts, only: init_lorentz_layouts
    use gs2_layouts, only: lz_lo
    use gs2_layouts, only: ig_idx, ik_idx, ie_idx, is_idx, it_idx
    implicit none

    real, intent (in), optional :: vnmult_target

    integer :: ig, il, ilz, it, ik, ie, is, je
    real, dimension (nlambda+1) :: aa, bb, cc, dd, hh
    real, dimension (max(2*nlambda,2*ng2+1)) :: a1, b1
    real :: slb0, slb1, slb2, slbl, slbr, vn, ee, vnh, vnc
    logical :: first_time = .true.

    if (first_time) then
       vnmult(1) = max(1.0, vnmult(1))
       first_time = .false.
    end if

    call init_lorentz_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    call init_lorentz_redistribute

    if (.not.allocated(glz)) then
       allocate (glz(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       glz = 0.0
       if (heating) then
          allocate (glzc(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
          glzc = 0.0
       end if
    end if

    if (.not.allocated(c1)) then
       allocate (c1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (betaa(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (ql   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (d1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       allocate (h1   (max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
    endif

    c1 = 0.0 ; betaa = 0.0 ; ql = 0.0 ; d1 = 0.0 ; h1 = 0.0

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
! no trapped particles if je == 0 
       if (je == 0) then  
          do il = 2, ng2-1
             ! slb = xi = v_par/v
             slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))   ! xi_{j-1}
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))     ! xi_j
             slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1)))   ! xi_{j+1}

!             write (*,*) il,' xi = ',slb1
             slbl = (slb1 + slb0)/2.0  ! xi(j-1/2)
             slbr = (slb1 + slb2)/2.0  ! xi(j+1/2)

             ee = 0.25*e(ie,is)*(1+slb1**2) &
                  / (bmag(ig)*spec(is)%zstm)**2 &
                  * kperp2(ig,it,ik)*cfac

             ! coefficients for tridiagonal matrix:
             cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
             bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt

             ! coefficients for entropy heating calculation
             dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
             hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          end do

! boundary at xi = 1
          slb0 = 1.0
          slb1 = sqrt(abs(1.0-bmag(ig)*al(1)))
          slb2 = sqrt(abs(1.0-bmag(ig)*al(2)))

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

          ee = 0.25*e(ie,is)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          cc(1) = -vn*code_dt*(-1.0 - slbr)/(slb2-slb1)
          aa(1) = 0.0
          bb(1) = 1.0 - (aa(1) + cc(1)) + ee*vn*code_dt

          dd(1) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(1) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)

! boundary at xi = 0
          il = ng2
          slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))
          slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
          slb2 = -slb1

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

          ee = 0.25*e(ie,is)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac

          cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
          aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt

          dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)

! start to fill in the arrays for the tridiagonal
          a1(:ng2) = aa(:ng2)
          b1(:ng2) = bb(:ng2)
          c1(:ng2,ilz) = cc(:ng2)

          d1(:ng2,ilz) = dd(:ng2)
          h1(:ng2,ilz) = hh(:ng2)

! assuming symmetry in xi, fill in the rest of the arrays.
          a1(ng2+1:2*ng2) = cc(ng2:1:-1)
          b1(ng2+1:2*ng2) = bb(ng2:1:-1)
          c1(ng2+1:2*ng2,ilz) =aa(ng2:1:-1)

          d1(ng2+1:2*ng2,ilz) = dd(ng2:1:-1)
          h1(ng2+1:2*ng2,ilz) = hh(ng2:1:-1)

          betaa(1,ilz) = 1.0/b1(1)
          do il = 1, 2*ng2-1
             ql(il+1,ilz) = a1(il+1)*betaa(il,ilz)
             betaa(il+1,ilz) = 1.0/(b1(il+1)-ql(il+1,ilz)*c1(il,ilz))
          end do

          ql(1,ilz) = 0.0
          ql(2*ng2+1:,ilz) = 0.0
          c1(2*ng2+1:,ilz) = 0.0
          betaa(2*ng2+1:,ilz) = 0.0

          d1(2*ng2+1:,ilz) = 0.0
          h1(2*ng2+1:,ilz) = 0.0
       else
          do il = 2, je-1
             slb0 = sqrt(abs(1.0 - bmag(ig)*al(il-1)))
             slb1 = sqrt(abs(1.0 - bmag(ig)*al(il)))
             slb2 = sqrt(abs(1.0 - bmag(ig)*al(il+1)))

             slbl = (slb1 + slb0)/2.0
             slbr = (slb1 + slb2)/2.0

             ee = 0.25*e(ie,is)*(1+slb1**2) &
                  / (bmag(ig)*spec(is)%zstm)**2 &
                  * kperp2(ig,it,ik)*cfac

             cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
             bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt

             dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
             hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          end do

          slb0 = 1.0
          slb1 = sqrt(abs(1.0-bmag(ig)*al(1)))
          slb2 = sqrt(abs(1.0-bmag(ig)*al(2)))

          slbr = (slb1 + slb2)/2.0

          ee = 0.25*e(ie,is)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac

          cc(1) = -vn*code_dt*(-1.0 - slbr)/(slb2-slb1)
          aa(1) = 0.0
          bb(1) = 1.0 - (aa(1) + cc(1)) + ee*vn*code_dt

          dd(1) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(1) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)

          il = je
          slb0 = sqrt(abs(1.0-bmag(ig)*al(il-1)))
          slb1 = 0.
          slb2 = -slb0                                                        

          ee = 0.25*e(ie,is)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

! Are cc(il) and aa(il) missing a factor of 2 here? I think it should be:
! cc(il) = -0.5*vn*code_dt*(1.0-slbl*slbl)/slb0/slb0...MAB
! was originally cc(il)-0.5*vn*code_dt*(1.0-slbl*slbl)/slbl/slb0
!          cc(il) = -0.5*vn*code_dt*(1.0-slbl*slbl)/slb0/slb0  ! NEW LINE
          cc(il) = -0.5*vn*code_dt*(1.0-slbl*slbl)/slbl/slb0   ! OLD LINE
          aa(il) = cc(il)
          bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt

          dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)

          a1(:je) = aa(:je)
          b1(:je) = bb(:je)
          c1(:je,ilz) = cc(:je)

          d1(:je,ilz) = dd(:je)
          h1(:je,ilz) = hh(:je)

          a1(je+1:2*je-1) = cc(je-1:1:-1)
          b1(je+1:2*je-1) = bb(je-1:1:-1)
          c1(je+1:2*je-1,ilz) = aa(je-1:1:-1)

          d1(je+1:2*je-1,ilz) = dd(je-1:1:-1)
          h1(je+1:2*je-1,ilz) = hh(je-1:1:-1)

          betaa(1,ilz) = 1.0/b1(1)
          do il = 1, 2*je-2
             ql(il+1,ilz) = a1(il+1)*betaa(il,ilz)
             betaa(il+1,ilz) = 1.0/(b1(il+1)-ql(il+1,ilz)*c1(il,ilz))
          end do

          ql(1,ilz) = 0.0
          c1(2*je:,ilz) = 0.0
          betaa(2*je:,ilz) = 0.0
          d1(2*je:,ilz) = 0.0
          h1(2*je:,ilz) = 0.0
       end if
    end do

  end subroutine init_lorentz

  subroutine init_lorentz_error

    use mp, only: proc0
    use le_grids, only: jend, al, ng2, nlambda
    use theta_grid, only: ntgrid, bmag
    implicit none
    
    integer :: je, te, ig, il, ip, ij, im
    real :: slb0, slb1, slb2, slbr, slbl
    real, dimension (:), allocatable :: slb
    real, dimension (:,:), allocatable :: dprod
    real, dimension (:,:,:), allocatable :: dlcoef, d2lcoef

    allocate(slb(2*nlambda))
    allocate(dprod(2*nlambda,max(ng2,2*(nlambda-ng2))))

    allocate (dlcoef(-ntgrid:ntgrid,2*nlambda-ng2,max(ng2,2*(nlambda-ng2))))
    allocate (d2lcoef(-ntgrid:ntgrid,2*nlambda-ng2,max(ng2,2*(nlambda-ng2))))
    allocate (dtot(-ntgrid:ntgrid,2*nlambda-ng2,max(ng2,2*(nlambda-ng2))))
    allocate (fdf(-ntgrid:ntgrid,nlambda), fdb(-ntgrid:ntgrid,nlambda))

    dlcoef = 1.0; d2lcoef = 0.0; dtot = 0.0
    fdf = 0.0; fdb = 0.0; slb = 0.0

    do ig=-ntgrid,ntgrid
       je = jend(ig)
       
       if (je == 0) then            ! no trapped particles

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
       end if

! compute coefficients (dlcoef) multipyling first derivative of h
       do ip=1,ng2
          do il=1,ng2
             if (il == ip) then
                dlcoef(ig,il,ip) = 0.0
                do ij=1,ng2
                   if (ij /= ip) dlcoef(ig,il,ip) = dlcoef(ig,il,ip) + 1/(slb(il)-slb(ij))
                end do
             else
                do ij=1,ng2
                   if (ij /= ip .and. ij /= il) then
                      dlcoef(ig,il,ip) = dlcoef(ig,il,ip)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                   end if
                end do
                dlcoef(ig,il,ip) = dlcoef(ig,il,ip)/(slb(ip)-slb(il))
             end if
             dlcoef(ig,il,ip) = -2.0*slb(il)*dlcoef(ig,il,ip)
          end do
       end do

       dprod = 2.0

! compute coefficients (d2lcoef) multiplying second derivative of h
       do ip=1,ng2
          do il=1,ng2
             if (il == ip) then
                do ij=1,ng2
                   if (ij /= ip) then
                      do im=1,ng2
                         if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                      end do
                   end if
                end do
             else
                do ij=1,ng2
                   if (ij /= il .and. ij /= ip) then
                      dprod(il,ip) = dprod(il,ip)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                   end if
                end do

                do ij=1,ng2
                   if (ij /= ip .and. ij /= il) then
                      d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./(slb(il)-slb(ij))
                   end if
                end do
                d2lcoef(ig,il,ip) = dprod(il,ip)*d2lcoef(ig,il,ip)/(slb(ip)-slb(il))
             end if
             d2lcoef(ig,il,ip) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip)
          end do
       end do

       if (je /= 0) then      ! have to handle trapped particles
          te = 2*je-ng2-1
          slb(je+1:te) = -slb(je-1:ng2+1:-1)          

          do ip=1,te-ng2
             do il=ng2+1,te
                if (il == ip+ng2) then
                   dlcoef(ig,il,ip) = 0.0
                   do ij=ng2+1,te
                      if (ij /= ip+ng2) dlcoef(ig,il,ip) = dlcoef(ig,il,ip) + 1/(slb(il)-slb(ij))
                   end do
                else
                   do ij=ng2+1,te
                      if (ij /= ip+ng2 .and. ij /= il) then
                         dlcoef(ig,il,ip) = dlcoef(ig,il,ip)*(slb(il)-slb(ij))/(slb(ip+ng2)-slb(ij))
                      end if
                   end do
                   dlcoef(ig,il,ip) = dlcoef(ig,il,ip)/(slb(ip+ng2)-slb(il))
                end if
                dlcoef(ig,il,ip) = -2.0*slb(il)*dlcoef(ig,il,ip)
             end do
          end do
          
          do ip=1,te-ng2
             do il=ng2+1,te
                if (il == ip+ng2) then
                   do ij=ng2+1,te
                      if (ij /= ip+ng2) then
                         do im=ng2+1,te
                            if (im /= ip+ng2 .and. im /= ij) d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                         end do
                      end if
                   end do
                else
                   do ij=ng2+1,te
                      if (ij /= il .and. ij /= ip+ng2) then
                         dprod(il,ip) = dprod(il,ip)*(slb(il)-slb(ij))/(slb(ip+ng2)-slb(ij))
                      end if
                   end do
                   
                   do ij=ng2+1,te
                      if (ij /= ip+ng2 .and. ij /= il) then
                         d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./(slb(il)-slb(ij))
                      end if
                   end do
                   d2lcoef(ig,il,ip) = dprod(il,ip)*d2lcoef(ig,il,ip)/(slb(ip+ng2)-slb(il))
                end if
                d2lcoef(ig,il,ip) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip)
             end do
          end do
       end if
    end do

    dtot = dlcoef + d2lcoef

    deallocate(slb,dprod,dlcoef,d2lcoef)
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
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: g, gold, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar, phinew, bparnew
    complex, dimension (:,:,:), allocatable :: gc1, gc2, gc3
    integer, optional, intent (in) :: diagnostics

    integer :: ik, ie, is, iglo

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

    select case (collision_model_switch)
    case (collision_model_lorentz,collision_model_lorentz_test)
       if (present(diagnostics)) then
          gc3 = g
          call solfp_lorentz (g, gc1, gc2, diagnostics)
       else
          call solfp_lorentz (g, gc1, gc2)
       end if

       if (conserve_momentum) call conserve_mom (g, g1)

       if (scatter_energy) call solfp_escatter (g)

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

    if (adjust) call g_adjust (g, phi, bpar, -fphi, -fbpar)
    if (adjust) call g_adjust (gold, phi, bpar, -fphi, -fbpar)
    
  end subroutine solfp1

  subroutine conserve_mom (g, g1)

    use mp, only: proc0
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: e, al, integrate_moment
    use dist_fn_arrays, only: aj0, aj1, vpa
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y0

    integer :: ig, isgn, iglo, ik, ie, il, is, it, all = 1

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
! v0 = nu V_perp
             g1(ig,isgn,iglo) = vnew_ss(ik,ie,is)*e(ie,is)*al(il)*aj1(ig,iglo) &
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
             g1(ig,isgn,iglo) = vnew_ss(ik,ie,is)*vpa(ig,isgn,iglo)*aj0(ig,iglo) &
                  * g(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (g1, v1y0, all)    ! v1y0

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
             g(ig,isgn,iglo) = g(ig,isgn,iglo) + 2.0*z0(ig,isgn,iglo)*v0y0(ig,it,ik,is) &
                  - 2.0*z1(ig,isgn,iglo) * v1y0(ig,it,ik,is)
          end do
       end do
    end do

    deallocate (v0y0, v1y0)

  end subroutine conserve_mom

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
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, jend, lintegrate, ng2, al
    use gs2_layouts, only: g_lo, gint_lo, lz_lo
    use gs2_layouts, only: ig_idx, ik_idx, il_idx, is_idx, it_idx, ie_idx
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use file_utils, only: open_output_file, close_output_file
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gc, gh
    integer, optional, intent (in) :: diagnostics

    complex, dimension (max(2*nlambda,2*ng2+1)) :: delta
    complex :: fac
    integer :: iglo, igint, ilz, ig, ik, il, is, je, it

    integer, save :: tmp_unit
    logical :: first = .true.
    
    call prof_entering ("solfp_lorentz", "collisions")

!    call check_g ('beg', g)

    call gather (lorentz_map, g, glz)

!    if (first) then
!       call open_output_file (tmp_unit,".tmp")       
!       do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
!          ig = ig_idx(lz_lo,ilz)
!          ie = ie_idx(lz_lo,ilz)
!          if (ig == 0 .and. ie == 1) then
!             je = 2*jend(ig)
!             if (je==0) je=2*ng2+1
!             do il=1,je
!                write(tmp_unit,*) jend(ig), ng2, je, il, al(il), glz(il, ilz)
!             end do
!          end if
!       end do
!       do iglo = g_lo%llim_proc, g_lo%ulim_proc
!          ie = ie_idx(g_lo,iglo)
!          il = il_idx(g_lo,iglo)
!          if (ie == 1) then
!             write(tmp_unit,*) il, al(il), g(0,1,iglo), g(0,2,iglo)
!          end if
!       end do
!       call close_output_file (tmp_unit)
!       first = .false.
!    end if

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
       ig = ig_idx(lz_lo,ilz)
       ik = ik_idx(lz_lo,ilz)
       is = is_idx(lz_lo,ilz)
       if (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) cycle
       je = 2*jend(ig)

       if (je == 0) then
          je = 2*ng2+1
       end if

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
!
! bug fixed 4.14.03
! was only a problem with collisions but with no trapped particles
! because of an array index out of bounds 
! Overall, this was a rare bug.
!
       if (jend(ig) /= 0) glz(je,ilz) = glz(jend(ig),ilz)

    end do

!    call check_glz ('end', glz)
    call scatter (lorentz_map, glz, g)

!    call check_g ('end', g)

    call prof_leaving ("solfp_lorentz", "collisions")
  end subroutine solfp_lorentz

  subroutine solfp_escatter (g)
    use mp, only: proc0
    use species, only: spec, nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use le_grids, only: negrid, integrate_moment
    use gs2_layouts, only: is_idx, e_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather, scatter
    use file_utils, only: open_output_file, close_output_file
    use egrid, only: zeroes, x0

    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx, il_idx
    use gs2_layouts, only: idx_local, proc_id
    use le_grids, only: forbid
    use mp, only: send, receive

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g

    complex, dimension (negrid) :: delta
    integer :: ielo, ie, is

    call gather (escatter_map, g, gesc)

    ! solve for gesc row by row
    do ielo = e_lo%llim_proc, e_lo%ulim_proc
!       is = is_idx(e_lo,ielo)
       is = 1

       if (spec(is)%vnewk < 2.0*epsilon(0.0)) cycle

       delta(1) = gesc(1,ielo)
       do ie = 1, negrid-1
          delta(ie+1) = gesc(ie+1,ielo) - eql(ie+1)*delta(ie)
       end do
       
       gesc(negrid+1,ielo) = 0.0
       do ie = negrid, 1, -1
          gesc(ie,ielo) = (delta(ie) - ec1(ie)*gesc(ie+1,ielo))*ebetaa(ie)
       end do

    end do

    call scatter (escatter_map, gesc, g)

  end subroutine solfp_escatter

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

