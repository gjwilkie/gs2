module dist_fn
  use init_g, only: ginit, tstart
  use redistribute, only: redist_type
  implicit none

  public :: init_dist_fn
  public :: timeadv
  public :: getfieldeq, getfieldeq0, getan
  public :: flux, neoclassical_flux
  public :: ginit
  public :: init_intcheck, init_vortcheck, init_fieldcheck
  public :: intcheck, vortcheck, fieldcheck
  public :: finish_intcheck, finish_vortcheck, finish_fieldcheck
  public :: t0, omega0, gamma0, thetas, k0, nperiod_guard
  public :: tstart
  public :: reset_init

  private

  ! knobs
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp  ! (nspec)
  real :: gridfac, apfac, driftknob, poisfac
  real :: t0, omega0, gamma0, thetas, k0
  real :: phi_ext, afilter
  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_fieldlineavg = 2, &
       adiabatic_option_yavg = 3
  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_zero = 2, source_option_sine = 3, &
       source_option_test1 = 4, source_option_phiext_full = 5, &
       source_option_test2_full = 6, &
       source_option_convect_full = 7
  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_connected = 3, &
       boundary_option_alternate_zero = 4, &
       boundary_option_conn_via_guards = 5
  logical :: use_shmem_in_connected
  integer :: nperiod_guard

  ! internal arrays

  integer, dimension (:), allocatable :: ittp
  ! (-ntgrid:ntgrid)

  real, dimension (:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:,:,:), allocatable :: wdriftttp
  ! (-ntgrid:ntgrid,naky,ntheta0,negrid,nspec) replicated

  real, dimension (:,:,:), allocatable :: wstar
  ! (naky,negrid,nspec) replicated

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  ! (-ntgrid:ntgrid,naky,ntheta0) replicated

  complex, dimension (:,:), allocatable :: a, b, r, ainv
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:), allocatable :: gridfac1
  ! (-ntgrid:ntgrid,naky,ntheta0)

  complex, dimension (:,:,:), allocatable :: g0
  ! (-ntgrid:ntgrid,2, -g-layout-)

  complex, dimension (:,:,:), allocatable, save :: g_nl
  ! (-ntgrid:ntgrid,2, -g-layout-)

  ! connected bc

  integer, dimension (:,:), allocatable :: itleft, itright
  ! (naky,ntheta0)

  type :: connections_type
     integer :: iproc_left,  iglo_left
     integer :: iproc_right, iglo_right
  end type connections_type

  type (connections_type), dimension (:), allocatable :: connections
  ! (-g-layout-)

  ! connected bc: shmem only
  integer :: shmem_buff_size
  integer, dimension (:), allocatable :: iglo_llim
  ! (0:nproc-1)

  ! connected bc: mp only
  integer :: n_unconnected
  integer, dimension (:), allocatable :: iglo_unconnected
  type :: connected_bc_type
     integer :: n1, n2
     integer, dimension (:), pointer :: iglo1, iglo2
     integer, dimension (:), pointer :: ibc1_in, ibc2_in
     integer, dimension (:), pointer :: ibc1_out, ibc2_out
     integer, dimension (:), pointer :: nin1_start, nin1_end
     integer, dimension (:), pointer :: nin2_start, nin2_end
     integer, dimension (:), pointer :: nout1, nout2
  end type connected_bc_type
  integer :: npass_connected
  type (connected_bc_type), dimension (:), allocatable :: connected
  complex, dimension (:,:), allocatable :: bcleft_out, bcright_out
  complex, dimension (:), allocatable :: bcleft_in, bcright_in

  ! connected-via-guards only
  integer :: lslo, lshi, ldlo, ldhi
  integer :: rslo, rshi, rdlo, rdhi
  type (redist_type) :: gc_from_left, gc_from_right

  type :: guard_connections_type
     integer :: n_left, n_right
     integer, dimension (:), pointer :: iglo_left, iglo_right
  end type guard_connections_type

  type (guard_connections_type), dimension (:), allocatable :: gc
  integer :: ntgrid_noguard

  logical :: initialized = .false.
  logical :: initializing = .true.

contains

  subroutine init_dist_fn
    use mp, only: proc0
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0
    use le_grids, only: init_le_grids, nlambda, negrid
    use run_parameters, only: init_run_parameters
    use collisions, only: init_collisions
    use gs2_layouts, only: init_dist_fn_layouts
    use init_g, only: init_init_g
    implicit none

    if (initialized) return
    initialized = .true.

    call init_species
    call init_theta_grid
    call init_kt_grids
    call init_le_grids
    call init_run_parameters
    call init_collisions
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
    call init_init_g

    call read_parameters
    call allocate_arrays
    call init_vpar
    call init_wdrift
    call init_wstar
    call init_bessel
    call init_invert_rhs
    call init_fieldeq

  end subroutine init_dist_fn

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use theta_grid, only: nperiod
    use init_g, only: init_g_k0 => k0
    use text_options
    use species, only: nspec
    use mp, only: proc0, broadcast
    use shmem, only: shmem_available
    implicit none
    type (text_option), dimension (8), parameter :: sourceopts = &
         (/ text_option('default', source_option_full), &
            text_option('full', source_option_full), &
            text_option('zero', source_option_zero), &
            text_option('sine', source_option_sine), &
            text_option('test1', source_option_test1), &
            text_option('phiext_full', source_option_phiext_full), &
            text_option('test2_full', source_option_test2_full), &
            text_option('convect_full', source_option_convect_full) /)
    character(20) :: source_option

    type (text_option), dimension (8), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('kperiod=1', boundary_option_self_periodic), &
            text_option('connected', boundary_option_connected), &
            text_option('alternate-zero', boundary_option_alternate_zero), &
            text_option('connected-via-guards', &
                        boundary_option_conn_via_guards) /)
    character(20) :: boundary_option

    type (text_option), dimension (6), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg) /)
    character(30) :: adiabatic_option
            
    namelist /dist_fn_knobs/ boundary_option, gridfac, apfac, driftknob, &
         use_shmem_in_connected, nperiod_guard, poisfac, adiabatic_option
    namelist /source_knobs/ t0, omega0, gamma0, &
!         thetas, k0, phi_ext, source_option, thetas, afilter
! D. Ernst 4/8/99: delete 2nd thetas 
           thetas, k0, phi_ext, source_option, afilter
    integer :: ierr

    if (proc0) then
       boundary_option = 'default'
       adiabatic_option = 'default'
       poisfac = 0.0
       gridfac = 5e4
       apfac = 1.0
       driftknob = 1.0
       use_shmem_in_connected = .true.
       nperiod_guard = 1
       t0 = 100.0
       omega0 = 0.0
       gamma0 = 0.0
       thetas = 1.0
       k0 = init_g_k0
       phi_ext = 0.0
       afilter = 0.0
       source_option = 'default'
       read (unit=input_unit("dist_fn_knobs"), nml=dist_fn_knobs)
       read (unit=input_unit("source_knobs"), nml=source_knobs)

       ierr = error_unit()
       call get_option_value &
            (boundary_option, boundaryopts, boundary_option_switch, &
            ierr, "boundary_option in dist_fn_knobs")
       call get_option_value &
            (source_option, sourceopts, source_option_switch, &
            ierr, "source_option in source_knobs")
       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")

       use_shmem_in_connected = use_shmem_in_connected .and. shmem_available

       if (boundary_option_switch == boundary_option_conn_via_guards) then
          if (nperiod_guard > ((nperiod-1)*2+1)/3) then
             nperiod_guard = ((nperiod-1)*2+1)/3
             write (error_unit(), *) "WARNING: nperiod_guard is too large ", &
                  "for nperiod=", nperiod, &
                  ", setting nperiod_guard to ", nperiod_guard
             print *, "WARNING: nperiod_guard is too large ", &
                  "for nperiod=", nperiod, &
                  ", setting nperiod_guard to ", nperiod_guard
          end if

          if (nperiod_guard <= 0) then
             boundary_option_switch = boundary_option_zero
          end if
       else
          nperiod_guard = 0
       end if
    end if
    if (.not.allocated(fexp)) allocate (fexp(nspec), bkdiff(nspec), bd_exp(nspec))
    if (proc0) call read_species_knobs

    call broadcast (boundary_option_switch)
    call broadcast (adiabatic_option_switch)
    call broadcast (gridfac)
    call broadcast (poisfac)
    call broadcast (apfac)
    call broadcast (driftknob)
    call broadcast (t0)
    call broadcast (omega0)
    call broadcast (gamma0)
    call broadcast (thetas)
    call broadcast (k0)
    call broadcast (phi_ext)
    call broadcast (afilter)
    call broadcast (source_option_switch)
    call broadcast (fexp)
    call broadcast (bkdiff)
    call broadcast (bd_exp)
    call broadcast (use_shmem_in_connected)
    call broadcast (nperiod_guard)
  end subroutine read_parameters

  subroutine read_species_knobs
    use species, only: nspec
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer :: is, unit

    do is = 1, nspec
       fexp(is) = (0.4,0.0)
       bkdiff(is) = 0.0
       bd_exp(is) = 0
       call get_indexed_namelist_unit (unit, "dist_fn_species_knobs", is)
       call fill_species_knobs (unit, fexp(is), bkdiff(is), bd_exp(is))
       close (unit=unit)
    end do
  end subroutine read_species_knobs

  subroutine fill_species_knobs (unit, fexp_out, bakdif_out, bd_exp_out)
    implicit none
    integer, intent (in) :: unit
    complex, intent (in out) :: fexp_out
    real, intent (in out) :: bakdif_out
    integer, intent (in out) :: bd_exp_out
    real :: fexpr, fexpi, bakdif, bd_exp
    namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

    fexpr = real(fexp_out)
    fexpi = aimag(fexp_out)
    bakdif = bakdif_out
    bd_exp = bd_exp_out
    read (unit=unit, nml=dist_fn_species_knobs)
    fexp_out = cmplx(fexpr,fexpi)
    bakdif_out = bakdif
    bd_exp_out = bd_exp
  end subroutine fill_species_knobs

  subroutine init_wdrift
    use species, only: nspec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda, e, al, jend
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo
    logical :: alloc = .true.

    if (alloc) allocate (ittp(-ntgrid:ntgrid))
    ittp = 0
    do ig = -ntgrid+1, ntgrid-1
       if (jend(ig) > 0 .and. jend(ig) <= nlambda) then
          if (1.0-al(jend(ig))*bmag(ig+1) < 2.0*epsilon(0.0) &
               .and. 1.0-al(jend(ig))*bmag(ig-1) < 2.0*epsilon(0.0)) &
          then
             ittp(ig) = jend(ig)
          end if
       end if
    end do

    if (alloc) allocate (wdrift(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    if (alloc) allocate (wdriftttp(-ntgrid:ntgrid,naky,ntheta0,negrid,nspec))
    wdrift = 0.  ; wdriftttp = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid
          wdrift(ig,iglo) &
               = wdrift_func(ig, ik_idx(g_lo,iglo), it_idx(g_lo,iglo), &
                                 il_idx(g_lo,iglo), ie_idx(g_lo,iglo), &
                                 is_idx(g_lo,iglo))
       end do
    end do
    wdriftttp = 0.0
    do ig = -ntgrid, ntgrid
       if (ittp(ig) == 0) cycle
       do ik = 1, naky
          do it = 1, ntheta0
             do ie = 1, negrid
                do is = 1, nspec
                   wdriftttp(ig,ik,it,ie,is) &
                        = wdrift_func(ig,ik,it,ittp(ig),ie,is)*driftknob
                end do
             end do
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid-1
          wdrift(ig,iglo) = 0.5*(wdrift(ig,iglo) + wdrift(ig+1,iglo))*driftknob
       end do
    end do

    alloc = .false.

  end subroutine init_wdrift

  function wdrift_func (ig, ik, it, il, ie, is)
    use theta_grid, only: theta, bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: shat
    use kt_grids, only: aky, theta0, akx
    use le_grids, only: e, al
    use run_parameters, only: delt, wunits
    implicit none
    real :: wdrift_func
    integer, intent (in) :: ig, ik, it, il, ie, is

    if (aky(ik) == 0.0) then
       wdrift_func = akx(it)/shat &
                    *(cvdrift0(ig)*e(ie,is)*(1.0 - al(il)*bmag(ig)) &
                      + gbdrift0(ig)*0.5*e(ie,is)*al(il)*bmag(ig)) &
                     *delt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(ik,it)*cvdrift0(ig)) &
                        *e(ie,is)*(1.0 - al(il)*bmag(ig)) &
                      + (gbdrift(ig) + theta0(ik,it)*gbdrift0(ig)) &
                        *0.5*e(ie,is)*al(il)*bmag(ig)) &
                     *delt*wunits(ik)
    end if
  end function wdrift_func

  subroutine init_vpar
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2
    use species, only: spec, zstm
    use theta_grid, only: ntgrid, delthet, bmag, gradpar
    use kt_grids, only: aky
    use le_grids, only: e, al
    use run_parameters, only: delt, tunits
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo, ik, is
    real :: al1, e1
    

    if (.not.allocated(vpa)) then
       allocate (vpa(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vperp2(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpar(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    vpa = 0. ; vpac = 0. ; vperp2 = 0. ; vpar = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       al1 = al(il_idx(g_lo,iglo))
       e1 = e(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))

       vpa(:,1,iglo) = sqrt(e1*max(0.0, 1.0 - al1*bmag))
       vpa(:,2,iglo) = - vpa(:,1,iglo)
       vperp2(:,iglo) = bmag*al1*e1

       where (1.0 - al1*bmag < 100.0*epsilon(0.0))
          vpa(:,1,iglo) = 0.0
          vpa(:,2,iglo) = 0.0
       end where

       where (1.0 - al1*0.5*(bmag(-ntgrid:ntgrid-1)+bmag(-ntgrid+1:ntgrid)) &
              < 0.0)
          vpac(-ntgrid:ntgrid-1,1,iglo) = 1.0
          vpac(-ntgrid:ntgrid-1,2,iglo) = -1.0
       elsewhere
          vpac(-ntgrid:ntgrid-1,1,iglo) = &
              0.5*(vpa(-ntgrid:ntgrid-1,1,iglo) + vpa(-ntgrid+1:ntgrid,1,iglo))
          vpac(-ntgrid:ntgrid-1,2,iglo) = &
              0.5*(vpa(-ntgrid:ntgrid-1,2,iglo) + vpa(-ntgrid+1:ntgrid,2,iglo))
       end where
       vpac(ntgrid,:,iglo) = 0.0

       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       vpar(-ntgrid:ntgrid-1,1,iglo) = &
            zstm(spec,is)*tunits(ik)*delt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,1,iglo)
       vpar(-ntgrid:ntgrid-1,2,iglo) = &
            zstm(spec,is)*tunits(ik)*delt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,2,iglo)
       vpar(ntgrid,:,iglo) = 0.0
    end do
  end subroutine init_vpar

  subroutine init_wstar
    use species, only: nspec, spec
    use kt_grids, only: naky
    use le_grids, only: negrid, e
    use run_parameters, only: delt, wunits
    implicit none
    integer :: ik, ie, is

    if(.not.allocated(wstar)) allocate (wstar(naky,negrid,nspec))

    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             wstar(ik,ie,is) = delt*wunits(ik) &
                  *(spec(is)%fprim+spec(is)%tprim*(e(ie,is)-1.5))
          end do
       end do
    end do
  end subroutine init_wstar

  subroutine init_bessel
    use dist_fn_arrays, only: aj0, aj1, kperp2, aj0f, aj1f
    use species, only: nspec, spec, smz
    use theta_grid, only: ntgrid, bmag, gds2, gds21, gds22, shat
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    use le_grids, only: e, al
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo
    real :: arg
    logical :: done = .false.

    if (done) return
    done = .true.

    allocate (kperp2(-ntgrid:ntgrid,naky,ntheta0))
    do it = 1, ntheta0
       do ik = 1, naky
          if (aky(ik) == 0.0) then
             kperp2(:,ik,it) = akx(it)*akx(it)*gds22/(shat*shat)
          else
             kperp2(:,ik,it) = aky(ik)*aky(ik) &
                  *(gds2 + 2.0*theta0(ik,it)*gds21 &
                    + theta0(ik,it)*theta0(ik,it)*gds22)
          end if
       end do
    end do

    allocate (aj0(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj0f(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1f(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    aj0 = 0. ; aj1 = 0.; aj0f = 0. ; aj1f = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          arg = smz(spec,is)*sqrt(e(ie,is)*al(il)/bmag(ig)*kperp2(ig,ik,it))
          aj0(ig,iglo) = j0(arg)
          aj1(ig,iglo) = j1(arg)
          aj0f(ig,iglo) = j0(arg)*exp(-afilter**2*kperp2(ig,ik,it))
          aj1f(ig,iglo) = j1(arg)*exp(-afilter**2*kperp2(ig,ik,it))
       end do
    end do

  end subroutine init_bessel

  subroutine init_invert_rhs
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2
    use species, only: nspec, spec, tz
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda, ng2, forbid
    use constants
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is
    real :: wd, wdttp, vp
    real, dimension (-ntgrid:ntgrid) :: th_tmp

    if (.not.allocated(a)) then
       allocate (a(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (b(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (r(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (ainv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    a = 0. ; b = 0. ; r = 0. ; ainv = 0.
    
    th_tmp(-ntgrid) = -pi/2.
    do ig = -ntgrid+1, ntgrid
       th_tmp(ig) = th_tmp(ig-1) + pi/real(2*ntgrid)
    enddo
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid-1
          wd = wdrift(ig,iglo)
          wdttp = wdriftttp(ig,ik,it,ie,is)
          vp = vpar(ig,1,iglo)
          
          ainv(ig,iglo) &
               = 1.0/(1.0 + bkdiff(is)*sin(th_tmp(ig))**bd_exp(is) &
               + (1.0-fexp(is))*tz(spec,is)*(zi*wd + 2.0*vp))
          r(ig,iglo) &
               = (1.0 - bkdiff(is)*sin(th_tmp(ig))**bd_exp(is) &
               + (1.0-fexp(is))*tz(spec,is)*(zi*wd - 2.0*vp)) &
               *ainv(ig,iglo)
          a(ig,iglo) &
               = 1.0 + bkdiff(is)*sin(th_tmp(ig))**bd_exp(is) &
               + fexp(is)*tz(spec,is)*(-zi*wd - 2.0*vp)
          b(ig,iglo) &
               = 1.0 - bkdiff(is)*sin(th_tmp(ig))**bd_exp(is) &
               + fexp(is)*tz(spec,is)*(-zi*wd + 2.0*vp)
          
          if (nlambda > ng2) then
             ! zero out forbidden regions
             if (forbid(ig,il) .or. forbid(ig+1,il)) then
                r(ig,iglo) = 0.0
                ainv(ig,iglo) = 0.0
             end if
             
             ! ???? mysterious mucking around at lower bounce point
             if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
                ainv(ig,iglo) = 1.0 + ainv(ig,iglo)
             end if
             
             ! ???? mysterious mucking around with totally trapped particles
             if (il == ittp(ig)) then
                ainv(ig,iglo) = 1.0/(1.0 + zi*(1.0-fexp(is))*tz(spec,is)*wdttp)
                a(ig,iglo) = 1.0 - zi*fexp(is)*tz(spec,is)*wdttp
                r(ig,iglo) = 0.0
             end if
          end if
       end do
    end do

    select case (boundary_option_switch)
    case (boundary_option_connected,boundary_option_conn_via_guards)
       call init_connected_bc
    case default
       !nothing
    end select

    initializing = .false.
  end subroutine init_invert_rhs

  subroutine init_connected_bc
    use file_utils, only: error_unit
    use theta_grid, only: ntgrid, nperiod, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use le_grids, only: ng2, nlambda
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id, idx_local
    use mp, only: proc0, iproc, nproc, broadcast, max_allreduce
    use constants
    use redistribute, only: index_list_type, init_fill, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_left, from_right
    type (index_list_type), dimension(0:nproc-1) :: to_right, from_left
    integer, dimension (0:nproc-1) :: nn_from_left, nn_to_right
    integer, dimension (0:nproc-1) :: nn_from_right, nn_to_left
    integer, dimension (3) :: to_low, from_low
    integer :: ik, it, il, ie, is, iglo, it0, itl, itr, jshift0
    integer :: ip, ipleft, ipright, ipass, ipass_1, ipass_2
    integer :: iglo_left, iglo_right, ipp, i
    integer :: n, n1, n2, isign, ig
    integer, dimension (:,:), allocatable :: nbcin1, nbcin2
    logical :: done = .false.

    if (done) return
    done = .true.

    ntgrid_noguard = ntheta/2 + (nperiod-nperiod_guard-1)*ntheta

    if (naky > 1 .and. ntheta0 > 1) then
       jshift0 = int((theta(ntgrid_noguard)-theta(-ntgrid_noguard)) &
                     /(theta0(2,2)-theta0(2,1)) + 0.1)
    else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) /= 0.0) then
       jshift0 = int((theta(ntgrid_noguard)-theta(-ntgrid_noguard)) &
                     /(theta0(1,2)-theta0(1,1)) + 0.1)
    else
       jshift0 = 1
    end if

    allocate (itleft(naky,ntheta0), itright(naky,ntheta0))
    itleft(1,:) = -1
    itright(1,:) = -1
    do ik = 1, naky
       do it = 1, ntheta0
          if (it > (ntheta0+1)/2) then
             it0 = it - ntheta0 - 1
          else
             it0 = it - 1
          end if

          if (ik == 1) then
             if (aky(ik) /= 0.0 .and. naky == 1) then
                itl = it0 + jshift0
                itr = it0 - jshift0
             else
                itl = ntheta0
                itr = ntheta0
             end if
          else
             itl = it0 + (ik-1)*jshift0
             itr = it0 - (ik-1)*jshift0
          end if

          if (itl >= 0 .and. itl < (ntheta0+1)/2) then
             itleft(ik,it) = itl + 1
          else if (itl + ntheta0 + 1 > (ntheta0+1)/2 &
             .and. itl + ntheta0 + 1 <= ntheta0) then
             itleft(ik,it) = itl + ntheta0 + 1
          else
             itleft(ik,it) = -1
          end if

          if (itr >= 0 .and. itr < (ntheta0+1)/2) then
             itright(ik,it) = itr + 1
          else if (itr + ntheta0 + 1 > (ntheta0+1)/2 &
             .and. itr + ntheta0 + 1 <= ntheta0) then
             itright(ik,it) = itr + ntheta0 + 1
          else
             itright(ik,it) = -1
          end if
       end do
    end do

    allocate (connections(g_lo%llim_proc:g_lo%ulim_alloc))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       if (nlambda > ng2 .and. il >= ng2+2) then
          connections(iglo)%iproc_left = -1
          connections(iglo)%iglo_left = -1
          connections(iglo)%iproc_right = -1
          connections(iglo)%iglo_right = -1
       else
          if (itleft(ik,it) < 0) then
             connections(iglo)%iproc_left = -1
             connections(iglo)%iglo_left = -1
          else
             connections(iglo)%iproc_left &
                  = proc_id(g_lo,idx(g_lo,ik,itleft(ik,it),il,ie,is))
             connections(iglo)%iglo_left &
                  = idx(g_lo,ik,itleft(ik,it),il,ie,is)
          end if
          if (itright(ik,it) < 0) then
             connections(iglo)%iproc_right = -1
             connections(iglo)%iglo_right = -1
          else
             connections(iglo)%iproc_right &
                  = proc_id(g_lo,idx(g_lo,ik,itright(ik,it),il,ie,is))
             connections(iglo)%iglo_right &
                  = idx(g_lo,ik,itright(ik,it),il,ie,is)
          end if
       end if
    end do

    select case (boundary_option_switch)
    case (boundary_option_connected)
       if (use_shmem_in_connected) then
          allocate (iglo_llim(0:nproc-1))
          iglo_llim(iproc) = g_lo%llim_proc
          do ip = 0, nproc-1
             call broadcast (iglo_llim(ip), ip)
          end do
          shmem_buff_size = g_lo%ulim_proc - g_lo%llim_proc + 1
          call max_allreduce (shmem_buff_size)
       else
          ! get n_unconnected, iglo_unconnected
          n_unconnected = 0
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             if (itleft(ik,it) < 0 .and. itright(ik,it) < 0) then
                n_unconnected = n_unconnected + 1
             end if
          end do
          if (n_unconnected > 0) allocate (iglo_unconnected(n_unconnected))
          n_unconnected = 0
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             if (itleft(ik,it) < 0 .and. itright(ik,it) < 0) then
                n_unconnected = n_unconnected + 1
                iglo_unconnected(n_unconnected) = iglo
             end if
          end do

          ! get npass_nconnected
          npass_connected = 1
          do iglo = g_lo%llim_world, g_lo%ulim_world
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             if (itleft(ik,it) < 0) then
                call get_ipass_2 (iglo, ipass_2)
                npass_connected = max(npass_connected, ipass_2)
             end if
          end do

          allocate (connected(npass_connected))
          do ipass = 1, npass_connected
             connected(ipass)%n1 = 0
             connected(ipass)%n2 = 0
             allocate (connected(ipass)%nin1_start(0:nproc-1))
             allocate (connected(ipass)%nin1_end(0:nproc-1))
             connected(ipass)%nin1_start = 1
             connected(ipass)%nin1_end = 0
             allocate (connected(ipass)%nin2_start(0:nproc-1))
             allocate (connected(ipass)%nin2_end(0:nproc-1))
             connected(ipass)%nin2_start = 1
             connected(ipass)%nin2_end = 0
             allocate (connected(ipass)%nout1(0:nproc-1))
             connected(ipass)%nout1 = 0
             allocate (connected(ipass)%nout2(0:nproc-1))
             connected(ipass)%nout2 = 0
          end do

          ! get connected(:)%n1, connected(:)%n2
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             call get_left_connection (iglo, iglo_left, ipleft)
             if (ipleft < 0) then
                ipass_1 = 1
                connected(ipass_1)%n1 = connected(ipass_1)%n1 + 1
             else if (ipleft /= iproc) then
                call get_ipass_1 (iglo, ipass_1)
                connected(ipass_1)%n1 = connected(ipass_1)%n1 + 1
                connected(ipass_1)%nin1_end(ipleft) &
                     = connected(ipass_1)%nin1_end(ipleft) + 1
             end if

             call get_right_connection (iglo, iglo_right, ipright)
             if (ipright < 0) then
                ipass_2 = 1
                connected(ipass_2)%n2 = connected(ipass_2)%n2 + 1
             else if (ipright /= iproc) then
                call get_ipass_2 (iglo, ipass_2)
                connected(ipass_2)%n2 = connected(ipass_2)%n2 + 1
                connected(ipass_2)%nin2_end(ipright) &
                     = connected(ipass_2)%nin2_end(ipright) + 1
             end if
          end do

          do ipass = 1, npass_connected
             n1 = connected(ipass)%n1
             n2 = connected(ipass)%n2
             allocate (connected(ipass)%iglo1(n1))
             allocate (connected(ipass)%iglo2(n2))
             allocate (connected(ipass)%ibc1_in(n1))
             allocate (connected(ipass)%ibc2_in(n2))
             allocate (connected(ipass)%ibc1_out(n1))
             allocate (connected(ipass)%ibc2_out(n2))
             connected(ipass)%n1 = 0
             connected(ipass)%n2 = 0
          end do

          ! get connected(:)%iglo1(:), connected(:)%iglo2(:)
          ! get connected(:)%ibc1_out(:), connected(:)%ibc2_out(:)
          ! get connected(:)%nout1(:), connected(:)%nout2(:)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             call get_left_connection (iglo, iglo_left, ipleft)
             call get_right_connection (iglo, iglo_right, ipright)

             if (ipleft /= iproc) then
                if (ipleft < 0) then
                   ipass_1 = 1
                else
                   call get_ipass_1 (iglo, ipass_1)
                end if
                n1 = connected(ipass_1)%n1 + 1
                connected(ipass_1)%n1 = n1
                connected(ipass_1)%iglo1(n1) = iglo
                connected(ipass_1)%nout1(ipright) &
                     = connected(ipass_1)%nout1(ipright) + 1
                connected(ipass_1)%ibc1_out(n1) &
                     = connected(ipass_1)%nout1(ipright)
             end if

             if (ipright /= iproc) then
                if (ipright < 0) then
                   ipass_2 = 1
                else
                   call get_ipass_2 (iglo, ipass_2)
                end if
                n2 = connected(ipass_2)%n2 + 1
                connected(ipass_2)%n2 = n2
                connected(ipass_2)%iglo2(n2) = iglo
                connected(ipass_2)%nout2(ipleft) &
                     = connected(ipass_2)%nout2(ipleft) + 1
                connected(ipass_2)%ibc2_out(n2) &
                     = connected(ipass_2)%nout2(ipleft)
             end if
          end do

          ! get connected(:)%nin1_start(:), connected%nin1_end(:)
          ! get connected(:)%nin2_start(:), connected%nin2_end(:)
          do iglo = g_lo%llim_world, g_lo%ulim_world
             ip = proc_id(g_lo,iglo)
             if (ip == iproc) cycle

             call get_left_connection (iglo, iglo_left, ipleft)
             call get_right_connection (iglo, iglo_right, ipright)

             if (ipleft == iproc) then
                call get_ipass_1 (iglo_left, ipass_1)
                connected(ipass_1)%nin1_end(ip) &
                     = connected(ipass_1)%nin1_end(ip) + 1
             end if

             if (ipright == iproc) then
                call get_ipass_2 (iglo_right, ipass_2)
                connected(ipass_2)%nin2_end(ip) &
                     = connected(ipass_2)%nin2_end(ip) + 1
             end if
          end do
          do ipass = 1, npass_connected
             do ip = 1, nproc-1
                connected(ipass)%nin1_start(ip) &
                     = connected(ipass)%nin1_end(ip-1) + 1
                connected(ipass)%nin1_end(ip) &
                     = connected(ipass)%nin1_end(ip) &
                       + connected(ipass)%nin1_start(ip) - 1
                connected(ipass)%nin2_start(ip) &
                     = connected(ipass)%nin2_end(ip-1) + 1
                connected(ipass)%nin2_end(ip) &
                     = connected(ipass)%nin2_end(ip) &
                       + connected(ipass)%nin2_start(ip) - 1
             end do
          end do

          ! get connected(:)%ibc1_in(:), connected(:)%ibc2_in(:)
          allocate (nbcin1(npass_connected,0:nproc-1))
          allocate (nbcin2(npass_connected,0:nproc-1))
          do ipass = 1, npass_connected
             nbcin1(ipass,:) = connected(ipass)%nin1_start(:) - 1
             nbcin2(ipass,:) = connected(ipass)%nin2_start(:) - 1
          end do
          do iglo = g_lo%llim_world, g_lo%ulim_world
             ip = proc_id(g_lo,iglo)
             if (ip == iproc) cycle

             call get_left_connection (iglo, iglo_left, ipleft)
             call get_right_connection (iglo, iglo_right, ipright)

             if (ipleft == iproc) then
                call get_ipass_1 (iglo_left, ipass_1)
                nbcin1(ipass_1,ip) = nbcin1(ipass_1,ip) + 1
                do
                   call get_left_connection (iglo_left, i, ipp)
                   if (ipp /= ipleft) exit
                   iglo_left = i
                end do
                do i = 1, connected(ipass_1)%n1
                   if (iglo_left == connected(ipass_1)%iglo1(i)) then
                      connected(ipass_1)%ibc1_in(i) = nbcin1(ipass_1,ip)
                      exit
                   end if
                end do
             end if

             if (ipright == iproc) then
                call get_ipass_2 (iglo_right, ipass_2)
                nbcin2(ipass_2,ip) = nbcin2(ipass_2,ip) + 1
                do
                   call get_right_connection (iglo_right, i, ipp)
                   if (ipp /= ipright) exit
                   iglo_right = i
                end do
                do i = 1, connected(ipass_2)%n2
                   if (iglo_right == connected(ipass_2)%iglo2(i)) then
                      connected(ipass_2)%ibc2_in(i) = nbcin2(ipass_2,ip)
                   end if
                end do
             end if
          end do

          ! the following doesn't work:
          !allocate (bcleft_out(maxval(connected%nout1),0:nproc-1))
          !allocate (bcright_out(maxval(connected%nout2),0:nproc-1))
          ! replace with the following:
          n1 = 0
          n2 = 0
          do ipass = 1, npass_connected
             n1 = max(n1,maxval(connected(ipass)%nout1))
             n2 = max(n2,maxval(connected(ipass)%nout2))
          end do
          allocate (bcleft_out(n1,0:nproc-1))
          allocate (bcright_out(n2,0:nproc-1))
          
          allocate (bcleft_in(maxval(nbcin1)))
          allocate (bcright_in(maxval(nbcin2)))
          
          deallocate (nbcin1, nbcin2)
       end if
    case (boundary_option_conn_via_guards)
       ldlo = -ntgrid
       ldhi = -ntgrid + ntheta*nperiod_guard - 1
       lslo = -ntgrid + ntheta*nperiod_guard + 1
       lshi = -ntgrid + 2*ntheta*nperiod_guard

       rdhi = ntgrid
       rdlo = ntgrid - ntheta*nperiod_guard + 1
       rshi = ntgrid - ntheta*nperiod_guard - 1
       rslo = ntgrid - 2*ntheta*nperiod_guard

       nn_from_left = 0
       nn_to_right = 0

       nn_from_right = 0
       nn_to_left = 0

       do iglo = g_lo%llim_world, g_lo%ulim_world
          ip = proc_id(g_lo,iglo)
          call get_left_connection (iglo, iglo_left, ipleft)
          call get_right_connection (iglo, iglo_right, ipright)

          if (ip == iproc .and. ipleft >= 0) then
             do isign = 1, 2
                do ig = lslo, lshi
                   nn_from_left(ipleft) = nn_from_left(ipleft) + 1
                end do
             end do
          end if
          if (ipleft == iproc) then 
             do isign = 1, 2
                do ig = rdlo, rdhi
                   nn_to_right(ip) = nn_to_right(ip) + 1
                end do
             end do
          endif

          if (ip == iproc .and. ipright >= 0 .and. ipright /= iproc) then
             do isign = 1, 2
                do ig = rslo, rshi
                   nn_from_right(ipright) = nn_from_right(ipright) +1
                end do
             end do
          end if
          if (ipright == iproc .and. ip /= iproc) then
             do isign = 1, 2
                do ig = ldlo, ldhi
                   nn_to_left(ip) = nn_to_left(ip) + 1
                end do
             end do
          end if
       end do

       do ip = 0, nproc-1
          if (nn_from_left(ip) > 0) then
             allocate (from_left(ip)%first(nn_from_left(ip)))
             allocate (from_left(ip)%second(nn_from_left(ip)))
             allocate (from_left(ip)%third(nn_from_left(ip)))
          endif
          if (nn_to_right(ip) > 0) then
             allocate (to_right(ip)%first(nn_to_right(ip)))
             allocate (to_right(ip)%second(nn_to_right(ip)))
             allocate (to_right(ip)%third(nn_to_right(ip)))
          endif
          if (nn_from_right(ip) > 0) then
             allocate (from_right(ip)%first(nn_from_right(ip)))
             allocate (from_right(ip)%second(nn_from_right(ip)))
             allocate (from_right(ip)%third(nn_from_right(ip)))
          endif
          if (nn_to_left(ip) > 0) then
             allocate (to_left(ip)%first(nn_to_left(ip)))
             allocate (to_left(ip)%second(nn_to_left(ip)))
             allocate (to_left(ip)%third(nn_to_left(ip)))
          endif
       end do

       nn_from_left = 0
       nn_to_right = 0

       nn_from_right = 0
       nn_to_left = 0

       do iglo = g_lo%llim_world, g_lo%ulim_world
          ip = proc_id(g_lo,iglo)
          call get_left_connection (iglo, iglo_left, ipleft)
          call get_right_connection (iglo, iglo_right, ipright)

          if (ip == iproc .and. ipleft >= 0) then
             do isign = 1, 2
                do ig = lslo, lshi
                   n = nn_from_left(ipleft) + 1
                   nn_from_left(ipleft) = n
                   from_left(ipleft)%first(n) = ig
                   from_left(ipleft)%second(n) = isign
                   from_left(ipleft)%third(n) = iglo
                end do
             end do
          end if
          if (ipleft == iproc) then
             do isign = 1, 2
                do ig = rdlo, rdhi
                   n = nn_to_right(ip) + 1
                   nn_to_right(ip) = n
                   to_right(ip)%first(n) = ig
                   to_right(ip)%second(n) = isign
                   to_right(ip)%third(n) = iglo_left
                end do
             end do
          end if

          if (ip == iproc .and. ipright >= 0 .and. ipright /= iproc) then
             do isign = 1, 2
                do ig = rslo, rshi
                   n = nn_from_right(ipright) + 1
                   nn_from_right(ipright) = n
                   from_right(ipright)%first(n) = ig
                   from_right(ipright)%second(n) = isign
                   from_right(ipright)%third(n) = iglo
                end do
             end do
          end if
          if (ipright == iproc .and. ip /= iproc) then
             do isign = 1, 2
                do ig = ldlo, ldhi
                   n = nn_to_left(ip) + 1
                   nn_to_left(ip) = n
                   to_left(ip)%first(n) = ig
                   to_left(ip)%second(n) = isign
                   to_left(ip)%third(n) = iglo_right
                end do
             end do
          end if
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       to_low (1) = -ntgrid
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       call init_fill ( gc_from_left, 'c', to_low, to_right, from_low, from_left)
       call init_fill (gc_from_right, 'c', to_low, to_left, from_low, from_right)
       
       call delete_list (from_left)
       call delete_list (to_right)
       
       call delete_list (from_right)
       call delete_list (to_left)

    end select

  end subroutine init_connected_bc

  subroutine get_left_connection (iglo, iglo_left, iproc_left)
    use gs2_layouts, only: g_lo, proc_id, idx
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_left, iproc_left
    integer :: ik, it, il, ie, is

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    
    if (itleft(ik,it) < 0) then
       iglo_left = -1
       iproc_left = -1
       return
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    iglo_left = idx(g_lo,ik,itleft(ik,it),il,ie,is)
    iproc_left = proc_id(g_lo,iglo_left)
  end subroutine get_left_connection

  subroutine get_right_connection (iglo, iglo_right, iproc_right)
    use gs2_layouts, only: g_lo, proc_id, idx
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_right, iproc_right
    integer :: ik, it, il, ie, is

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    
    if (itright(ik,it) < 0) then
       iglo_right = -1
       iproc_right = -1
       return
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    iglo_right = idx(g_lo,ik,itright(ik,it),il,ie,is)
    iproc_right = proc_id(g_lo,iglo_right)
  end subroutine get_right_connection

  subroutine get_ipass_1 (iglo, ipass_1)
    use gs2_layouts, only: g_lo, proc_id
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: ipass_1
    integer :: iglo1, iglo1next, ip, ipnext

    iglo1 = iglo
    ipass_1 = 1
    ip = proc_id(g_lo,iglo)
    do
       call get_left_connection (iglo1, iglo1next, ipnext)
       if (ipnext < 0) return
       if (ipnext /= ip) then
          ip = ipnext
          ipass_1 = ipass_1 + 1
       end if
       iglo1 = iglo1next
    end do
  end subroutine get_ipass_1

  subroutine get_ipass_2 (iglo, ipass_2)
    use gs2_layouts, only: g_lo, proc_id
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: ipass_2
    integer :: iglo1, iglo1next, ip, ipnext

    iglo1 = iglo
    ipass_2 = 1
    ip = proc_id(g_lo,iglo)
    do
       call get_right_connection (iglo1, iglo1next, ipnext)
       if (ipnext < 0) return
       if (ipnext /= ip) then
          ip = ipnext
          ipass_2 = ipass_2 + 1
       end if
       iglo1 = iglo1next
    end do
  end subroutine get_ipass_2

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo
    implicit none
    logical :: alloc = .true.

    if (alloc) then
       allocate (g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g_nl(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif

    g = 0. ; gnew = 0. ; g0 = 0.;  g_nl = 0.

    alloc = .false.
  end subroutine allocate_arrays

  subroutine timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep, dt_cfl)
    use theta_grid, only: ntgrid
    use collisions, only: solfp1
    use dist_fn_arrays, only: gnew
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    real :: dt_cfl

    call invert_rhs (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
    call solfp1 (gnew, g0, phinew, aparnew, aperpnew)
  end subroutine timeadv

  subroutine get_source_term &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, &
        isgn, iglo, sourcefac, source)
    use dist_fn_arrays, only: aj0, aj1, vperp2, vpar, vpac, g
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: aky, theta0
    use le_grids, only: nlambda, ng2, lmax, anon, e, negrid
    use species, only: spec, zstm, stm, tz, nspec
    use run_parameters, only: fphi, fapar, faperp, wunits, delt, tunits
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    integer, intent (in) :: isgn, iglo
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:), intent (out) :: source

    integer :: ig, ik, it, il, ie, is
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    phigavg  = (fexp(is)*phi(:,ik,it)   + (1.0-fexp(is))*phinew(:,ik,it)) &
                *aj0(:,iglo)*fphi &
             + (fexp(is)*aperp(:,ik,it) + (1.0-fexp(is))*aperpnew(:,ik,it))&
                *aj1(:,iglo)*faperp*2.0*vperp2(:,iglo)*tz(spec,is)
    apargavg = (fexp(is)*apar(:,ik,it)  + (1.0-fexp(is))*aparnew(:,ik,it)) &
                *aj0(:,iglo)*fapar

    ! source term in finite difference equations
    select case (source_option_switch)
    case (source_option_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if       
    case(source_option_phiext_full)
       if (il <= lmax) then
          call set_source
          if (istep > 0 .and. aky(ik) < epsilon(0.0)) then
             source(:ntgrid-1) = source(:ntgrid-1) &
                  - zi*anon(ie,is)*wdrift(:ntgrid-1,iglo)*2.0*phi_ext*sourcefac
          end if
       else
          source = 0.0
       end if
    case(source_option_test2_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if
       if (aky(ik) == 0.0 .and. istep > 0) then
          source(:ntgrid-1) = source(:ntgrid-1) &
               + tunits(ik)*delt &
               *(aj0(:ntgrid-1,iglo)**2 + aj0(-ntgrid+1:,iglo)**2) &
               *(e(ie,is)-1.5)*sourcefac &
               *exp(-(theta(:ntgrid-1)/thetas)**2)
       end if
    case(source_option_convect_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if
       if (aky(ik) /= 0.0 .and. istep > 0) then
          source(:ntgrid-1) = sourcefac &
               *exp(cmplx((theta(:ntgrid-1)-theta0(ik,it))**2, &
               k0*theta(:ntgrid-1)))
       end if       
    case (source_option_zero)
       source = 0.0
    case (source_option_sine)
       source(:ntgrid-1) = tunits(ik)*delt &
            *(sin(theta(:ntgrid-1)) + sin(theta(-ntgrid+1:)))
    case (source_option_test1)
       do ig = -ntgrid, ntgrid-1
          source(ig) = -zi*(wdrift_func(ig,ik,it,il,ie,is) &
               + wdrift_func(ig+1,ik,it,il,ie,is)) &
               *sourcefac
       end do
    end select

    if (il <= lmax) then
       if (isgn == 1) then
          do ig = -ntgrid, ntgrid-1
             source(ig) = source(ig) &
                  + b(ig,iglo)*g(ig,1,iglo) + a(ig,iglo)*g(ig+1,1,iglo)
          end do
       else
          do ig = -ntgrid, ntgrid-1
             source(ig) = source(ig) &
                  + a(ig,iglo)*g(ig,2,iglo) + b(ig,iglo)*g(ig+1,2,iglo)
          end do
       end if
    end if

    source(ntgrid) = source(-ntgrid)

    ! special source term for totally trapped particles
    if (nlambda > ng2 .and. isgn == 2) then
       do ig = -ntgrid, ntgrid
          if (il /= ittp(ig)) cycle
          source(ig) &
               = g(ig,2,iglo)*a(ig,iglo) &
                 - anon(ie,is)*zi*wdriftttp(ig,ik,it,ie,is)*phigavg(ig) &
                 + zi*wstar(ik,ie,is)*phigavg(ig)
       end do
       if (source_option_switch == source_option_phiext_full .and.  &
            aky(ik) < epsilon(0.0) .and. istep > 0) then
          do ig = -ntgrid, ntgrid
             if (il /= ittp(ig)) cycle             
             source(ig) = source(ig) - zi*anon(ie,is)* &
                  wdriftttp(ig,ik,it,ie,is)*phi_ext*sourcefac          
          end do
       endif
    end if

  contains

    subroutine set_source
      complex :: apar_p, apar_m, phi_p, phi_m
      real, dimension(:,:), allocatable, save :: ufac
      integer :: i_e, i_s
      logical :: first = .true.

      if (first) then
         first = .false.
         allocate (ufac(negrid, nspec))
         do i_e = 1, negrid
            do i_s = 1, nspec
               ufac(i_e, i_s) = (2.0*spec(i_s)%uprim &
                    + spec(i_s)%uprim2*e(i_e,i_s)**(1.5)*sqrt(pi)/4.0)
            end do
         end do
      endif

      do ig = -ntgrid, ntgrid-1
         phi_m = phigavg(ig+1)-phigavg(ig)
         phi_p = phigavg(ig+1)+phigavg(ig)
         apar_p = apargavg(ig+1)+apargavg(ig)
         apar_m = aparnew(ig+1,ik,it)+aparnew(ig,ik,it) & 
              -apar(ig+1,ik,it)-apar(ig,ik,it)

         source(ig) = anon(ie,is)*(-2.0*vpar(ig,isgn,iglo)*phi_m &
              -zstm(spec,is)*vpac(ig,isgn,iglo) &
              *(aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m &
              -zi*wdrift(ig,iglo)*phi_p) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*delt*wunits(ik)*ufac(ie,is)) &
              *(phi_p - apar_p*stm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik))
      end do
    end subroutine set_source

  end subroutine get_source_term

  subroutine invert_rhs_1 &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, &
        iglo, sourcefac)
    use dist_fn_arrays, only: g, gnew
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, ng2, lmax, forbid
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac

    integer :: ig, ik, it, il, ie, isgn, is
    integer :: ilmin
    complex :: beta1
    complex, dimension (-ntgrid:ntgrid,2) :: source, g1, g2
    logical :: kperiod_flag

    call prof_entering ("invert_rhs_1", "dist_fn")

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    do isgn = 1, 2
       call get_source_term (phi, apar, aperp, phinew, aparnew, aperpnew, &
            istep, isgn, iglo, sourcefac, source(:,isgn))
    end do

    ! gnew is the inhomogeneous solution
    gnew(:,:,iglo) = 0.0

    ! g1 is the homogeneous solution
    g1 = 0.0
    kperiod_flag = boundary_option_switch == boundary_option_self_periodic &
         .or. aky(ik) == 0.0
    if (kperiod_flag) then
       if (il <= ng2+1) then
          g1(-ntgrid,1) = 1.0
          g1( ntgrid,2) = 1.0
       end if
    end if

    ! g2 is the initial condition for the homogeneous solution
    g2 = 0.0
    ! initialize to 1.0 at upper bounce point
    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
       do ig=-ntgrid,ntgrid-1
          if (forbid(ig+1,il).and..not.forbid(ig,il)) g2(ig,2) = 1.0
       end do
    end if

    ! time advance vpar < 0 inhomogeneous part
    do ig = ntgrid-1, -ntgrid, -1
       gnew(ig,2,iglo) &
            = -gnew(ig+1,2,iglo)*r(ig,iglo) + ainv(ig,iglo)*source(ig,2)
    end do

    if (kperiod_flag) then
       ilmin = 1
    else
       ilmin = ng2 + 2
    end if

    ! time advance vpar < 0 homogeneous part
    if (il >= ilmin) then
       do ig = ntgrid-1, -ntgrid, -1
          g1(ig,2) = -g1(ig+1,2)*r(ig,iglo) + g2(ig,2)
       end do
    end if

    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
       ! match boundary conditions at lower bounce point
       ! ???? do some mysterious mucking around with source
       do ig = -ntgrid, ntgrid-1
          if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
             g2(ig,1) = g1(ig+1,2)
             source(ig,1) = gnew(ig+1,2,iglo)
          end if
       end do
    end if

    ! time advance vpar > 0 inhomogeneous part
    if (il <= lmax) then
       do ig = -ntgrid, ntgrid-1
          gnew(ig+1,1,iglo) &
               = -gnew(ig,1,iglo)*r(ig,iglo) + ainv(ig,iglo)*source(ig,1)
       end do
    end if

    ! balancing totally trapped particles
    do ig = -ntgrid, ntgrid
       if (il == ittp(ig)) then
          if (forbid(ig,il)) then
             gnew(ig,1,iglo) = 0.0
          else
             gnew(ig,1,iglo) = gnew(ig,2,iglo)
          end if
       end if
    end do

    ! time advance vpar > 0 homogeneous part
    if (il >= ilmin) then
       do ig = -ntgrid, ntgrid - 1
          g1(ig+1,1) = -g1(ig,1)*r(ig,iglo) + g2(ig,1)
       end do
    end if

    ! add correct amount of homogeneous solution
    if (kperiod_flag .and. il <= ng2+1) then
       beta1 = (gnew(ntgrid,1,iglo) - gnew(-ntgrid,1,iglo)) &
            /(1.0 - g1(ntgrid,1))
       gnew(:,1,iglo) = gnew(:,1,iglo) + beta1*g1(:,1)
       beta1 = (gnew(-ntgrid,2,iglo) - gnew(ntgrid,2,iglo)) &
            /(1.0 - g1(-ntgrid,2))
       gnew(:,2,iglo) = gnew(:,2,iglo) + beta1*g1(:,2)
    end if

    ! add correct amount of homogeneous solution
    if (il >= ng2+2 .and. il <= lmax) then
       beta1 = 0.0
       do ig = ntgrid-1, -ntgrid, -1
          if (ittp(ig) == il) cycle
          if (forbid(ig,il)) then
             beta1 = 0.0
          else if (forbid(ig+1,il)) then
             beta1 = (gnew(ig,1,iglo) - gnew(ig,2,iglo))/(1.0 - g1(ig,1))
          end if
          gnew(ig,:,iglo) = gnew(ig,:,iglo) + beta1*g1(ig,:)
       end do
    end if

    ! zero out spurious gnew outside trapped boundary
    where (forbid(:,il))
       gnew(:,1,iglo) = 0.0
       gnew(:,2,iglo) = 0.0
    end where

    call prof_leaving ("invert_rhs_1", "dist_fn")
  end subroutine invert_rhs_1

  subroutine invert_rhs_untrapped_1 &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, &
        iglo, sourcefac, bc)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac
    complex, intent (in) :: bc

    integer :: ig
    complex, dimension (-ntgrid:ntgrid) :: source, g1
    integer, parameter :: isgn = 1

    call prof_entering ("invert_rhs_untrapped_1", "dist_fn")

    call get_source_term (phi, apar, aperp, phinew, aparnew, aperpnew, &
         istep, isgn, iglo, sourcefac, source)

    gnew(:,isgn,iglo) = 0.0
    gnew(-ntgrid,isgn,iglo) = bc

    do ig = -ntgrid, ntgrid-1
       gnew(ig+1,1,iglo) &
            = -gnew(ig,1,iglo)*r(ig,iglo) + ainv(ig,iglo)*source(ig)
    end do

    call prof_leaving ("invert_rhs_untrapped_1", "dist_fn")

  end subroutine invert_rhs_untrapped_1

  subroutine invert_rhs_untrapped_2 &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, &
        iglo, sourcefac, bc)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac
    complex, intent (in) :: bc

    integer :: ig
    complex, dimension (-ntgrid:ntgrid) :: source, g1
    integer, parameter :: isgn = 2

    call prof_entering ("invert_rhs_untrapped_2", "dist_fn")

    call get_source_term (phi, apar, aperp, phinew, aparnew, aperpnew, &
         istep, isgn, iglo, sourcefac, source)

    gnew(:,isgn,iglo) = 0.0
    gnew(ntgrid,isgn,iglo) = bc

    do ig = ntgrid-1, -ntgrid, -1
       gnew(ig,2,iglo) &
            = -gnew(ig+1,2,iglo)*r(ig,iglo) + ainv(ig,iglo)*source(ig)
    end do

    call prof_leaving ("invert_rhs_untrapped_2", "dist_fn")

  end subroutine invert_rhs_untrapped_2

  subroutine invert_rhs_connected_bc_mp &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, sourcefac)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use mp, only: iproc, nproc, send, receive
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    complex, intent (in) :: sourcefac

    integer :: i, iglo, ipass
    integer :: idp, ipto, ipfrom, iadp
    complex :: bc

    ! do unconnected
    do i = 1, n_unconnected
       iglo = iglo_unconnected(i)
       call invert_rhs_1 (phi, apar, aperp, phinew, aparnew, aperpnew, &
            istep, iglo, sourcefac)
    end do

    ! set bcin to zero
    bcleft_in = 0.0
    bcright_in = 0.0

    do ipass = 1, npass_connected
       do i = 1, connected(ipass)%n1
          iglo = connected(ipass)%iglo1(i)
          bc = bcleft_in(connected(ipass)%ibc1_in(i))
          do
             call invert_rhs_untrapped_1 &
                  (phi, apar, aperp, phinew, aparnew, aperpnew, &
                   istep, iglo, sourcefac, bc)
             bc = gnew(ntgrid,1,iglo)
             if (connections(iglo)%iproc_right == iproc) then
                iglo = connections(iglo)%iglo_right
                cycle
             else if (connections(iglo)%iproc_right >= 0) then
                bcleft_out(connected(ipass)%ibc1_out(i), &
                           connections(iglo)%iproc_right) = bc
                exit
             else
                exit
             end if
          end do
       end do
       do i = 1, connected(ipass)%n2
          iglo = connected(ipass)%iglo2(i)
          bc = bcright_in(connected(ipass)%ibc2_in(i))
          do
             call invert_rhs_untrapped_2 &
                  (phi, apar, aperp, phinew, aparnew, aperpnew, &
                   istep, iglo, sourcefac, bc)
             bc = gnew(-ntgrid,2,iglo)
             if (connections(iglo)%iproc_left == iproc) then
                iglo = connections(iglo)%iglo_left
                cycle
             else if (connections(iglo)%iproc_left >= 0) then
                bcright_out(connected(ipass)%ibc2_out(i), &
                            connections(iglo)%iproc_left) = bc
                exit
             else
                exit
             end if
          end do
       end do

       ! send/receive
       do idp = 1, nproc-1
          ipto = mod(iproc + idp, nproc)
          ipfrom = mod(iproc + nproc - idp, nproc)
          iadp = min(idp, nproc - idp)

          if (mod(iproc/iadp,2) == 0) then
             ! send +
             if (connected(ipass)%nout1(ipto) > 0) then
                call send (bcleft_out(:connected(ipass)%nout1(ipto),ipto), &
                           ipto, idp)
             end if
             ! send -
             if (connected(ipass)%nout2(ipto) > 0) then
                call send (bcright_out(:connected(ipass)%nout2(ipto),ipto), &
                           ipto, idp)
             end if

             ! receive +
             if (connected(ipass)%nin1_start(ipfrom) &
                  <= connected(ipass)%nin1_end(ipfrom)) &
             then
                call receive (bcleft_in(connected(ipass)%nin1_start(ipfrom) &
                                        :connected(ipass)%nin1_end(ipfrom)), &
                                        ipfrom, idp)
             end if
             ! receive -
             if (connected(ipass)%nin2_start(ipfrom) &
                  <= connected(ipass)%nin2_end(ipfrom)) &
             then
                call receive (bcright_in(connected(ipass)%nin2_start(ipfrom) &
                                         :connected(ipass)%nin2_end(ipfrom)), &
                                         ipfrom, idp)
             end if
          else
             ! receive +
             if (connected(ipass)%nin1_start(ipfrom) &
                  <= connected(ipass)%nin1_end(ipfrom)) &
             then
                call receive (bcleft_in(connected(ipass)%nin1_start(ipfrom) &
                                        :connected(ipass)%nin1_end(ipfrom)), &
                                        ipfrom, idp)
             end if
             ! receive -
             if (connected(ipass)%nin2_start(ipfrom) &
                  <= connected(ipass)%nin2_end(ipfrom)) &
             then
                call receive (bcright_in(connected(ipass)%nin2_start(ipfrom) &
                                         :connected(ipass)%nin2_end(ipfrom)), &
                                         ipfrom, idp)
             end if

             ! send +
             if (connected(ipass)%nout1(ipto) > 0) then
                call send (bcleft_out(:connected(ipass)%nout1(ipto),ipto), &
                           ipto, idp)
             end if
             ! send -
             if (connected(ipass)%nout2(ipto) > 0) then
                call send (bcright_out(:connected(ipass)%nout2(ipto),ipto), &
                           ipto, idp)
             end if
          end if
       end do
    end do
  end subroutine invert_rhs_connected_bc_mp

  subroutine invert_rhs_with_guards &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, &
        iglo, sourcefac)
    use dist_fn_arrays, only: g, gnew
    use theta_grid, only: ntgrid, theta
    use le_grids, only: nlambda, ng2, lmax, forbid
    use kt_grids, only: aky, theta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac

    integer :: ig, ik, it, il, ie, isgn, is
    integer :: ilmin
    complex :: beta1
    complex, dimension (-ntgrid:ntgrid,2) :: source, g1, g2
    logical :: kperiod_flag
    integer :: ntgl, ntgr

    call prof_entering ("invert_rhs_1", "dist_fn")

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    do isgn = 1, 2
       call get_source_term (phi, apar, aperp, phinew, aparnew, aperpnew, &
            istep, isgn, iglo, sourcefac, source(:,isgn))
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.false. .and. istep == 1) then!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ig = -ntgrid, ntgrid
       write (9, "(19(x,e12.6))") theta(ig), aky(ik), theta0(ik,it), &
            source(ig,1), source(ig,2), phi(ig,ik,it), phinew(ig,ik,it), &
            theta(ig) - theta0(ik,it)
    end do
    write (9, "()")
    end if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! gnew is the inhomogeneous solution
    gnew(:,:,iglo) = 0.0

    ! g1 is the homogeneous solution
    g1 = 0.0
    kperiod_flag = boundary_option_switch == boundary_option_self_periodic &
         .or. aky(ik) == 0.0

    if (istep == 0 .and. .not. kperiod_flag) then
       ntgl = -ntgrid
       ntgr = ntgrid
    else
       if (itleft(ik,it) < 0 .or. kperiod_flag) then
          ntgl = -ntgrid_noguard
       else
          ntgl = -ntgrid
          !gnew(-ntgrid,1,iglo) = g(-ntgrid,1,iglo)
       end if
       if (itright(ik,it) < 0 .or. kperiod_flag) then
          ntgr = ntgrid_noguard
       else
          ntgr = ntgrid
          !gnew(ntgrid,2,iglo) = g(ntgrid,2,iglo)
       end if
    end if

    if (kperiod_flag) then
       if (il <= ng2+1) then
          g1(ntgl,1) = 1.0
          g1(ntgr,2) = 1.0
       end if
    end if

    ! g2 is the initial condition for the homogeneous solution
    g2 = 0.0
    ! initialize to 1.0 at upper bounce point
    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
       do ig=ntgl,ntgr-1
          if (forbid(ig+1,il).and..not.forbid(ig,il)) g2(ig,2) = 1.0
       end do
    end if

    ! time advance vpar < 0 inhomogeneous part
    do ig = ntgr-1, ntgl, -1
       gnew(ig,2,iglo) &
            = -gnew(ig+1,2,iglo)*r(ig,iglo) + ainv(ig,iglo)*source(ig,2)
    end do

    if (kperiod_flag) then
       ilmin = 1
    else
       ilmin = ng2 + 2
    end if

    ! time advance vpar < 0 homogeneous part
    if (il >= ilmin) then
       do ig = ntgr-1, ntgl, -1
          g1(ig,2) = -g1(ig+1,2)*r(ig,iglo) + g2(ig,2)
       end do
    end if

    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
       ! match boundary conditions at lower bounce point
       ! ???? do some mysterious mucking around with source
       do ig = ntgl, ntgr-1
          if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
             g2(ig,1) = g1(ig+1,2)
             source(ig,1) = gnew(ig+1,2,iglo)
          end if
       end do
    end if

    ! time advance vpar > 0 inhomogeneous part
    if (il <= lmax) then
       do ig = ntgl, ntgr-1
          gnew(ig+1,1,iglo) &
               = -gnew(ig,1,iglo)*r(ig,iglo) + ainv(ig,iglo)*source(ig,1)
       end do
    end if

    ! balancing totally trapped particles
    do ig = ntgl, ntgr
       if (il == ittp(ig)) then
          if (forbid(ig,il)) then
             gnew(ig,1,iglo) = 0.0
          else
             gnew(ig,1,iglo) = gnew(ig,2,iglo)
          end if
       end if
    end do

    ! time advance vpar > 0 homogeneous part
    if (il >= ilmin) then
       do ig = ntgl, ntgr-1
          g1(ig+1,1) = -g1(ig,1)*r(ig,iglo) + g2(ig,1)
       end do
    end if

    ! add correct amount of homogeneous solution
    if (kperiod_flag .and. il <= ng2+1) then
       beta1 = (gnew(ntgr,1,iglo) - gnew(ntgl,1,iglo)) &
            /(1.0 - g1(ntgr,1))
       gnew(:,1,iglo) = gnew(:,1,iglo) + beta1*g1(:,1)
       beta1 = (gnew(ntgl,2,iglo) - gnew(ntgr,2,iglo)) &
            /(1.0 - g1(ntgl,2))
       gnew(:,2,iglo) = gnew(:,2,iglo) + beta1*g1(:,2)
    end if

    ! add correct amount of homogeneous solution
    if (il >= ng2+2 .and. il <= lmax) then
       beta1 = 0.0
       do ig = ntgr-1, ntgl, -1
          if (ittp(ig) == il) cycle
          if (forbid(ig,il)) then
             beta1 = 0.0
          else if (forbid(ig+1,il)) then
             beta1 = (gnew(ig,1,iglo) - gnew(ig,2,iglo))/(1.0 - g1(ig,1))
          end if
          gnew(ig,:,iglo) = gnew(ig,:,iglo) + beta1*g1(ig,:)
       end do
    end if

    ! zero out spurious gnew outside trapped boundary
    where (forbid(:,il))
       gnew(:,1,iglo) = 0.0
       gnew(:,2,iglo) = 0.0
    end where

    call prof_leaving ("invert_rhs_1", "dist_fn")
  end subroutine invert_rhs_with_guards

  subroutine invert_rhs_conn_via_guards &
       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, sourcefac)
    use dist_fn_arrays, only: g
    use theta_grid, only: ntgrid, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use gs2_layouts, only: g_lo
    use mp, only: iproc, nproc, send, receive
    use redistribute, only: fill
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep
    complex, intent (in) :: sourcefac

    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) &
         :: phig, aparg, aperpg, phinewg, aparnewg, aperpnewg
    integer :: ik, it, itl, itr
    integer :: iglo

    phig = phi
    aparg = apar
    aperpg = aperp
    phinewg = phinew
    aparnewg = aparnew
    aperpnewg = aperpnew
    if (istep > 0) then
       do it = 1, ntheta0
          do ik = 1, naky
             itl = itleft(ik,it)
             itr = itright(ik,it)

             if (itl >= 0) then
                phig     (ldlo:ldhi,ik,it) = phi     (rslo:rshi,ik,itl)
                aparg    (ldlo:ldhi,ik,it) = apar    (rslo:rshi,ik,itl)
                aperpg   (ldlo:ldhi,ik,it) = aperp   (rslo:rshi,ik,itl)
                phinewg  (ldlo:ldhi,ik,it) = phinew  (rslo:rshi,ik,itl)
                aparnewg (ldlo:ldhi,ik,it) = aparnew (rslo:rshi,ik,itl)
                aperpnewg(ldlo:ldhi,ik,it) = aperpnew(rslo:rshi,ik,itl)
             else
                phig     (ldlo:ldhi,ik,it) = 0.0
                aparg    (ldlo:ldhi,ik,it) = 0.0
                aperpg   (ldlo:ldhi,ik,it) = 0.0
                phinewg  (ldlo:ldhi,ik,it) = 0.0
                aparnewg (ldlo:ldhi,ik,it) = 0.0
                aperpnewg(ldlo:ldhi,ik,it) = 0.0
             end if

             if (itr >= 0) then
                phig     (rdlo:rdhi,ik,it) = phi     (lslo:lshi,ik,itr)
                aparg    (rdlo:rdhi,ik,it) = apar    (lslo:lshi,ik,itr)
                aperpg   (rdlo:rdhi,ik,it) = aperp   (lslo:lshi,ik,itr)
                phinewg  (rdlo:rdhi,ik,it) = phinew  (lslo:lshi,ik,itr)
                aparnewg (rdlo:rdhi,ik,it) = aparnew (lslo:lshi,ik,itr)
                aperpnewg(rdlo:rdhi,ik,it) = aperpnew(lslo:lshi,ik,itr)
             else
                phig     (rdlo:rdhi,ik,it) = 0.0
                aparg    (rdlo:rdhi,ik,it) = 0.0
                aperpg   (rdlo:rdhi,ik,it) = 0.0
                phinewg  (rdlo:rdhi,ik,it) = 0.0
                aparnewg (rdlo:rdhi,ik,it) = 0.0
                aperpnewg(rdlo:rdhi,ik,it) = 0.0
             end if

          end do
       end do

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          if (connections(iglo)%iproc_left < 0) then
             g(ldlo:ldhi,1,iglo) = 0.0
             g(ldlo:ldhi,2,iglo) = 0.0
          end if
          if (connections(iglo)%iproc_right < 0) then
             g(rdlo:rdhi,1,iglo) = 0.0
             g(rdlo:rdhi,2,iglo) = 0.0
          end if
       end do

       call fill (gc_from_left, g, g)
       call fill (gc_from_right, g, g)

    else
       ! no guard copies when initializing response matrix
       ! (plus, making the copies would make the matrix singular)
       ! g = 0 when initializing response matrix
    end if

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       call invert_rhs_with_guards &
            (phig, aparg, aperpg, phinewg, aparnewg, aperpnewg, &
            istep, iglo, sourcefac)
    end do
  end subroutine invert_rhs_conn_via_guards

  subroutine invert_rhs (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, ng2
    use run_parameters, only: delt
    use gs2_layouts, only: g_lo, il_idx
    use gs2_time, only: stime
    use constants
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep

    integer :: iglo, il

    real :: time
    complex :: sourcefac

    call prof_entering ("invert_rhs", "dist_fn")

    time = stime()
    if (time > t0) then
       sourcefac = exp(-zi*omega0*time+gamma0*time)
    else
       sourcefac = (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time+gamma0*time)
    end if

    select case (boundary_option_switch)
    case (boundary_option_connected)
!       if (use_shmem_in_connected) then
!          call invert_rhs_connected_bc_shmem &
!               (phi, apar, aperp, phinew, aparnew, aperpnew, istep, sourcefac)
!       else
          call invert_rhs_connected_bc_mp &
               (phi, apar, aperp, phinew, aparnew, aperpnew, istep, sourcefac)
!       end if
    case (boundary_option_conn_via_guards)
       call invert_rhs_conn_via_guards &
            (phi, apar, aperp, phinew, aparnew, aperpnew, istep, sourcefac)
    case (boundary_option_alternate_zero)
       if (nlambda > ng2) then
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             il = il_idx(g_lo,iglo)
             if (il >= ng2+2) then
                call invert_rhs_1 &
                     (phi, apar, aperp, phinew, aparnew, aperpnew, &
                      istep, iglo, sourcefac)
             else
                call invert_rhs_untrapped_1 &
                     (phi, apar, aperp, phinew, aparnew, aperpnew, &
                      istep, iglo, sourcefac, cmplx(0.0,0.0))
                call invert_rhs_untrapped_2 &
                     (phi, apar, aperp, phinew, aparnew, aperpnew, &
                      istep, iglo, sourcefac, cmplx(0.0,0.0))
             end if
          end do
       else
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             call invert_rhs_untrapped_1 &
                  (phi, apar, aperp, phinew, aparnew, aperpnew, &
                   istep, iglo, sourcefac, cmplx(0.0,0.0))
             call invert_rhs_untrapped_2 &
                  (phi, apar, aperp, phinew, aparnew, aperpnew, &
                   istep, iglo, sourcefac, cmplx(0.0,0.0))
          end do
       end if
    case default
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          call invert_rhs_1 (phi, apar, aperp, phinew, aparnew, aperpnew, &
               istep, iglo, sourcefac)
       end do
    end select

    call prof_leaving ("invert_rhs", "dist_fn")
  end subroutine invert_rhs

  subroutine getan (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2, aj0, aj1, gnew
    use dist_fn_arrays, only: aj0f, aj1f
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species
    use run_parameters, only: beta
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt

    integer :: ig, ik, it, isgn

    call prof_entering ("getan", "dist_fn")
    do isgn = 1, 2
       g0(:,isgn,:) = aj0f*gnew(:,isgn,:)
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, antot)

    do isgn = 1, 2
       g0(:,isgn,:) = aj0f*vpa(:,isgn,:)*gnew(:,isgn,:)
    end do
    wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
    call integrate_species (g0, wgt, antota)

    do isgn = 1, 2
       g0(:,isgn,:) = aj1f*vperp2*gnew(:,isgn,:)
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g0, wgt, antotp)
    call prof_leaving ("getan", "dist_fn")
  end subroutine getan

  subroutine init_fieldeq
    use dist_fn_arrays, only: aj0, aj1, vperp2, kperp2
    use species, only: nspec, spec, stm, has_electron_species
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: anon, integrate_species
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use run_parameters, only: tite
    implicit none
    integer :: iglo, isgn
    integer :: ik, it, ie, is
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: tot
    real, dimension (nspec) :: wgt
    logical :: done = .false.

    if (done) return
    done = .true.

    allocate (gridfac1(-ntgrid:ntgrid,naky,ntheta0))
    gridfac1 = 1.0
    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       ! nothing
    case (boundary_option_connected)
       do it = 1, ntheta0
          do ik = 1, naky
             if (aky(ik) == 0.0) cycle
             if (itleft(ik,it) < 0) gridfac1(-ntgrid,ik,it) = gridfac
             if (itright(ik,it) < 0) gridfac1(ntgrid,ik,it) = gridfac
          end do
       end do
    case (boundary_option_conn_via_guards)
       do it = 1, ntheta0
          do ik = 1, naky
             if (aky(ik) == 0.0) cycle
             if (itleft(ik,it) < 0) gridfac1(-ntgrid_noguard,ik,it) = gridfac
             if (itright(ik,it) < 0) gridfac1(ntgrid_noguard,ik,it) = gridfac
          end do
       end do
    case default
       do ik = 1, naky
          if (aky(ik) == 0.0) cycle
          gridfac1(-ntgrid,ik,:) = gridfac
          gridfac1(ntgrid,ik,:) = gridfac
       end do
    end select

    allocate (gamtot(-ntgrid:ntgrid,naky,ntheta0))
    allocate (gamtot1(-ntgrid:ntgrid,naky,ntheta0))
    allocate (gamtot2(-ntgrid:ntgrid,naky,ntheta0))
    if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       allocate (gamtot3(-ntgrid:ntgrid,naky,ntheta0))
    endif
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = (1.0 - aj0(:,iglo)**2)*anon(ie,is)
       end do
    end do
    wgt = spec%z*spec%z*spec%dens/spec%temp
    call integrate_species (g0, wgt, tot)
    gamtot = real(tot) + kperp2*poisfac
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(:,iglo)*aj1(:,iglo) &
               *2.0*vperp2(:,iglo)*anon(ie,is)
       end do
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj1(:,iglo)**2*2.0*vperp2(:,iglo)**2*anon(ie,is)
       end do
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g0, wgt, tot)
    gamtot2 = real(tot)

! adiabatic electrons 
    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_yavg) then
          do ik = 1, naky
             if (aky(ik) > epsilon(0.0)) gamtot(:,ik,:) = gamtot(:,ik,:) + tite
          end do
       elseif (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          gamtot  = gamtot + tite
          gamtot3 = (gamtot-tite) / gamtot
       else
          gamtot = gamtot + tite ! fieldlineavg broken right now !!??
       endif
    endif
  end subroutine init_fieldeq

  subroutine getfieldeq1 (phi, apar, aperp, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use dist_fn_arrays, only: kperp2
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky, theta0
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    real, allocatable, dimension(:,:), save :: fl_avg, awgt
    integer :: ik, it
    logical :: first = .true.
    
    if (first) allocate (fl_avg(naky, ntheta0))
    fl_avg = 0.

    if (.not. has_electron_species(spec) &
         .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then

       if (first) then 
          allocate (awgt(naky, ntheta0))
          awgt = 0.
          do it = 1, ntheta0
             do ik = 1, naky
                if (aky(ik) > epsilon(0.0)) cycle
                awgt(ik,it) = 1.0/sum(delthet*jacob*gamtot3(:,ik,it))
             end do
          end do
       endif

       do it = 1, ntheta0
          do ik = 1, naky
             fl_avg(ik,it) = tite*sum(delthet*jacob*antot(:,ik,it))*awgt(ik,it)
          end do
       end do
    end if
    
    do it = 1, ntheta0
       do ik = 1, naky
          fieldeq(:,ik,it) &
               = antot(:,ik,it) + aperp(:,ik,it)*gamtot1(:,ik,it) &
                 - gamtot(:,ik,it)*gridfac1(:,ik,it)*phi(:,ik,it) &
                 + fl_avg(ik,it)
          fieldeqa(:,ik,it) &
               = antota(:,ik,it) &
                 - kperp2(:,ik,it)*gridfac1(:,ik,it)*apar(:,ik,it)
          fieldeqp(:,ik,it) &
               = (antotp(:,ik,it) + aperp(:,ik,it)*gamtot2(:,ik,it) &
                  + 0.5*phi(:,ik,it)*gamtot1(:,ik,it))*beta*apfac/bmag**2 &
                 + aperp(:,ik,it)*gridfac1(:,ik,it)
       end do
    end do

    first = .false.

  end subroutine getfieldeq1

  subroutine getfieldeq (phi, apar, aperp, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: antot, antota, antotp

    call getan (antot, antota, antotp)
    call getfieldeq1 (phi, apar, aperp, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)
  end subroutine getfieldeq

  subroutine getfieldeq0 (phi, apar, aperp, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: antot, antota, antotp

    integer :: ik, it, itl, itr

    call getan (antot, antota, antotp)
    if (.false.) then
    if (boundary_option_switch == boundary_option_conn_via_guards) then
       do it = 1, ntheta0
          do ik = 1, naky
             itl = itleft(ik,it)
             itr = itright(ik,it)
             if (itl >= 0) then
                antot (ldlo:ldhi,ik,it) = antot (rslo:rshi,ik,itl)
                antota(ldlo:ldhi,ik,it) = antota(rslo:rshi,ik,itl)
                antotp(ldlo:ldhi,ik,it) = antotp(rslo:rshi,ik,itl)
             else
                antot (ldlo:ldhi,ik,it) = 0.0
                antota(ldlo:ldhi,ik,it) = 0.0
                antotp(ldlo:ldhi,ik,it) = 0.0
             end if
             if (itr >= 0) then
                antot (rdlo:rdhi,ik,it) = antot (lslo:lshi,ik,itr)
                antota(rdlo:rdhi,ik,it) = antota(lslo:lshi,ik,itr)
                antotp(rdlo:rdhi,ik,it) = antotp(lslo:lshi,ik,itr)
             else
                antot (rdlo:rdhi,ik,it) = 0.0
                antota(rdlo:rdhi,ik,it) = 0.0
                antotp(rdlo:rdhi,ik,it) = 0.0
             end if
          end do
       end do
    end if
    end if
    call getfieldeq1 (phi, apar, aperp, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)
  end subroutine getfieldeq0

  pure function j0 (x)
! A&S, p. 369, 9.4
    implicit none
    real, intent (in) :: x
    real :: j0
    real, parameter, dimension (7) :: a = &
         (/ 1.0000000, -2.2499997, 1.2656208, -0.3163866, &
            0.0444479, -0.0039444, 0.0002100 /)
    real, parameter, dimension (7) :: b = &
         (/  0.79788456, -0.00000770, -0.00552740, -0.00009512, &
             0.00137237, -0.00072805,  0.00014476 /)
    real, parameter, dimension (7) :: c = &
         (/ -0.78539816, -0.04166397, -0.00003954,  0.00262573, &
            -0.00054125, -0.00029333,  0.00013558 /)
    real :: y

    if (x <= 3.0) then
       y = (x/3.0)**2
       j0 = a(1)+y*(a(2)+y*(a(3)+y*(a(4)+y*(a(5)+y*(a(6)+y*a(7))))))
    else
       y = 3.0/x
       j0 = (b(1)+y*(b(2)+y*(b(3)+y*(b(4)+y*(b(5)+y*(b(6)+y*b(7))))))) &
            *cos(x+c(1)+y*(c(2)+y*(c(3)+y*(c(4)+y*(c(5)+y*(c(6)+y*c(7))))))) &
            /sqrt(x)
    end if
  end function j0

  pure function j1 (x)
! A&S, p. 370, 9.4 j1 = 1/x J_1(x)
    implicit none
    real, intent (in) :: x
    real :: j1
    real, parameter, dimension (7) :: a = &
         (/  0.50000000, -0.56249985,  0.21093573, -0.03954289, &
             0.00443319, -0.00031761,  0.00001109 /)
    real, parameter, dimension (7) :: b = &
         (/  0.79788456,  0.00000156,  0.01659667,  0.00017105, &
            -0.00249511,  0.00113653,  0.00020033 /)
    real, parameter, dimension (7) :: c = &
         (/ -2.35619449,  0.12499612,  0.00005650,  -0.00637879, &
             0.00074348,  0.00079824, -0.00029166 /)
    real :: y

    if (x <= 3.0) then
       y = (x/3.0)**2
       j1 = a(1)+y*(a(2)+y*(a(3)+y*(a(4)+y*(a(5)+y*(a(6)+y*a(7))))))
    else
       y = 3.0/x
       j1 = (b(1)+y*(b(2)+y*(b(3)+y*(b(4)+y*(b(5)+y*(b(6)+y*b(7))))))) &
            *cos(x+c(1)+y*(c(2)+y*(c(3)+y*(c(4)+y*(c(5)+y*(c(6)+y*c(7))))))) &
            /x**1.5
    end if
  end function j1

  subroutine flux (phi, apar, aperp, &
       pflux, qheat, vflux, pmflux, qmheat, vmflux, anorm)
    use species, only: nspec, spec, stm
    use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, e
    use dist_fn_arrays, only: g, aj0, vpac, vpa
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use mp, only: proc0
    use run_parameters, only: woutunits
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    real, dimension (:,:,:), intent (out) :: pflux, qheat, vflux
    real, dimension (:,:,:), intent (out) :: pmflux, qmheat, vmflux
    real, intent (out) :: anorm
    real, dimension (-ntgrid:ntgrid,naky) :: dnorm
    integer :: ig, ik, is, isgn
    integer :: iglo

    integer :: ng

    ng = ntgrid_noguard

    if (proc0) then
       pflux = 0.0;   qheat = 0.0;   vflux = 0.0
       pmflux = 0.0;  qmheat = 0.0;  vmflux = 0.0
    end if

    do ik = 1, naky
       do ig = -ng, ng
          dnorm(ig,ik) = delthet(ig)*grho(ig)/gradpar(ig)/bmag(ig) &
                                    *woutunits(ik)/sqrt(2.0)
       end do
    end do

    anorm = sum(spread(dnorm(-ng:ng,:),3,ntheta0) &
         *abs(phi(-ng:ng,:,:))**2 + abs(apar(-ng:ng,:,:))**2)

    do isgn = 1, 2
       g0(:,isgn,:) = g(:,isgn,:)*aj0
    end do
    call get_flux (phi, pflux, anorm, dnorm)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
    end do
    call get_flux (phi, qheat, anorm, dnorm)

    do isgn = 1, 2
       g0(:,isgn,:) = g(:,isgn,:)*aj0*vpac(:,isgn,:)
    end do
    call get_flux (phi, vflux, anorm, dnorm)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) &
               = g(:,isgn,iglo)*aj0(:,iglo)*stm(spec,is)*vpa(:,isgn,iglo)
       end do
    end do
    call get_flux (apar, pmflux, anorm, dnorm)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
    end do
    call get_flux (phi, qmheat, anorm, dnorm)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) &
               = g(:,isgn,iglo)*aj0(:,iglo)*stm(spec,is) &
                 *vpa(:,isgn,iglo)*vpac(:,isgn,iglo)
       end do
    end do
    call get_flux (apar, vmflux, anorm, dnorm)
  end subroutine flux

  subroutine get_flux (fld, flx, anorm, dnorm)
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, aky
    use le_grids, only: integrate
    use gs2_layouts, only: geint_lo, ik_idx, it_idx, is_idx, proc_id, idx_local
    use mp, only: proc0, send, receive, barrier
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (out) :: flx
    real :: anorm
    real, dimension (-ntgrid:,:) :: dnorm
    complex, &
         dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc) :: &
         geint
    real :: tmp
    integer :: igeint
    integer :: ik, it, is

    integer :: ng

    ng = ntgrid_noguard

    if (abs(anorm) < 2.0*epsilon(0.0)) then
       if (proc0) flx = 0.0
       return
    end if
    call integrate (g0, geint)
    do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       tmp = sum(-aimag(geint(-ng:ng,igeint)*conjg(fld(-ng:ng,ik,it))) &
                  *dnorm(-ng:ng,ik)*aky(ik))/anorm
       geint(0,igeint) = tmp
    end do
    do igeint = geint_lo%llim_world, geint_lo%ulim_world
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       call barrier
       if (proc0) then
          if (idx_local(geint_lo,igeint)) then
             flx(ik,it,is) = real(geint(0,igeint))
          else
             call receive (flx(ik,it,is), proc_id(geint_lo,igeint), igeint)
          end if
       else if (idx_local(geint_lo,igeint)) then
          call send (real(geint(0,igeint)), 0, igeint)
       end if
    end do
  end subroutine get_flux

  subroutine neoclassical_flux (pflux, qflux)
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, gradpar, bmag, delthet
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid, e, integrate
    use dist_fn_arrays, only: g
    use gs2_layouts, only: g_lo, geint_lo, ik_idx, it_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id, idx_local
    use mp, only: proc0, send, receive, broadcast
    use run_parameters, only: woutunits, delt
    use constants
    implicit none
    complex, dimension (:,:,:), intent (out) :: pflux, qflux
    real, dimension (-ntgrid:ntgrid,naky) :: dnorm
    complex, dimension (naky,ntheta0,nspec) :: anorm
    complex, &
         dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc) :: &
         geint
    integer :: ig, ik, it, ie, is
    integer :: iglo, igeint
    complex :: x

    integer :: ng

    ng = ntgrid_noguard

    if (proc0) then
       pflux = 0.0
       qflux = 0.0
    end if

    do ik = 1, naky
       do ig = -ntgrid, ntgrid
          dnorm(ig,ik) = 1.0/gradpar(ig)/bmag(ig)*woutunits(ik)/sqrt(2.0)
       end do
    end do
    dnorm(-ntgrid,:) = 0.5*dnorm(-ntgrid,:)
    dnorm(ntgrid,:)  = 0.5*dnorm(ntgrid,:)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g0(:,1,iglo) = dnorm(:,ik)
       g0(:,2,iglo) = dnorm(:,ik)
    end do
    call integrate (g0, geint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             igeint = idx(geint_lo,ik,it,is)
             if (idx_local(geint_lo,igeint)) then
                anorm(ik,it,is) &
                     = sum(geint(-ng:ng-1,igeint)*delthet(-ng:ng-1))
             end if
             call broadcast (anorm(ik,it,is), proc_id(geint_lo,igeint))
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       g0(:,1,iglo) = dnorm(:,ik)*zi*wdrift(:,iglo)/delt*g(:,1,iglo)
       g0(:,2,iglo) = dnorm(:,ik)*zi*wdrift(:,iglo)/delt*g(:,2,iglo)
    end do
    call integrate (g0, geint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             igeint = idx(geint_lo,ik,it,is)
             if (proc0) then
                if (idx_local(geint_lo,igeint)) then
                   pflux(ik,it,is) &
                        = sum(geint(-ng:ng-1,igeint) &
                              *delthet(-ng:ng-1))/anorm(ik,it,is)
                else
                   call receive (pflux(ik,it,is), proc_id(geint_lo,igeint))
                end if
             else if (idx_local(geint_lo,igeint)) then
                x = sum(geint(-ng:ng-1,igeint)*delthet(-ng:ng-1)) &
                     /anorm(ik,it,is)
                call send (x, 0)
             end if
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g0(:,:,iglo) = g0(:,:,iglo)*e(ie,is)**(-1.5)
    end do
    call integrate (g0, geint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             igeint = idx(geint_lo,ik,it,is)
             if (proc0) then
                if (idx_local(geint_lo,igeint)) then
                   qflux(ik,it,is) &
                        = sum(geint(-ng:ng-1,igeint) &
                              *delthet(-ng:ng-1))/anorm(ik,it,is)
                else
                   call receive (qflux(ik,it,is), proc_id(geint_lo,igeint))
                end if
             else if (idx_local(geint_lo,igeint)) then
                x = sum(geint(-ng:ng-1,igeint)*delthet(-ng:ng-1)) &
                     /anorm(ik,it,is)
                call send (x, 0)
             end if
          end do
       end do
    end do
  end subroutine neoclassical_flux

  subroutine init_fieldcheck
  end subroutine init_fieldcheck

  subroutine finish_fieldcheck
  end subroutine finish_fieldcheck

  subroutine fieldcheck (phi, apar, aperp)
    use file_utils, only: open_output_file, close_output_file
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, aky, ntheta0, theta0
    use dist_fn_arrays, only: kperp2
    use mp, only: proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0)::fieldeq,fieldeqa,fieldeqp
    integer :: ig, ik, it, unit

    call getfieldeq (phi, apar, aperp, fieldeq, fieldeqa, fieldeqp)

    if (proc0) then
       call open_output_file (unit, ".fieldcheck")
       do it = 1, ntheta0
          do ik = 1, naky
             do ig = -ntgrid, ntgrid
                write (unit,*) "j,theta,ik,aky,it,theta0 ", &
                     ig, theta(ig), ik, aky(ik), it, theta0(ik,it)
                write (unit,"(8(1pe12.5,x))") phi(ig,ik,it),fieldeq(ig,ik,it),&
                     gamtot(ig,ik,it)
                write (unit,"(8(1pe12.5,x))") &
                     kperp2(ig,ik,it)*gridfac1(ig,ik,it)*apar(ig,ik,it), &
                     fieldeqa(ig,ik,it), gamtot1(ig,ik,it)
                write (unit,"(8(1pe12.5,x))") &
                     gridfac1(ig,ik,it)*aperp(ig,ik,it), fieldeqp(ig,ik,it), &
                     gamtot2(ig,ik,it)
             end do
          end do
       end do
       call close_output_file (unit)
     end if
  end subroutine fieldcheck

  subroutine init_vortcheck
  end subroutine init_vortcheck

  subroutine finish_vortcheck
  end subroutine finish_vortcheck

  subroutine vortcheck (phi, apar, aperp)
    use file_utils, only: open_output_file, close_output_file
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use le_grids, only: e, al, anon, integrate_species
    use kt_grids, only: naky, aky, ntheta0, theta0
    use dist_fn_arrays, only: aj0, aj1, vperp2, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use run_parameters, only: delt, fphi, faperp
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    real, dimension (nspec) :: wgt
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0,2) :: apchk
    integer :: iglo, ig, ik, it, il, ie, is
    real :: temp, z
    real, dimension (-ntgrid:ntgrid) :: vplus, cvtot, gbtot, delwd
    integer :: unit
    
    apchk = 0.0

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       temp = spec(is)%temp
       z = spec(is)%z
       vplus = e(ie,is)*al(il)*bmag
       cvtot = cvdrift + theta0(ik,it)*cvdrift0
       gbtot = gbdrift + theta0(ik,it)*gbdrift0
       delwd = temp/z*delt*(gbtot-cvtot)*vplus/2.0

       g0(:,1,iglo) &
         = (gnew(:,1,iglo) &
            +anon(ie,is)*2.0*vperp2(:,iglo)*aj1(:,iglo)*faperp*aperp(:,ik,it) &
            +fphi*z*anon(ie,is)*phi(:,ik,it)*aj0(:,iglo)/temp) &
          *delwd
    end do
    wgt = spec%dens*spec%z
    call integrate_species (g0, wgt, apchk(:,:,:,1))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       temp = spec(is)%temp
       z = spec(is)%z

       g0(:,1,iglo) = aperp(:,ik,it)*faperp*vperp2(:,iglo)*temp/z
    end do
    wgt = spec%dens*spec%z
    call integrate_species (g0, wgt, apchk(:,:,:,2))

    call open_output_file (unit, ".vortcheck")
    do it = 1, ntheta0
       do ik = 1, naky
          write (unit,*) 'aky=',aky(ik), ' theta0=',theta0(ik,it)
          do ig = -ntgrid, ntgrid
             write (unit,*) apchk(ig,ik,it,1), apchk(ig,ik,it,2)
          end do
       end do
    end do
    call close_output_file (unit)
  end subroutine vortcheck

  subroutine init_intcheck
  end subroutine init_intcheck

  subroutine finish_intcheck
  end subroutine finish_intcheck

  subroutine intcheck
    use le_grids, only: le_intcheck => intcheck
    implicit none
    call le_intcheck (g0)
  end subroutine intcheck

  subroutine reset_init

    initializing  = .true.
    initialized = .false.

  end subroutine reset_init

end module dist_fn


!  subroutine invert_rhs_connected_bc_shmem &
!       (phi, apar, aperp, phinew, aparnew, aperpnew, istep, sourcefac)
!    use dist_fn_arrays, only: gnew
!    use theta_grid, only: ntgrid
!    use gs2_layouts, only: g_lo
!    use mp, only: iproc, barrier
!    use shmem
!    implicit none
!    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
!    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
!    integer, intent (in) :: istep
!    complex, intent (in) :: sourcefac
!
!    complex, dimension (shmem_buff_size,2) :: shmem_buff
!    logical, dimension (shmem_buff_size,2) :: shmem_buff_ready
!    !DIR$ SYMMETRIC shmem_buff, shmem_buff_ready
!
!    logical, dimension (g_lo%llim_proc:g_lo%ulim_alloc,2) :: done
!    integer :: iglo
!    integer :: llim, ulim, llim_new, ulim_new
!    integer :: ipto, ibuffto, ibufffrom
!
!    call barrier
!    shmem_buff_ready = .false.
!    call barrier
!
!    done = .false.
!    llim_new = g_lo%ulim_proc + 1
!    ulim_new = g_lo%llim_proc - 1
!    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!       if (connections(iglo)%iglo_left < 0 &
!            .and. connections(iglo)%iglo_right < 0) &
!       then
!          call invert_rhs_1 (phi, apar, aperp, phinew, aparnew, aperpnew, &
!               istep, iglo, sourcefac)
!          done(iglo,:) = .true.
!          llim_new = min(llim_new,iglo+1)
!          ulim_new = max(ulim_new,iglo-1)
!       else if (connections(iglo)%iglo_left < 0) then
!          call do_invert_1_1 (iglo, cmplx(0.0,0.0))
!       else if (connections(iglo)%iglo_right < 0) then
!          call do_invert_1_2 (iglo, cmplx(0.0,0.0))
!       end if
!    end do
!   
!    do
!       llim = llim_new
!       ulim = ulim_new
!       if (llim > ulim) exit
!       llim_new = ulim + 1
!       ulim_new = llim - 1
!
!       do iglo = llim, ulim
!          if (done(iglo,1) .and. done(iglo,2)) cycle
!          ibufffrom = 1 + iglo - iglo_llim(iproc)
!          if (.not. done(iglo,1) .and. shmem_buff_ready(ibufffrom,1)) then
!             call do_invert_1_1 (iglo, shmem_buff(ibufffrom,1))
!          end if
!          if (.not. done(iglo,2) .and. shmem_buff_ready(ibufffrom,2)) then
!             call do_invert_1_2 (iglo, shmem_buff(ibufffrom,2))
!          end if
!          if (done(iglo,1) .and. done(iglo,2)) then
!             llim_new = min(llim_new,iglo+1)
!             ulim_new = max(ulim_new,iglo-1)
!          end if
!       end do
!
!       llim = llim_new
!       ulim = ulim_new
!       if (llim > ulim) exit
!       llim_new = ulim + 1
!       ulim_new = llim - 1
!
!       do iglo = ulim, llim, -1
!          if (done(iglo,1) .and. done(iglo,2)) cycle
!          ibufffrom = 1 + iglo - iglo_llim(iproc)
!          if (.not. done(iglo,1) .and. shmem_buff_ready(ibufffrom,1)) then
!             call do_invert_1_1 (iglo, shmem_buff(ibufffrom,1))
!          end if
!          if (.not. done(iglo,2) .and. shmem_buff_ready(ibufffrom,2)) then
!             call do_invert_1_2 (iglo, shmem_buff(ibufffrom,2))
!          end if
!          if (done(iglo,1) .and. done(iglo,2)) then
!             llim_new = min(llim_new,iglo+1)
!             ulim_new = max(ulim_new,iglo-1)
!          end if
!       end do
!    end do
!
!  contains
!    subroutine do_invert_1_1 (iglo, bc_in)
!      implicit none
!      integer, intent (in) :: iglo
!      complex, intent (in) :: bc_in
!
!      integer :: ipto, ibuffto
!
!      call invert_rhs_untrapped_1 &
!           (phi, apar, aperp, phinew, aparnew, aperpnew, istep, iglo, &
!            sourcefac, bc_in)
!      done(iglo,1) = .true.
!      ipto = connections(iglo)%iproc_right
!      if (ipto == iproc) then
!         ibuffto = 1 + connections(iglo)%iglo_right - iglo_llim(ipto)
!         shmem_buff(ibuffto,1) = gnew(ntgrid,1,iglo)
!         shmem_buff_ready(ibuffto,1) = .true.
!      else
!         ibuffto = 1 + connections(iglo)%iglo_right - iglo_llim(ipto)
!         call shmem_complex_put &
!              (shmem_buff(ibuffto,1), gnew(ntgrid,1,iglo), 1, ipto)
!         call shmem_logical_put &
!              (shmem_buff_ready(ibuffto,1), .true., 1, ipto)
!      end if
!    end subroutine do_invert_1_1
!
!    subroutine do_invert_1_2 (iglo, bc_in)
!      implicit none
!      integer, intent (in) :: iglo
!      complex, intent (in) :: bc_in
!
!      integer :: ipto, ibuffto
!
!      call invert_rhs_untrapped_2 &
!           (phi, apar, aperp, phinew, aparnew, aperpnew, istep, iglo, &
!            sourcefac, bc_in)
!      done(iglo,2) = .true.
!      ipto = connections(iglo)%iproc_left
!      if (ipto == iproc) then
!         ibuffto = 1 + connections(iglo)%iglo_left - iglo_llim(ipto)
!         shmem_buff(ibuffto,2) = gnew(-ntgrid,2,iglo)
!         shmem_buff_ready(ibuffto,2) = .true.
!      else
!         ibuffto = 1 + connections(iglo)%iglo_left - iglo_llim(ipto)
!         call shmem_complex_put &
!              (shmem_buff(ibuffto,2), gnew(-ntgrid,2,iglo), 1, ipto)
!         call shmem_logical_put &
!              (shmem_buff_ready(ibuffto,2), .true., 1, ipto)
!      end if
!    end subroutine do_invert_1_2
!  end subroutine invert_rhs_connected_bc_shmem
!
