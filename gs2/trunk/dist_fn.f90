module dist_fn
  use init_g, only: ginit
  implicit none

  public :: init_dist_fn
  public :: timeadv
  public :: getfieldeq
  public :: flux, neoclassical_flux
  public :: ginit
  public :: init_intcheck, init_vortcheck, init_fieldcheck
  public :: intcheck, vortcheck, fieldcheck
  public :: finish_intcheck, finish_vortcheck, finish_fieldcheck
  public :: t0, omega0, gamma0, gamma_damp, gamma_damp_e, thetas

  private

  ! knobs
  integer :: kperiod
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  real :: gridfac, apfac, driftknob, poisfac
  real :: t0, omega0, gamma0, gamma_damp, gamma_damp_e, thetas
  real :: phi_ext
  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_zero = 2, source_option_sine = 3, &
       source_option_test1 = 4, source_option_phiext_full = 5, &
       source_option_test2_full = 6

  ! internal arrays

  integer, dimension (:), allocatable :: ittp
  ! (-ntgrid:ntgrid)

  real, dimension (:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:,:,:), allocatable :: wdriftttp
  ! (-ntgrid:ntgrid,naky,ntheta0,negrid,nspec) replicated

  real, dimension (:,:,:), allocatable :: wstar
  ! (naky,negrid,nspec) replicated

  complex, dimension (:,:), allocatable :: a, b, r, ainv
  ! (-ntgrid:ntgrid, -g-layout-)

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2
  ! (-ntgrid:ntgrid,naky,ntheta0) replicated

  real, dimension (:,:), allocatable :: gridfac1
  ! (-ntgrid:ntgrid,naky)

  complex, dimension (:,:,:), allocatable :: g0
  ! (-ntgrid:ntgrid,2, -g-layout-)

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
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_species
    call init_theta_grid
    call init_kt_grids
    call init_le_grids
    call init_run_parameters
    call init_collisions
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

    if (proc0) call read_parameters
    call broadcast_parameters
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
    use text_options
    implicit none
    type (text_option), dimension (7), parameter :: sourceopts = &
         (/ text_option('default', source_option_full), &
            text_option('full', source_option_full), &
            text_option('zero', source_option_zero), &
            text_option('sine', source_option_sine), &
            text_option('test1', source_option_test1), &
            text_option('phiext_full', source_option_phiext_full), &
            text_option('test2_full', source_option_test2_full) /)
    character(20) :: source_option
    namelist /dist_fn_knobs/ kperiod, gridfac, apfac, driftknob, poisfac
    namelist /source_knobs/ t0, omega0, gamma0, gamma_damp, gamma_damp_e, &
         thetas, phi_ext, source_option, thetas
    integer :: ierr
    kperiod = 0
    poisfac = 0.
    gridfac = 5e4
    apfac = 1.0
    driftknob = 1.0
    t0 = 100.0
    omega0 = 0.0
    gamma0 = 0.0
    gamma_damp = 0.0
    gamma_damp_e = 0.0
    thetas = 1.0
    phi_ext = 0.0
    source_option = 'default'
    read (unit=input_unit("dist_fn_knobs"), nml=dist_fn_knobs)
    read (unit=input_unit("source_knobs"), nml=source_knobs)

    ierr = error_unit()
    call get_option_value &
         (source_option, sourceopts, source_option_switch, &
          ierr, "source_option in source_knobs")

    call read_species_knobs
  end subroutine read_parameters

  subroutine read_species_knobs
    use species, only: nspec
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer :: is, unit

    allocate (fexp(nspec), bkdiff(nspec))
    do is = 1, nspec
       fexp(is) = (0.4,0.0)
       bkdiff(is) = 0.0
       call get_indexed_namelist_unit (unit, "dist_fn_species_knobs", is)
       call fill_species_knobs (unit, fexp(is), bkdiff(is))
       close (unit=unit)
    end do
  end subroutine read_species_knobs

  subroutine fill_species_knobs (unit, fexp_out, bakdif_out)
    implicit none
    integer, intent (in) :: unit
    complex, intent (in out) :: fexp_out
    real, intent (in out) :: bakdif_out
    real :: fexpr, fexpi, bakdif
    namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif

    fexpr = real(fexp_out)
    fexpi = aimag(fexp_out)
    bakdif = bakdif_out
    read (unit=unit, nml=dist_fn_species_knobs)
    fexp_out = cmplx(fexpr,fexpi)
    bakdif_out = bakdif
  end subroutine fill_species_knobs

  subroutine broadcast_parameters
    use mp, only: proc0, broadcast
    use species, only: nspec
    implicit none

    call broadcast (kperiod)
    call broadcast (gridfac)
    call broadcast (poisfac)
    call broadcast (apfac)
    call broadcast (driftknob)
    call broadcast (t0)
    call broadcast (phi_ext)
    call broadcast (omega0)
    call broadcast (gamma0)
    call broadcast (gamma_damp)
    call broadcast (gamma_damp_e)
    call broadcast (thetas)
    call broadcast (source_option_switch)

    if (.not. proc0) allocate (fexp(nspec), bkdiff(nspec))
    call broadcast (fexp)
    call broadcast (bkdiff)
  end subroutine broadcast_parameters

  subroutine init_wdrift
    use species, only: nspec
    use theta_grid, only: ntgrid, bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda, e, al, jend
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo

    allocate (ittp(-ntgrid:ntgrid))
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

    allocate (wdrift(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (wdriftttp(-ntgrid:ntgrid,naky,ntheta0,negrid,nspec))
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
                      + (gbdrift(ig) + theta0(ik,it)*gbdrift(ig)) &
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

    allocate (vpa(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (vpac(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (vperp2(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (vpar(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))

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

    allocate (wstar(naky,negrid,nspec))

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
    use dist_fn_arrays, only: aj0, aj1, kperp2
    use species, only: nspec, spec, smz
    use theta_grid, only: ntgrid, bmag, gds2, gds21, gds22, shat
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    use le_grids, only: e, al
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo
    real :: arg

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

    allocate (aj0(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))

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

    allocate (a(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (b(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (r(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (ainv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_proc))

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
               = 1.0/(1.0 + bkdiff(is) &
                   + (1.0-fexp(is))*tz(spec,is)*(zi*wd + 2.0*vp))
          r(ig,iglo) &
               = (1.0 - bkdiff(is) &
                  + (1.0-fexp(is))*tz(spec,is)*(zi*wd - 2.0*vp)) &
                 *ainv(ig,iglo)
          a(ig,iglo) &
               = 1.0 + bkdiff(is) + fexp(is)*tz(spec,is)*(-zi*wd - 2.0*vp)
          b(ig,iglo) &
               = 1.0 - bkdiff(is) + fexp(is)*tz(spec,is)*(-zi*wd + 2.0*vp)

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
  end subroutine init_invert_rhs

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo
    implicit none

    allocate (g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (gnew(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
    allocate (g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
  end subroutine allocate_arrays

  subroutine timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
    use theta_grid, only: ntgrid
    use collisions, only: solfp1
    use dist_fn_arrays, only: gnew
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep

    call invert_rhs (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
    call solfp1 (gnew, g0, phinew, aparnew, aperpnew)
  end subroutine timeadv

  subroutine invert_rhs (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2, aj0, aj1, g, gnew
    use species, only: nspec, spec, stm, zstm, tz
    use theta_grid, only: ntgrid, theta
    use le_grids, only: nlambda, ng2, lmax, anon, forbid, e
    use kt_grids, only: aky
    use run_parameters, only: fphi, fapar, faperp, delt, wunits, tunits, delt
    use constants
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, aperpnew
    integer, intent (in) :: istep

    integer :: iglo
    integer :: ig, ik, it, il, ie, isgn, is
    integer :: ilmin
    complex :: beta1
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg
    complex, dimension (-ntgrid:ntgrid,2) :: source, g1, g2

    logical :: kperiod_flag
    real :: time, sourcefac1
    complex :: sourcefac2, sourcefac

    time = (istep-1)*delt
    if (time > t0) then
       sourcefac1 = 1.0
    else
       sourcefac1 = 0.5 - 0.5*cos(pi*time/t0)
    end if
    sourcefac2 = exp(-zi*omega0*time+gamma0*time)
    sourcefac = sourcefac1*sourcefac2

    call prof_entering ("invert_rhs", "dist_fn")
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       ! gyroaveraged phi and apar
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
             do isgn = 1, 2
                do ig = -ntgrid, ntgrid-1
                   source(ig,isgn) = anon(ie,is) &
                        *(-2.0*vpar(ig,isgn,iglo) &
                            *(phigavg(ig+1) - phigavg(ig)) &
                          -zstm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                            *(aj0(ig+1,iglo) + aj0(ig,iglo))*0.5 &
                            *(aparnew(ig+1,ik,it) + aparnew(ig,ik,it) &
                              - apar(ig+1,ik,it) - apar(ig,ik,it)) &
                          -zi*wdrift(ig,iglo)*(phigavg(ig+1) + phigavg(ig))) &
                      + zi*(wstar(ik,ie,is) &
                            + vpac(ig,isgn,iglo)*delt*wunits(ik) &
                              *(2.0*spec(is)%uprim &
			+ spec(is)%uprim2*e(ie,is)**(1.5)*sqrt(pi)/4.0)) &
                          *(phigavg(ig+1) + phigavg(ig) &
                            - stm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                              *(apargavg(ig+1) + apargavg(ig)))
                end do
             end do
          else
             source = 0.0
          end if
       case (source_option_zero)
          source = 0.0
       case (source_option_sine)
          source(:ntgrid-1,1) = tunits(ik)*delt &
               *(sin(theta(:ntgrid-1)) + sin(theta(-ntgrid+1:)))
          source(:ntgrid-1,2) = source(:ntgrid-1,1)
       case (source_option_test1)
          do ig = -ntgrid, ntgrid-1
             source(ig,1) = -zi*(wdrift_func(ig,ik,it,il,ie,is) &
                  + wdrift_func(ig+1,ik,it,il,ie,is)) &
                  *sourcefac
          end do
          source(:ntgrid-1,2) = source(:ntgrid-1,1)
       case (source_option_phiext_full)
          if (il <= lmax) then
             if (istep > 0 .and. aky(ik) == 0.0) then
                do isgn = 1, 2
                   do ig = -ntgrid, ntgrid-1
                      source(ig,isgn) = anon(ie,is) &
                           *(-2.0*vpar(ig,isgn,iglo) &
                               *(phigavg(ig+1) - phigavg(ig)) &
                             -zstm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                               *(aj0(ig+1,iglo) + aj0(ig,iglo))*0.5 &
                               *(aparnew(ig+1,ik,it) + aparnew(ig,ik,it) &
                                 - apar(ig+1,ik,it) - apar(ig,ik,it)) &
                             -zi*wdrift(ig,iglo) &
                              *(phigavg(ig+1) + phigavg(ig) &
                                + 2.0*phi_ext*sourcefac)) &
                         + zi*(wstar(ik,ie,is) &
                               + vpac(ig,isgn,iglo)*delt*wunits(ik) &
                                 *(2.0*spec(is)%uprim &
			+ spec(is)%uprim2*e(ie,is)**(1.5)*sqrt(pi)/4.0)) &
                             *(phigavg(ig+1) + phigavg(ig) &
                               + 2.0*phi_ext*sourcefac &
                               - stm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                                 *(apargavg(ig+1) + apargavg(ig)))
                   end do
                end do
             else
                do isgn = 1, 2
                   do ig = -ntgrid, ntgrid-1
                      source(ig,isgn) = anon(ie,is) &
                           *(-2.0*vpar(ig,isgn,iglo) &
                               *(phigavg(ig+1) - phigavg(ig)) &
                             -zstm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                               *(aj0(ig+1,iglo) + aj0(ig,iglo))*0.5 &
                               *(aparnew(ig+1,ik,it) + aparnew(ig,ik,it) &
                                 - apar(ig+1,ik,it) - apar(ig,ik,it)) &
                             -zi*wdrift(ig,iglo)*(phigavg(ig+1)+phigavg(ig))) &
                         + zi*(wstar(ik,ie,is) &
                               + vpac(ig,isgn,iglo)*delt*wunits(ik) &
                                 *(2.0*spec(is)%uprim &
			+ spec(is)%uprim2*e(ie,is)**(1.5)*sqrt(pi)/4.0)) &
                             *(phigavg(ig+1) + phigavg(ig) &
                               - stm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                                 *(apargavg(ig+1) + apargavg(ig)))
                   end do
                end do
             end if
          else
             source = 0.0
          end if
       case (source_option_test2_full)
          if (il <= lmax) then
             do isgn = 1, 2
                do ig = -ntgrid, ntgrid-1
                   source(ig,isgn) = anon(ie,is) &
                        *(-2.0*vpar(ig,isgn,iglo) &
                            *(phigavg(ig+1) - phigavg(ig)) &
                          -zstm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                            *(aj0(ig+1,iglo) + aj0(ig,iglo))*0.5 &
                            *(aparnew(ig+1,ik,it) + aparnew(ig,ik,it) &
                              - apar(ig+1,ik,it) - apar(ig,ik,it)) &
                          -zi*wdrift(ig,iglo)*(phigavg(ig+1) + phigavg(ig))) &
                      + zi*(wstar(ik,ie,is) &
                            + vpac(ig,isgn,iglo)*wunits(ik) &
                              *2.0*spec(is)%uprim*delt) &
                          *(phigavg(ig+1) + phigavg(ig) &
                            - stm(spec,is)*vpac(ig,isgn,iglo)*wunits(ik) &
                              *(apargavg(ig+1) + apargavg(ig)))
                end do
             end do
          else
             source = 0.0
          end if
          if (aky(ik) == 0.0 .and. istep > 0) then
             source(:ntgrid-1,1) = source(:ntgrid-1,1) &
                  + tunits(ik)*delt &
                    *(aj0(:ntgrid-1,iglo)**2 + aj0(-ntgrid+1:,iglo)**2) &
                    *(e(ie,is)-1.5)*sourcefac &
                    *exp(-(theta(:ntgrid-1)/thetas)**2)
             source(:ntgrid-1,2) = source(:ntgrid-1,2) &
                  + tunits(ik)*delt &
                    *(aj0(:ntgrid-1,iglo)**2 + aj0(-ntgrid+1:,iglo)**2) &
                    *(e(ie,is)-1.5)*sourcefac &
                    *exp(-(theta(:ntgrid-1)/thetas)**2)
          end if
       end select

       if (il <= lmax) then
          do ig = -ntgrid, ntgrid-1
             source(ig,1) = source(ig,1) &
                  + b(ig,iglo)*g(ig,1,iglo) + a(ig,iglo)*g(ig+1,1,iglo)
             source(ig,2) = source(ig,2) &
                  + a(ig,iglo)*g(ig,2,iglo) + b(ig,iglo)*g(ig+1,2,iglo)
          end do
       end if

       source(ntgrid,:) = source(-ntgrid,:)

       ! special source term for totally trapped particles
       if (nlambda > ng2) then
          do ig = -ntgrid, ntgrid
             if (il /= ittp(ig)) cycle
             source(ig,2) &
                  = g(ig,2,iglo)*a(ig,iglo) &
                    - anon(ie,is)*(zi*wdriftttp(ig,ik,it,ie,is)*phigavg(ig)) &
                    + zi*wstar(ik,ie,is)*phigavg(ig)
          end do
       end if

       ! gnew is the inhomogeneous solution
       gnew(:,:,iglo) = 0.0

       ! g1 is the homogeneous solution
       g1 = 0.0
       kperiod_flag = kperiod == 1 .or. aky(ik) == 0.0
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

       ! numerical knob games
       gnew(:,:,iglo) = gnew(:,:,iglo) &
            /(1.0 + gamma_damp*delt/(1.0 + gamma_damp_e*e(ie,is)))
    end do
    call prof_leaving ("invert_rhs", "dist_fn")
  end subroutine invert_rhs

  subroutine getan (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2, aj0, aj1, gnew
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
       g0(:,isgn,:) = aj0*gnew(:,isgn,:)
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g0, wgt, antot)

    do isgn = 1, 2
       g0(:,isgn,:) = aj0*vpa(:,isgn,:)*gnew(:,isgn,:)
    end do
    wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
    call integrate_species (g0, wgt, antota)

    do isgn = 1, 2
       g0(:,isgn,:) = aj1*vperp2*gnew(:,isgn,:)
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g0, wgt, antotp)
    call prof_leaving ("getan", "dist_fn")
  end subroutine getan

  subroutine init_fieldeq
    use dist_fn_arrays, only: aj0, aj1, vperp2, kperp2
    use species, only: nspec, spec, stm, electron_species
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: anon, integrate_species
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use run_parameters, only: teti
    implicit none
    integer :: iglo, isgn
    integer :: ik, ie, is
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: tot
    real, dimension (nspec) :: wgt

    allocate (gridfac1(-ntgrid:ntgrid,naky))
    gridfac1 = 1.0
    if (kperiod /= 1) then
       do ik = 1, naky
          if (aky(ik) == 0.0) cycle
          gridfac1(-ntgrid,ik) = gridfac
          gridfac1(ntgrid,ik) = gridfac
       end do
    end if

    allocate (gamtot(-ntgrid:ntgrid,naky,ntheta0))
    allocate (gamtot1(-ntgrid:ntgrid,naky,ntheta0))
    allocate (gamtot2(-ntgrid:ntgrid,naky,ntheta0))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g0(:,isgn,iglo) = (1.0 - aj0(:,iglo)**2)*anon(ie,is)
       end do
    end do
    wgt = spec%z*spec%z*spec%dens/spec%temp
    call integrate_species (g0, wgt, tot)
    gamtot = real(tot) + poisfac*kperp2
    
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

    if (.not. any(spec%type == electron_species)) then
       ! adiabatic electrons
       do ik = 1, naky
          if (aky(ik) == 0.0) cycle
          gamtot = gamtot + teti
       end do
    end if
  end subroutine init_fieldeq

  subroutine getfieldeq1 (phi, apar, aperp, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use dist_fn_arrays, only: kperp2
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0, aky, theta0
    use run_parameters, only: beta
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    integer :: ik, it

    do it = 1, ntheta0
       do ik = 1, naky
          fieldeq(:,ik,it) &
               = antot(:,ik,it) + aperp(:,ik,it)*gamtot1(:,ik,it) &
                 - gamtot(:,ik,it)*gridfac1(:,ik)*phi(:,ik,it)
          fieldeqa(:,ik,it) &
               = antota(:,ik,it) &
                 - kperp2(:,ik,it)*gridfac1(:,ik)*apar(:,ik,it)
          fieldeqp(:,ik,it) &
               = (antotp(:,ik,it) + aperp(:,ik,it)*gamtot2(:,ik,it) &
                  + 0.5*phi(:,ik,it)*gamtot1(:,ik,it))*beta*apfac/bmag**2 &
                 + aperp(:,ik,it)*gridfac1(:,ik)
       end do
    end do
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

  pure function j0 (x)
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
       pflux, qheat, vflux, pmflux, qmheat, vmflux)
    use species, only: nspec, spec, stm
    use theta_grid, only: ntgrid, bmag, gradpar
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: nlambda, negrid, e
    use dist_fn_arrays, only: g, aj0, vpac, vpa
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use mp, only: proc0
    use run_parameters, only: woutunits
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, aperp
    real, dimension (:,:,:), intent (out) :: pflux, qheat, vflux
    real, dimension (:,:,:), intent (out) :: pmflux, qmheat, vmflux
    real, dimension (-ntgrid:ntgrid,naky) :: dnorm
    real :: anorm
    integer :: ig, ik, is, isgn
    integer :: iglo

    if (proc0) then
       pflux = 0.0;   qheat = 0.0;   vflux = 0.0
       pmflux = 0.0;  qmheat = 0.0;  vmflux = 0.0
    end if

    do ik = 1, naky
       do ig = -ntgrid, ntgrid
          dnorm(ig,ik) = 1.0/gradpar(ig)/bmag(ig)*woutunits(ik)/sqrt(2.0)
       end do
    end do
    dnorm(-ntgrid,:) = 0.5*dnorm(-ntgrid,:)
    dnorm(ntgrid,:)  = 0.5*dnorm(ntgrid,:)

    anorm = sum(spread(dnorm,3,ntheta0)*abs(phi)**2 + abs(apar)**2)

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
    use kt_grids, only: naky, ntheta0
    use le_grids, only: integrate
    use gs2_layouts, only: geint_lo, ik_idx, it_idx, is_idx, proc_id, idx_local
    use mp, only: proc0, send, receive, barrier
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (out) :: flx
    real :: anorm
    real, dimension (-ntgrid:,:) :: dnorm
    complex, &
         dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_proc) :: &
         geint
    real :: tmp
    integer :: igeint
    integer :: ik, it, is

    if (abs(anorm) < 2.0*epsilon(0.0)) then
       if (proc0) flx = 0.0
       return
    end if
    call integrate (g0, geint)
    do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       tmp = sum(-aimag(geint(:,igeint)*conjg(fld(:,ik,it)))*dnorm(:,ik))/anorm
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

  subroutine neoclassical_flux (pflux, qflux, istep)
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
    integer, intent (in) :: istep
    real, dimension (-ntgrid:ntgrid,naky) :: dnorm
    complex, dimension (naky,ntheta0,nspec) :: anorm
    complex, &
         dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_proc) :: &
         geint
    integer :: ig, ik, it, ie, is
    integer :: iglo, igeint
    complex :: x

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
                     = sum(geint(:ntgrid-1,igeint)*delthet(:ntgrid-1))
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
                        = sum(geint(:ntgrid-1,igeint) &
                              *delthet(:ntgrid-1))/anorm(ik,it,is)
                else
                   call receive (pflux(ik,it,is), proc_id(geint_lo,igeint))
                end if
             else if (idx_local(geint_lo,igeint)) then
                x = sum(geint(:ntgrid-1,igeint)*delthet(:ntgrid-1)) &
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
                        = sum(geint(:ntgrid-1,igeint) &
                              *delthet(:ntgrid-1))/anorm(ik,it,is)
                else
                   call receive (qflux(ik,it,is), proc_id(geint_lo,igeint))
                end if
             else if (idx_local(geint_lo,igeint)) then
                x = sum(geint(:ntgrid-1,igeint)*delthet(:ntgrid-1)) &
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
                     kperp2(ig,ik,it)*gridfac1(ig,ik)*apar(ig,ik,it), &
                     fieldeqa(ig,ik,it), gamtot1(ig,ik,it)
                write (unit,"(8(1pe12.5,x))") &
                     gridfac1(ig,ik)*aperp(ig,ik,it), fieldeqp(ig,ik,it), &
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
end module dist_fn
