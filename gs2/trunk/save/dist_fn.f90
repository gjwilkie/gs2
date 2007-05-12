module dist_fn
  use init_g, only: ginit
  use redistribute, only: redist_type
  implicit none
  public :: init_dist_fn
  public :: timeadv, get_stress, exb_shear
  public :: getfieldeq, getan, getfieldexp, getmoms
  public :: flux, neoclassical_flux, lambda_flux
  public :: ginit, get_epar, e_flux, get_heat
  public :: init_intcheck, init_vortcheck, init_fieldcheck
  public :: intcheck, vortcheck, fieldcheck
  public :: finish_intcheck, finish_vortcheck, finish_fieldcheck
  public :: t0, omega0, gamma0, thetas, k0, nperiod_guard, source0
  public :: reset_init, write_f, reset_physics
  public :: M_class, N_class, i_class, par_spectrum
  public :: l_links, r_links, itright, itleft, boundary
  public :: init_kperp2
  public :: get_dens_vel, get_jext !GGH

  private

  ! knobs
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp ! nspec
  real :: gridfac, apfac, driftknob, poisfac
  real :: t0, omega0, gamma0, thetas, k0, source0
  real :: phi_ext, afilter, kfilter, a_ext
  real :: aky_star, akx_star
  real :: D_kill, noise, cfac, wfb, g_exb

  integer :: adiabatic_option_switch!, heating_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4, &
       adiabatic_option_noJ = 5
!  integer, parameter ::  heating_option_cowley = 1
  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_zero = 2, source_option_sine = 3, &
       source_option_test1 = 4, source_option_phiext_full = 5, &
       source_option_test2_full = 6, source_option_cosine = 7, &
       source_option_convect_full = 8, source_option_hm_force = 9
  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_alternate_zero = 3, &
       boundary_option_linked = 4
  logical :: mult_imp, test, def_parity, even
  logical :: save_n, save_u, save_Tpar, save_Tperp, flr
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false., kill_grid, h_kill
  integer :: nperiod_guard
  
!! k_parallel filter items
!  real, dimension(:), allocatable :: work, tablekp
!  real :: scale
!  integer :: nwork, ntablekp

  ! internal arrays

  real, dimension (:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:,:,:), allocatable :: wdriftttp
  ! (-ntgrid:ntgrid,ntheta0,naky,negrid,nspec) replicated

  real, dimension (:,:,:), allocatable :: wstar
  ! (naky,negrid,nspec) replicated

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  complex, dimension (:,:), allocatable :: a, b, r, ainv
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:), allocatable :: gridfac1
  ! (-ntgrid:ntgrid,ntheta0,naky)

  complex, dimension (:,:,:), allocatable :: g0, g_h
  ! (-ntgrid:ntgrid,2, -g-layout-)

  complex, dimension (:,:,:), allocatable :: g_adj
  ! (N(links), 2, -g-layout-)

  complex, dimension (:,:,:), allocatable, save :: gnl_1, gnl_2
  ! (-ntgrid:ntgrid,2, -g-layout-)

  ! momentum conservation
  complex, dimension (:,:), allocatable :: g3int
  real, dimension (:,:,:), allocatable :: sq

  ! connected bc

  integer, dimension (:,:), allocatable :: itleft, itright
  ! (naky,ntheta0)

  type :: connections_type
     integer :: iproc_left,  iglo_left
     integer :: iproc_right, iglo_right
     logical :: neighbor
  end type connections_type

  type (connections_type), dimension (:), allocatable :: connections
  ! (-g-layout-)

  ! linked only
  type (redist_type), save :: gc_from_left, gc_from_right
  type (redist_type), save :: links_p, links_h
  type (redist_type), save :: wfb_p, wfb_h
  integer, dimension (:,:), allocatable :: l_links, r_links
  integer, dimension (:,:,:), allocatable :: n_links
  logical, dimension (:,:), allocatable :: save_h
  logical :: no_comm = .false.
  integer, dimension(:), allocatable :: M_class, N_class
  integer :: i_class

  logical :: initialized = .false.
  logical :: initializing = .true.

contains

  subroutine init_dist_fn
    use mp, only: proc0, finish_mp
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0, akx, aky
    use le_grids, only: init_le_grids, nlambda, negrid
    use run_parameters, only: init_run_parameters
    use collisions, only: init_collisions
    use gs2_layouts, only: init_dist_fn_layouts, init_gs2_layouts
    use nonlinear_terms, only: init_nonlinear_terms
    use init_g, only: init_init_g
    use hyper, only: init_hyper
    implicit none

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts
    call init_species
    call init_theta_grid
    call init_kt_grids
    call init_le_grids (accelerated_x, accelerated_v)
    call read_parameters

    if (test) then
       if (proc0) then
          write (*,*) 'nspecies = ',nspec
          write (*,*) 'nlambda = ', nlambda
          write (*,*) 'negrid = ',negrid
          write (*,*) 'ntheta0 = ',ntheta0
          write (*,*) 'naky = ',naky
       end if
       call finish_mp
       stop
    end if

    call init_run_parameters
    call init_kperp2
    call init_collisions ! needs to be after init_run_parameters
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
    call init_init_g 
    call init_nonlinear_terms 
    call allocate_arrays
    call init_vpar
    call init_wdrift
    call init_wstar
    call init_bessel
    call init_par_filter
    call init_invert_rhs
    call init_fieldeq
    call init_hyper

  end subroutine init_dist_fn

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid, only: nperiod, shat
    use init_g, only: init_g_k0 => k0
    use text_options
    use species, only: nspec
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (10), parameter :: sourceopts = &
         (/ text_option('default', source_option_full), &
            text_option('full', source_option_full), &
            text_option('zero', source_option_zero), &
            text_option('sine', source_option_sine), &
            text_option('cosine', source_option_cosine), &
            text_option('test1', source_option_test1), &
            text_option('hm', source_option_hm_force), &
            text_option('phiext_full', source_option_phiext_full), &
            text_option('test2_full', source_option_test2_full), &
            text_option('convect_full', source_option_convect_full) /)
    character(20) :: source_option

    type (text_option), dimension (8), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('kperiod=1', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked), &
            text_option('alternate-zero', boundary_option_alternate_zero) /)
    character(20) :: boundary_option

!    type (text_option), dimension (3), parameter :: heatingopts = &
!         (/ text_option('default', heating_option_qdh), &
!            text_option('qdh', heating_option_qdh), &
!            text_option('Cowley',  heating_option_cowley) /)
!    character(20) :: heating_option

    type (text_option), dimension (8), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
! eventually add in iphi00 = 0 option:
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg), &
            text_option('dimits', adiabatic_option_noJ) /)
    character(30) :: adiabatic_option
            
    namelist /dist_fn_knobs/ boundary_option, gridfac, apfac, driftknob, &
         nperiod_guard, poisfac, adiabatic_option, &
         kfilter, afilter, mult_imp, test, def_parity, even, wfb, &
         save_n, save_u, save_Tpar, save_Tperp, D_kill, noise, flr, cfac, &
         kill_grid, h_kill, g_exb
    
    namelist /source_knobs/ t0, omega0, gamma0, source0, &
           thetas, k0, phi_ext, source_option, a_ext, aky_star, akx_star
    integer :: ierr, is, in_file
    logical :: exist
    real :: bd
    logical :: done = .false.

    if (done) return
    done = .true.

    if (proc0) then
       cfac = 1.0
       flr = .true.
       save_n = .true.
       save_u = .true.
       save_Tpar = .true.
       save_Tperp = .true.
       boundary_option = 'default'
       adiabatic_option = 'default'
       poisfac = 0.0
       gridfac = 1.0  ! used to be 5.e4
       apfac = 1.0
       driftknob = 1.0
       t0 = 100.0
       source0 = 1.0
       omega0 = 0.0
       gamma0 = 0.0
       thetas = 1.0
       k0 = init_g_k0
       aky_star = 0.0
       akx_star = 0.0
       phi_ext = 0.0
       a_ext = 0.0
       afilter = 0.0
       kfilter = 0.0
       g_exb = 0.0
       h_kill = .true.
       D_kill = -10.0
       noise = -1.
       wfb = 1.
       mult_imp = .false.
       test = .false.
       def_parity = .false.
       even = .true.
       source_option = 'default'
       kill_grid = .false.
       in_file = input_unit_exist("dist_fn_knobs", exist)
       if (exist) read (unit=input_unit("dist_fn_knobs"), nml=dist_fn_knobs)
       in_file = input_unit_exist("source_knobs", exist)
       if (exist) read (unit=input_unit("source_knobs"), nml=source_knobs)

       if(abs(shat) <=  1.e-5) boundary_option = 'periodic'

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
!       call get_option_value &
!            (heating_option, heatingopts, heating_option_switch, &
!            ierr, "heating_option in dist_fn_knobs")

       nperiod_guard = 0
    end if
    if (.not.allocated(fexp)) allocate (fexp(nspec), bkdiff(nspec), bd_exp(nspec))
    if (proc0) call read_species_knobs

    call broadcast (cfac)
    call broadcast (kill_grid)
    call broadcast (flr)
    call broadcast (save_n)
    call broadcast (save_u)
    call broadcast (save_Tpar)
    call broadcast (save_Tperp)
    call broadcast (boundary_option_switch)
    call broadcast (adiabatic_option_switch)
!    call broadcast (heating_option_switch)
    call broadcast (gridfac)
    call broadcast (poisfac)
    call broadcast (apfac)
    call broadcast (driftknob)
    call broadcast (t0)
    call broadcast (source0)
    call broadcast (omega0)
    call broadcast (gamma0)
    call broadcast (thetas)
    call broadcast (k0)
    call broadcast (aky_star)
    call broadcast (akx_star)
    call broadcast (phi_ext)
    call broadcast (a_ext)
    call broadcast (D_kill)
    call broadcast (h_kill)
    call broadcast (g_exb)
    call broadcast (noise)
    call broadcast (afilter)
    call broadcast (kfilter)
    call broadcast (source_option_switch)
    call broadcast (fexp)
    call broadcast (bkdiff)
    call broadcast (bd_exp)
    call broadcast (nperiod_guard)
    call broadcast (mult_imp)
    call broadcast (test)
    call broadcast (def_parity)
    call broadcast (even)
    call broadcast (wfb)

    if (mult_imp) then
       ! nothing -- fine for linear runs, but not implemented nonlinearly
    else
! consistency check for bkdiff
       bd = bkdiff(1)
       do is = 1, nspec
          if (bkdiff(is) /= bd) then
             if (proc0) write(*,*) 'Forcing bkdiff for species ',is,' equal to ',bd
             if (proc0) write(*,*) 'If this is a linear run, and you want unequal bkdiff'
             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
             bkdiff(is) = bd
          endif
       end do
! consistency check for fexp
!       fe = fexp(1)
!       do is = 1, nspec
!          if (fexp(is) /= fe) then
!             if (proc0) write(*,*) 'Forcing fexp for species ',is,' equal to ',fe
!             if (proc0) write(*,*) 'If this is a linear run, and you want unequal fexp'
!             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
!             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
!             fexp(is) = fe
!          endif
!       end do
    end if

! consistency check for afilter
!    if (afilter /= 0.0) then
!       if (proc0) write(*,*) 'Forcing afilter = 0.0'
!       afilter = 0.0
!    end if

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
    integer :: bd_exp
    real :: fexpr, fexpi, bakdif
    namelist /dist_fn_species_knobs/ fexpr, fexpi, bakdif, bd_exp

    fexpr = real(fexp_out)
    fexpi = aimag(fexp_out)
    bakdif = bakdif_out
    bd_exp = bd_exp_out
    read (unit=unit, nml=dist_fn_species_knobs)
    fexp_out = cmplx(fexpr,fexpi)
    bd_exp_out = bd_exp
    bakdif_out = bakdif
  end subroutine fill_species_knobs

  subroutine init_wdrift
    use species, only: nspec
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use le_grids, only: negrid, nlambda, al, jend
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use dist_fn_arrays, only: ittp
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo, ierr
    logical :: alloc = .true.

! find totally trapped particles 
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
    if (alloc) allocate (wdriftttp(-ntgrid:ntgrid,ntheta0,naky,negrid,nspec))
    wdrift = 0.  ; wdriftttp = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid
          wdrift(ig,iglo) &
               = wdrift_func(ig, il_idx(g_lo,iglo), ie_idx(g_lo,iglo), &
                                 it_idx(g_lo,iglo), ik_idx(g_lo,iglo), &
                                 is_idx(g_lo,iglo))
       end do
    end do
    wdriftttp = 0.0
    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             do it = 1, ntheta0
! moved this loop inside. 4.10.99
                do ig = -ntgrid, ntgrid
                   if (ittp(ig) == 0) cycle
                   wdriftttp(ig,it,ik,ie,is) &
                        = wdrift_func(ig,ittp(ig),ie,it,ik,is)*driftknob
                end do
             end do
          end do
       end do
    end do

! This should be weighted by bakdif to be completely consistent
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid-1
          wdrift(ig,iglo) = 0.5*(wdrift(ig,iglo) + wdrift(ig+1,iglo))*driftknob
       end do
    end do

    alloc = .false.

  end subroutine init_wdrift

  function wdrift_func (ig, il, ie, it, ik, is)
    use theta_grid, only: bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use theta_grid, only: shat
    use kt_grids, only: aky, theta0, akx
    use le_grids, only: e, al
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    implicit none
    real :: wdrift_func
    integer, intent (in) :: ig, ik, it, il, ie, is

    if (aky(ik) == 0.0) then
       wdrift_func = akx(it)/shat &
                    *(cvdrift0(ig)*e(ie,is)*(1.0 - al(il)*bmag(ig)) &
                      + gbdrift0(ig)*0.5*e(ie,is)*al(il)*bmag(ig)) &
                     *code_dt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
                        *e(ie,is)*(1.0 - al(il)*bmag(ig)) &
                      + (gbdrift(ig) + theta0(it,ik)*gbdrift0(ig)) &
                        *0.5*e(ie,is)*al(il)*bmag(ig)) &
                     *code_dt*wunits(ik)
    end if
  end function wdrift_func

  subroutine init_vpar
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2
    use species, only: spec
    use theta_grid, only: ntgrid, delthet, bmag, gradpar
    use le_grids, only: e, al
    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo, ik, is, ie, il
    real :: al1, e1
    

    if (.not.allocated(vpa)) then
       allocate (vpa    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpac   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vperp2 (-ntgrid:ntgrid,  g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpar   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
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

! Where vpac /= 1, it could be weighted by bakdif for better consistency??
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
            spec(is)%zstm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,1,iglo)
       vpar(-ntgrid:ntgrid-1,2,iglo) = &
            spec(is)%zstm*tunits(ik)*code_dt &
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
    use run_parameters, only: wunits
    use gs2_time, only: code_dt
    implicit none
    integer :: ik, ie, is

    if(.not.allocated(wstar)) allocate (wstar(naky,negrid,nspec))

    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             wstar(ik,ie,is) = code_dt*wunits(ik) &
                  *(spec(is)%fprim+spec(is)%tprim*(e(ie,is)-1.5))
          end do
       end do
    end do
  end subroutine init_wstar

  subroutine init_bessel
    use dist_fn_arrays, only: aj0, aj1, kperp2
    use species, only: spec
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

    call init_kperp2

    allocate (aj0(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    aj0 = 0. ; aj1 = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          arg = spec(is)%smz*sqrt(e(ie,is)*al(il)/bmag(ig)*kperp2(ig,it,ik))
          aj0(ig,iglo) = j0(arg)
          aj1(ig,iglo) = j1(arg)
       end do
    end do

  end subroutine init_bessel

  subroutine init_kperp2
    use dist_fn_arrays, only: kperp2
    use species, only: spec
    use theta_grid, only: ntgrid, bmag, gds2, gds21, gds22, shat
    use kt_grids, only: naky, ntheta0, aky, theta0, akx
    implicit none
    integer :: ik, it
    logical :: done = .false.

    if (done) return
    done = .true.

    allocate (kperp2(-ntgrid:ntgrid,ntheta0,naky))
    do ik = 1, naky
       if (aky(ik) == 0.0) then
         do it = 1, ntheta0
             kperp2(:,it,ik) = akx(it)*akx(it)*gds22/(shat*shat)
          end do
       else
          do it = 1, ntheta0
             kperp2(:,it,ik) = aky(ik)*aky(ik) &
                  *(gds2 + 2.0*theta0(it,ik)*gds21 &
                  + theta0(it,ik)*theta0(it,ik)*gds22)
          end do
       end if
    end do

  end subroutine init_kperp2

  subroutine init_par_filter
    use theta_grid, only: ntgrid, nperiod
    use gs2_transforms, only: init_zf

    call init_zf (ntgrid, nperiod)

  end subroutine init_par_filter

  subroutine par_spectrum(an, an2)

    use gs2_transforms, only: kz_spectrum
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    complex, dimension(:,:,:) :: an, an2    
    integer :: it, ik
    real :: scale

    call kz_spectrum (an, an2, ntgrid, ntheta0, naky)
    scale = 1./real(4*ntgrid**2)
    an2 = an2*scale

  end subroutine par_spectrum

  subroutine init_invert_rhs
    use mp, only: proc0
    use dist_fn_arrays, only: vpa, vpar, vpac, ittp
    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: negrid, nlambda, ng2, forbid
    use constants
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is
    real :: wd, wdttp, vp, bd

    if (.not.allocated(a)) then
       allocate (a(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (b(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (r(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (ainv(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    a = 0. ; b = 0. ; r = 0. ; ainv = 0.
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid-1
          wd = wdrift(ig,iglo)
          wdttp = wdriftttp(ig,it,ik,ie,is)
          vp = vpar(ig,1,iglo)
          bd = bkdiff(is)

          ainv(ig,iglo) &
               = 1.0/(1.0 + bd &
               + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0+bd) + 2.0*vp))
          r(ig,iglo) &
               = (1.0 - bd &
               + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0-bd) - 2.0*vp)) &
               *ainv(ig,iglo)
          a(ig,iglo) &
               = 1.0 + bd &
               + fexp(is)*spec(is)%tz*(-zi*wd*(1.0+bd) - 2.0*vp)
          b(ig,iglo) &
               = 1.0 - bd &
               + fexp(is)*spec(is)%tz*(-zi*wd*(1.0-bd) + 2.0*vp)
          
          if (nlambda > ng2) then
             ! zero out forbidden regions
             if (forbid(ig,il) .or. forbid(ig+1,il)) then
                r(ig,iglo) = 0.0
                ainv(ig,iglo) = 0.0
             end if
             
             ! ???? mysterious mucking around at lower bounce point
             ! part of multiple trapped particle algorithm
             if (forbid(ig,il) .and. .not. forbid(ig+1,il)) then
                ainv(ig,iglo) = 1.0 + ainv(ig,iglo)
             end if
             
             ! ???? mysterious mucking around with totally trapped particles
             ! part of multiple trapped particle algorithm
             if (il == ittp(ig)) then
                ainv(ig,iglo) = 1.0/(1.0 + zi*(1.0-fexp(is))*spec(is)%tz*wdttp)
                a(ig,iglo) = 1.0 - zi*fexp(is)*spec(is)%tz*wdttp
                r(ig,iglo) = 0.0
             end if
          end if
       end do
    end do

    select case (boundary_option_switch)
    case (boundary_option_linked)
       call init_connected_bc
    case default
       if (.not. allocated(l_links)) then
          allocate (l_links(naky, ntheta0))
          allocate (r_links(naky, ntheta0))
          allocate (n_links(2, naky, ntheta0))
       end if
       l_links = 0;   r_links = 0;  n_links = 0
       
       i_class = 1
       if (.not. allocated(M_class)) then
          allocate (M_class(i_class))
          allocate (N_class(i_class))
       end if
       M_class = naky*ntheta0 ; N_class = 1
       
    end select

    initializing = .false.

  end subroutine init_invert_rhs

  subroutine init_connected_bc
    use theta_grid, only: ntgrid, nperiod, ntheta, theta
    use kt_grids, only: naky, ntheta0, aky, theta0
    use le_grids, only: ng2, nlambda
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id
    use mp, only: iproc, nproc, max_allreduce
    use constants
    use redistribute, only: index_list_type, init_fill, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_left, from_right
    type (index_list_type), dimension(0:nproc-1) :: to_right, from_left
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from_left, nn_to_right
    integer, dimension (0:nproc-1) :: nn_from_right, nn_to_left
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension (3) :: to_low, from_low, to_high, from_high
    integer :: ik, it, il, ie, is, iglo, it0, itl, itr, jshift0
    integer :: ip, ipleft, ipright
    integer :: iglo_left, iglo_right, i, j, k
    integer :: iglo_star, it_star, ncell
    integer :: n, isign, ig, n_links_max, nn_max
    integer :: ng
    integer, dimension(naky*ntheta0) :: n_k

    logical :: done = .false.

    if (done) return
    done = .true.

    ng = ntheta/2 + (nperiod-1)*ntheta

    if (naky > 1 .and. ntheta0 > 1) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,2)-theta0(1,2)) + 0.1)
    else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) /= 0.0) then
       jshift0 = int((theta(ng)-theta(-ng))/(theta0(2,1)-theta0(1,1)) + 0.1)
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
!                itl = ntheta0
!                itr = ntheta0
                itl = it0
                itr = it0
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

! if non-wfb trapped particle, no connections
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
       if (connections(iglo)%iproc_left >= 0 .or. &
            connections(iglo)%iproc_right >= 0) then
          connections(iglo)%neighbor = .true.
       else
          connections(iglo)%neighbor = .false.
       end if
    end do

    if (boundary_option_switch == boundary_option_linked) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc          
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          if (connections(iglo)%iglo_left >= 0 .and. aky(ik) /= 0.0) &
               save_h (1,iglo) = .true.
          if (connections(iglo)%iglo_right >= 0 .and. aky(ik) /= 0.0) &
               save_h (2,iglo) = .true.
! wfb (linked)
          if (nlambda > ng2               .and. &
               il == ng2+1                .and. &
               connections(iglo)%neighbor .and. &
               aky(ik) /= 0.0) &
               save_h (:,iglo) = .true.
       end do

       allocate (l_links(naky, ntheta0))
       allocate (r_links(naky, ntheta0))
       allocate (n_links(2, naky, ntheta0))

       n_links_max = 0
       do it = 1, ntheta0
          do ik = 1, naky
! count the links for each region
! l_links = number of links to the left
             l_links(ik, it) = 0
             it_star = it
             do 
                if (it_star == itleft(ik, it_star)) exit
                if (itleft(ik, it_star) >= 0) then
                   l_links(ik, it) = l_links(ik, it) + 1
                   it_star = itleft(ik, it_star)
                   if (l_links(ik, it) > 5000) then
! abort by deallocating twice
                      write(*,*) 'l_links error'
                      deallocate (l_links)
                      deallocate (l_links)
                   endif
                else
                   exit
                end if
             end do
! r_links = number of links to the right
             r_links(ik, it) = 0
             it_star = it
             do 
                if (it_star == itright(ik, it_star)) exit
                if (itright(ik, it_star) >= 0) then
                   r_links(ik, it) = r_links(ik, it) + 1
                   it_star = itright(ik, it_star)
                   if (r_links(ik, it) > 5000) then
! abort by deallocating twice
                      write(*,*) 'r_links error'
                      deallocate (r_links)
                      deallocate (r_links)
                   endif
                else
                   exit
                end if
             end do
! 'n_links' complex numbers are needed to specify bc for (ik, it) region
! ignoring wfb
! n_links(1,:,:) is for v_par > 0, etc.
             if (l_links(ik, it) == 0) then
                n_links(1, ik, it) = 0
             else 
                n_links(1, ik, it) = 2*l_links(ik, it) - 1
             end if

             if (r_links(ik, it) == 0) then
                n_links(2, ik, it) = 0
             else 
                n_links(2, ik, it) = 2*r_links(ik, it) - 1
             end if
             n_links_max = max(n_links_max, n_links(1,ik,it), n_links(2,ik,it))
          end do
       end do
! wfb
       if (n_links_max > 0) n_links_max = n_links_max + 3
       
! now set up communication pattern:
! excluding wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_star = iglo
          do j = 1, r_links(ik, it)
             call get_right_connection (iglo_star, iglo_right, ipright)
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1

             iglo_star = iglo_right
          end do
             
          iglo_star = iglo
          do j = 1, l_links(ik, it)
             call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             iglo_star = iglo_left
          end do

       end do
       
       nn_max = maxval(nn_to)
       call max_allreduce (nn_max)
       if (nn_max == 0) then
          no_comm = .true.
          goto 200
       end if

       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
          
          iglo_star = iglo
          do j = 1, r_links(ik, it)
             call get_right_connection (iglo_star, iglo_right, ipright)
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from(ipright)%first(n) = ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = j 
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
             iglo_star = iglo_right
          end do
             
          iglo_star = iglo
          do j = 1, l_links(ik, it)
             call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n
                from(ipleft)%first(n) = -ntgrid
                from(ipleft)%second(n) = 2
                from(ipleft)%third(n) = iglo
             end if
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = j 
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
             iglo_star = iglo_left
          end do
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       call init_fill (links_p, 'c', to_low, to_high, to, &
            from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)
       
! take care of wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
             
          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do
             
! v_par < 0:
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0: 
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from(ipright)%first(n) = ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = r_links(ik, it) + 1
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do
             
! v_par < 0: 
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n
                from(ipleft)%first(n) = -ntgrid
                from(ipleft)%second(n) = 2
                from(ipleft)%third(n) = iglo
             end if
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = l_links(ik, it) + 1
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       call init_fill (wfb_p, 'c', to_low, to_high, to, from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)
       
! n_links_max is typically 2 * number of cells in largest supercell
       allocate (g_adj(n_links_max, 2, g_lo%llim_proc:g_lo%ulim_alloc))

! now set up links_h:
! excluding wfb

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

! If this is not the first link in the chain, continue
          if (l_links(ik, it) > 0) then
! For each link to the right, do:
             iglo_star = iglo
             do j = 1, r_links(ik, it)
                call get_right_connection (iglo_star, iglo_right, ipright)
! sender
                if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
                if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
                iglo_star = iglo_right
             end do
          end if

          if (r_links(ik, it) > 0) then
! For each link to the left, do:
             iglo_star = iglo
             do j = 1, l_links(ik, it)
                call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
                if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
                if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
                iglo_star = iglo_left
             end do
          end if
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (nlambda > ng2 .and. il >= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          if (l_links(ik, it) > 0) then
! For each link to the right, do:
             iglo_star = iglo
             do j = 1, r_links(ik, it)
! get address of link
                call get_right_connection (iglo_star, iglo_right, ipright)
! sender
                if (ip == iproc) then
                   n = nn_from(ipright) + 1
                   nn_from(ipright) = n
                   from(ipright)%first(n) = ntgrid
                   from(ipright)%second(n) = 1
                   from(ipright)%third(n) = iglo
                end if
! receiver
                if (ipright == iproc) then
                   n = nn_to(ip) + 1
                   nn_to(ip) = n
                   to(ip)%first(n) = 2*l_links(ik, it) + j
                   to(ip)%second(n) = 1
                   to(ip)%third(n) = iglo_right
                end if
                iglo_star = iglo_right
             end do
          end if

          if (r_links(ik, it) > 0) then
! For each link to the left, do:
             iglo_star = iglo
             do j = 1, l_links(ik, it)
! get address of link
                call get_left_connection (iglo_star, iglo_left, ipleft)
! sender
                if (ip == iproc) then
                   n = nn_from(ipleft) + 1
                   nn_from(ipleft) = n
                   from(ipleft)%first(n) = -ntgrid
                   from(ipleft)%second(n) = 2   
                   from(ipleft)%third(n) = iglo
                end if
! receiver
                if (ipleft == iproc) then
                   n = nn_to(ip) + 1
                   nn_to(ip) = n
                   to(ip)%first(n) = 2*r_links(ik, it) + j
                   to(ip)%second(n) = 2
                   to(ip)%third(n) = iglo_left
                end if
                iglo_star = iglo_left
             end do
          end if
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       call init_fill (links_h, 'c', to_low, to_high, to, &
            from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)

! now take care of wfb (homogeneous part)

       nn_to = 0
       nn_from = 0 

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
! receiver
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do

! v_par < 0:
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
! receiver
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do
       
       do ip = 0, nproc-1
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do
       
       nn_from = 0
       nn_to = 0          

       do iglo = g_lo%llim_world, g_lo%ulim_world

          il = il_idx(g_lo,iglo)
          if (il /= ng2+1) cycle

          ip = proc_id(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          
          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle

          iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

! v_par > 0:
          call find_leftmost_link (iglo, iglo_right, ipright)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n
                from(ipright)%first(n) = ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if
! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = 2*ncell - r_links(ik, it)
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
             call get_right_connection (iglo_right, iglo_right, ipright)
          end do
 
! v_par < 0:
          call find_rightmost_link (iglo, iglo_left, ipleft)
          do j = 1, ncell
! sender
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n
                from(ipleft)%first(n) = -ntgrid
                from(ipleft)%second(n) = 2   
                from(ipleft)%third(n) = iglo
             end if
                
! receiver
             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to(ip)%first(n) = 2*ncell - l_links(ik, it)
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
             call get_left_connection (iglo_left, iglo_left, ipleft)
          end do
       end do

       from_low (1) = -ntgrid
       from_low (2) = 1
       from_low (3) = g_lo%llim_proc
       
       to_low (1) = 1
       to_low (2) = 1 
       to_low (3) = g_lo%llim_proc
       
       to_high(1) = n_links_max
       to_high(2) = 2
       to_high(3) = g_lo%ulim_alloc

       from_high(1) = ntgrid
       from_high(2) = 2
       from_high(3) = g_lo%ulim_alloc

       call init_fill (wfb_h, 'c', to_low, to_high, to, from_low, from_high, from)
       
       call delete_list (from)
       call delete_list (to)

200    continue

! Now set up class arrays for the implicit fields 
! i_class classes
! N_class(i) = number of linked cells for i_th class 
! M_class(i) = number of members in i_th class

! First count number of linked cells for each (kx, ky) 
       k = 1
       do it = 1, ntheta0
          do ik = 1, naky
             n_k(k) = 1 + l_links(ik, it) + r_links(ik, it)
             k = k + 1
          end do
       end do

! Count how many unique values of n_k there are.  This is the number 
! of classes.

! Sort: 
       do j = 1, naky*ntheta0-1
          do k = 1, naky*ntheta0-1                       
             if (n_k(k+1) < n_k(k)) then
                i = n_k(k)
                n_k(k) = n_k(k+1)
                n_k(k+1) = i
             end if
          end do
       end do

! Then count:
       i_class = 1
       do k = 1, naky*ntheta0-1
          if (n_k(k) == n_k(k+1)) cycle
          i_class = i_class + 1
       end do

! Allocate M, N:
       allocate (M_class(i_class))
       allocate (N_class(i_class))

! Initial values
       M_class = 1 ; N_class = 0

! Fill M, N arrays: 
       j = 1
       do k = 2, naky*ntheta0
          if (n_k(k) == n_k(k-1)) then
             M_class(j) = M_class(j) + 1
          else 
             N_class(j) = n_k(k-1)
             M_class(j) = M_class(j)/N_class(j)
             j = j + 1
          end if
       end do       
       j = i_class
       N_class(j) = n_k(naky*ntheta0)
       M_class(j) = M_class(j)/N_class(j)

! Check for consistency:
       
! j is number of linked cells in class structure
       j = 0
       do i = 1, i_class
          j = j + N_class(i)*M_class(i)
       end do

       if (j /= naky*ntheta0) then
          write(*,*) 'PE ',iproc,'has j= ',j,' k= ',naky*ntheta0,' : Stopping'
          stop
       end if

    end if

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

  subroutine allocate_arrays
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: g, gnew, gold
    use gs2_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    implicit none
    logical :: alloc = .true.

    if (alloc) then
       allocate (g    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
!       allocate (gold (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g0   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       if (nonlin) then
          allocate (gnl_1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gnl_2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          gnl_1 = 0. ; gnl_2 = 0.
       else
          allocate (gnl_1(1,2,1), gnl_2(1,2,1))
       end if
       if (boundary_option_switch == boundary_option_linked) then
          allocate (g_h(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          g_h = 0.
          allocate (save_h(2,g_lo%llim_proc:g_lo%ulim_alloc))
          save_h = .false.
       endif
    endif

    g = 0. ; gnew = 0. ; g0 = 0. !; gold = 0.

    alloc = .false.
  end subroutine allocate_arrays

  subroutine timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, mode)

    use theta_grid, only: ntgrid
    use collisions, only: solfp1
    use dist_fn_arrays, only: gnew, g, gold
    use nonlinear_terms, only: add_nonlinear_terms
    use hyper, only: hyper_diff
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, optional, intent (in) :: mode
    integer :: modep

    modep = 0
    if (present(mode)) modep = mode

    if (modep <= 0) then
       call add_nonlinear_terms (g, gnl_1, gnl_2, &
            phi, apar, bpar, istep, bkdiff(1), fexp(1))
       call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
       call hyper_diff (gnew, g0, phinew, bparnew)
       call kill (gnew, g0, phinew, bparnew)
       call solfp1 (gnew, g0, phinew, aparnew, bparnew)
       
       if (def_parity) then
          if (even) then
             gnew(-ntgrid:-1, 1,:) = gnew( ntgrid: 1:-1,2,:)
             gnew( 1: ntgrid, 1,:) = gnew(-1:-ntgrid:-1,2,:)
          else
             gnew( 1: ntgrid, 1,:) = -gnew(-1:-ntgrid:-1,2,:)
             gnew(-ntgrid:-1, 1,:) = -gnew( ntgrid: 1:-1,2,:)
          end if
       end if
       
    end if

  end subroutine timeadv

  subroutine exb_shear (g0, phi, apar, bpar)

    use gs2_layouts, only: ik_idx, it_idx, g_lo, is_kx_local, idx_local, idx
    use file_utils, only: error_unit
    use theta_grid, only: ntgrid
    use kt_grids, only: akx, aky, naky, ikx, ntheta0, box
    use le_grids, only: negrid, nlambda
    use species, only: nspec
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: kx_shift
    use gs2_time, only: code_dt

    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    integer, dimension(:), allocatable, save :: jump, ikx_indexed
    integer, dimension(1), save :: itmin
    integer :: ik, it, ie, is, il, ig, isgn, to_iglo, from_iglo, ierr
    real :: dkx, gdt
    logical :: exb_first = .true.
    logical :: kx_local

! If not in box configuration, return
    if (.not. box) then
       ierr = error_unit()
       g_exb = 0.
       write (ierr, fmt="('g_ExB set to zero unless box configuration is used.')")
       return
    end if

! If kx data is not local, no ExB shear will be included.

    call is_kx_local (negrid, nspec, nlambda, naky, ntheta0, kx_local)
    if (.not. kx_local) then
       ierr = error_unit()
       g_exb = 0.
       write (ierr, fmt="('Non-zero g_ExB not implemented for this layout.  g_ExB set to zero.')")
       return
    end if

! Initialize kx_shift, jump, idx_indexed

    if (exb_first) then
       exb_first = .false.
       allocate (kx_shift(naky), jump(naky))
       kx_shift = 0.
       jump = 0

       allocate (ikx_indexed(ntheta0))
       itmin = minloc (ikx)

       do it=itmin(1), ntheta0
          ikx_indexed (it+1-itmin(1)) = it
       end do

       do it=1,itmin(1)-1
          ikx_indexed (ntheta0 - itmin(1) + 1 + it)= it
       end do
    end if

    dkx = akx(2)
    
! BD: To do: Put the right timestep in here.
!
! For now, approximate Greg's dt == 1/2 (t_(n+1) - t_(n-1))
! with code_dt.  
!
! Note: at first time step, there is a difference of a factor of 2.
!

    gdt = code_dt

! kx_shift is a function of time.   Update it here:

    do ik=1, naky
       kx_shift(ik) = kx_shift(ik) - aky(ik)*g_exb*gdt
       jump(ik) = nint(kx_shift(ik)/dkx)
       kx_shift(ik) = kx_shift(ik) - jump(ik)*dkx
    end do

! BD: To do: Should save kx_shift array in restart file

    do ik = naky, 2, -1
       if (jump(ik) == 0) exit

       if (jump(ik) < 0) then

          if (fphi > epsilon(0.0)) then
             do it = 1, ntheta0 + jump(ik)
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = phi (ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = ntheta0 + jump(ik) + 1, ntheta0
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (fapar > epsilon(0.0)) then
             do it = 1, ntheta0 + jump(ik)
                do ig=-ntgrid,ntgrid
                   apar(ig,ikx_indexed(it),ik) = apar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = ntheta0 + jump(ik) + 1, ntheta0
                do ig=-ntgrid,ntgrid
                   apar (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (fbpar > epsilon(0.0)) then
             do it = 1, ntheta0 + jump(ik)
                do ig=-ntgrid,ntgrid
                   bpar(ig,ikx_indexed(it),ik) = bpar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = ntheta0 + jump(ik) + 1, ntheta0
                do ig=-ntgrid,ntgrid
                   bpar (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For some data layouts, there is no communication required.  Try such a case first.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (kx_local) then
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda
                      do it = 1, ntheta0 + jump(ik)
                         
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)
                         
                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = g0(ig,isgn,from_iglo)
                               end do
                            end do
                         end if
                      end do

                      do it = ntheta0 + jump(ik) + 1, ntheta0
                         
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         
                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = 0.
                               end do
                            end do
                         end if
                      end do
                   end do
                end do
             end do
          end if

       else  ! case for jump(ik) > 0

          if (fphi > epsilon(0.0)) then
             do it = ntheta0, 1+jump(ik), -1
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = phi (ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = jump(ik), 1, -1
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (fapar > epsilon(0.0)) then
             do it = ntheta0, 1+jump(ik), -1
                do ig=-ntgrid,ntgrid
                   apar(ig,ikx_indexed(it),ik) = apar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = jump(ik), 1, -1
                do ig=-ntgrid,ntgrid
                   apar(ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (fbpar > epsilon(0.0)) then
             do it = ntheta0, 1+jump(ik), -1
                do ig=-ntgrid,ntgrid
                   bpar(ig,ikx_indexed(it),ik) = bpar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = jump(ik), 1, -1
                do ig=-ntgrid,ntgrid
                   bpar(ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For some data layouts, there is no communication required.  Try such a case first.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (kx_local) then
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda
                      do it = ntheta0, 1+jump(ik), -1

                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)

                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = g0(ig,isgn,from_iglo)
                               end do
                            end do
                         end if
                      end do

                      do it = jump(ik), 1, -1

                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)

                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = 0.
                               end do
                            end do
                         end if
                      end do
                   end do
                end do
             end do
          end if

       end if

    end do
       
  end subroutine exb_shear

  subroutine kill (g0, g1, phi, bpar)
    
    use gs2_layouts, only: ik_idx, it_idx
    use kt_grids, only: akx, aky
    use theta_grid, only: ntgrid, bmag
    use run_parameters, only: fphi, fbpar, tnorm
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: kperp2
    use le_grids, only: integrate, geint2g, lintegrate, nlambda, al, gint2g
    use gs2_layouts, only: g_lo, geint_lo, gint_lo
    use constants
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0, g1
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    complex, dimension (:,:), allocatable :: g0eint, g1eint
    
    real, dimension (:,:), allocatable, save :: aintnorm
    real :: x
    integer :: iglo, ik, it, il, ige, ig, igint
    logical :: diff_first = .true.
    real :: chi_int ! chi_int needs to be calculated before it can be used.

    if (D_kill < 0.) return

    if (noise > 0.) then
       chi_int = 1.  ! placeholder
       D_kill = 0.5*sqrt(pi/2.)*noise*chi_int
    end if

    if (save_n) then
       allocate (g0eint(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
       allocate (g1eint(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))

       if (diff_first) then
          diff_first = .false.
          allocate (aintnorm(-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
          aintnorm = 0. ; g1 = 1.0
          call integrate (g1, g1eint)
          do ige = geint_lo%llim_proc, geint_lo%ulim_proc
             aintnorm(:,ige) = 1.0/real(g1eint(:,ige))
          end do
       end if
! get initial particle number
       call integrate (g0, g0eint)
    end if

    if (h_kill) call g_adjust (g0, phi, bpar, fphi, fbpar)

! kill higher moments of f
    if (kill_grid) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          g1(:,1,iglo) = g0(:,1,iglo) &
               * exp(-0.5*(akx(it)**2 + aky(ik)**2) * D_kill * code_dt / tnorm)
          g1(:,2,iglo) = g0(:,2,iglo) &
               * exp(-0.5*(akx(it)**2 + aky(ik)**2) * D_kill * code_dt / tnorm)
       end do
    else       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          g1(:,1,iglo) = g0(:,1,iglo) &
               * exp(-0.5*kperp2(:,it,ik) * D_kill * code_dt / tnorm)
          g1(:,2,iglo) = g0(:,2,iglo) &
               * exp(-0.5*kperp2(:,it,ik) * D_kill * code_dt / tnorm)
       end do
    end if

    if (h_kill) call g_adjust (g1, phi, bpar, -fphi, -fbpar)

! restore particle number
    if (save_n) then
       call integrate (g1, g1eint)
       do ige = geint_lo%llim_proc, geint_lo%ulim_proc
          g0eint(:,ige) = (g0eint(:,ige) - g1eint(:,ige)) &
               *aintnorm(:,ige)
       end do
       call geint2g (g0eint, g0)
       g0 = g0 + g1
    else
       g0 = g1
    endif

    if (save_n) deallocate (g0eint, g1eint)

  end subroutine kill

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

  subroutine get_source_term &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        isgn, iglo, sourcefac, source)
    use dist_fn_arrays, only: aj0, aj1, vperp2, vpar, vpac, g, ittp
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: aky, theta0, akx
    use le_grids, only: nlambda, ng2, lmax, anon, e, negrid
    use species, only: spec, nspec
    use run_parameters, only: fphi, fapar, fbpar, wunits, tunits
    use gs2_time, only: code_dt
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use run_parameters, only: tnorm
    use nonlinear_terms, only: nonlin
    use hyper, only: D_res
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: isgn, iglo
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:), intent (out) :: source
    real :: tfac, timep

    integer :: ig, ik, it, il, ie, is
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg!, bpargavg !GGH added bpargavg

!    call timer (0, 'get_source_term')

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    phigavg  = (fexp(is)*phi(:,it,ik)   + (1.0-fexp(is))*phinew(:,it,ik)) &
                *aj0(:,iglo)*fphi &
             + (fexp(is)*bpar(:,it,ik) + (1.0-fexp(is))*bparnew(:,it,ik))&
                *aj1(:,iglo)*fbpar*2.0*vperp2(:,iglo)*spec(is)%tz
    apargavg = (fexp(is)*apar(:,it,ik)  + (1.0-fexp(is))*aparnew(:,it,ik)) &
                *aj0(:,iglo)*fapar

! source term in finite difference equations
    select case (source_option_switch)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Default choice: solve self-consistent equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms + include external i omega_d * F_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_phiext_full)
       if (il <= lmax) then
          call set_source
          if (istep > 0 .and. aky(ik) < epsilon(0.0)) then
             source(:ntgrid-1) = source(:ntgrid-1) &
                  - zi*anon(ie,is)*wdrift(:ntgrid-1,iglo)*2.0*phi_ext*sourcefac &
                  *(spec(is)%fprim+spec(is)%tprim*(e(ie,is)-1.5))
          end if
       else
          source = 0.0
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms + include external i omega_d Phi * F_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_hm_force)
       if (il <= lmax) then
          call set_source
          if (istep > 0 .and. aky(int(aky_star)) == aky(ik) .and. akx(int(akx_star)) == akx(it)) &
               source(:ntgrid-1) = source(:ntgrid-1) - zi*phi_ext*sourcefac
       else 
          source = 0.0
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms + include external <Phi> * F_0 * other stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_test2_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if
       if (aky(ik) == 0.0 .and. istep > 0) then
          source(:ntgrid-1) = source(:ntgrid-1) &
               + tunits(ik)*code_dt &
               *(aj0(:ntgrid-1,iglo)**2 + aj0(-ntgrid+1:,iglo)**2) &
               *(e(ie,is)-1.5)*sourcefac &
               *exp(-(theta(:ntgrid-1)/thetas)**2)
       end if
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms for ky=0, something else for ky /= 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_convect_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if
       if (aky(ik) /= 0.0 .and. istep > 0) then
          source(:ntgrid-1) = sourcefac &
               *exp(cmplx((theta(:ntgrid-1)-theta0(it,ik))**2, &
               k0*theta(:ntgrid-1)))
       end if       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Include no source term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_zero)
       source = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! S = sin(theta)*f(omega)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_sine)
       source(:ntgrid-1) = tunits(ik)*code_dt &
            *(sin(theta(:ntgrid-1)) + sin(theta(-ntgrid+1:))) &
            *sourcefac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! S = -i*omega_d*f(omega)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_test1)
       do ig = -ntgrid, ntgrid-1
          source(ig) = -zi*(wdrift_func(ig,il,ie,it,ik,is) &
               + wdrift_func(ig+1,il,ie,it,ik,is)) &
               *sourcefac
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! S = cos(theta)*f(omega)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_cosine)
       source(:ntgrid-1) = tunits(ik)*code_dt &
            *(cos(theta(:ntgrid-1)) + cos(theta(-ntgrid+1:))) &
            *sourcefac
    end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do matrix multiplications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! special source term for totally trapped particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (source_option_switch == source_option_full .or. &
        source_option_switch == source_option_phiext_full .or. &
        source_option_switch == source_option_hm_force .or. &
        source_option_switch == source_option_test2_full .or. &
        source_option_switch == source_option_test1) then
       if (nlambda > ng2 .and. isgn == 2) then
          do ig = -ntgrid, ntgrid
             if (il /= ittp(ig)) cycle
             source(ig) &
                  = g(ig,2,iglo)*a(ig,iglo) &
                  - anon(ie,is)*zi*wdriftttp(ig,it,ik,ie,is)*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*phigavg(ig)
          end do

          if (source_option_switch == source_option_phiext_full .and.  &
               aky(ik) < epsilon(0.0) .and. istep > 0) then
             do ig = -ntgrid, ntgrid
                if (il /= ittp(ig)) cycle             
                source(ig) = source(ig) - zi*anon(ie,is)* &
                     wdriftttp(ig,it,ik,ie,is)*phi_ext*sourcefac &
                     *(spec(is)%fprim+spec(is)%tprim*(e(ie,is)-1.5))          
             end do
          endif

! add in nonlinear terms -- tfac normalizes the *amplitudes*.
          if (nonlin) then         
             tfac = 1./tnorm
             select case (istep)
             case (0)
                ! nothing
             case (1)
                do ig = -ntgrid, ntgrid
                   if (il /= ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*tfac*gnl_1(ig,isgn,iglo)
                end do
             case default
                do ig = -ntgrid, ntgrid
                   if (il /= ittp(ig)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*tfac*( &
                        1.5*gnl_1(ig,isgn,iglo) - 0.5*gnl_2(ig,isgn,iglo))
                end do
             end select
          end if

       end if
    end if

!    call timer (1, 'get_source_term')
  contains

    subroutine set_source

      complex :: apar_p, apar_m, phi_p, phi_m!, bpar_p !GGH added bpar_p
      real, dimension(:,:), allocatable, save :: ufac
      real :: bd, bdfac_p, bdfac_m
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

! try fixing bkdiff dependence
      bd = bkdiff(1)

      bdfac_p = 1.+bd*(3.-2.*real(isgn))
      bdfac_m = 1.-bd*(3.-2.*real(isgn))

      do ig = -ntgrid, ntgrid-1
         phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
         phi_m = phigavg(ig+1)-phigavg(ig)
         apar_p = apargavg(ig+1)+apargavg(ig)
         apar_m = aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
              -apar(ig+1,it,ik)-apar(ig,it,ik)
 
         source(ig) = anon(ie,is)*(-2.0*vpar(ig,isgn,iglo)*phi_m &
              -spec(is)%zstm*vpac(ig,isgn,iglo) &
              *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
              + D_res(it,ik)*apar_p) &
              -zi*wdrift(ig,iglo)*phi_p) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is)) &
              *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo)) 
      end do
        
! add in nonlinear terms -- tfac normalizes the *amplitudes*.
      if (nonlin) then         
         tfac = 1./tnorm
         select case (istep)
         case (0)
            ! nothing
         case (1)
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*tfac*gnl_1(ig,isgn,iglo)
            end do
         case default
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*tfac*( &
                    1.5*gnl_1(ig,isgn,iglo) - 0.5*gnl_2(ig,isgn,iglo))
            end do
         end select
      end if

    end subroutine set_source

  end subroutine get_source_term

  subroutine invert_rhs_1 &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        iglo, sourcefac)
    use dist_fn_arrays, only: gnew, ittp
    use run_parameters, only: eqzip, secondary, tertiary
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, ng2, lmax, forbid
    use kt_grids, only: aky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac

    integer :: ig, ik, it, il, ie, isgn, is
    integer :: ilmin
    complex :: beta1
    complex, dimension (-ntgrid:ntgrid,2) :: source, g1, g2
    logical :: kperiod_flag, speriod_flag
    integer :: ntgl, ntgr

    call prof_entering ("invert_rhs_1", "dist_fn")

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    if (eqzip) then
       if (secondary .and. ik == 2 .and. it == 1) return ! do not evolve primary mode
       if (tertiary .and. ik == 1) then
          if (it == 2 .or. it == ntheta0) return ! do not evolve periodic equilibrium
       end if
    end if

    do isgn = 1, 2
       call get_source_term (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, isgn, iglo, sourcefac, source(:,isgn))
    end do

    ! gnew is the inhomogeneous solution
    gnew(:,:,iglo) = 0.0

    ! g1 is the homogeneous solution
    g1 = 0.0

    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       kperiod_flag = .true.
    case (boundary_option_linked)
       kperiod_flag = .true.
       speriod_flag = aky(ik) == 0.0
    case default
       kperiod_flag = .false.
    end select

    kperiod_flag = kperiod_flag .or. aky(ik) == 0.0

    ntgl = -ntgrid
    ntgr =  ntgrid

! ng2+1 is WFB

    if (kperiod_flag) then
       if (il <= ng2+1) then
          g1(ntgl,1) = 1.0
          g1(ntgr,2) = 1.0
       end if
    end if
!!!
    if (il == ng2+1) then
       g1(ntgl,1) = wfb  ! wfb should be unity here; variable is for testing
       g1(ntgr,2) = wfb  ! wfb should be unity here; variable is for testing
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
!!!       ilmin = ng2 + 2
       ilmin = ng2 + 1
    end if

    ! time advance vpar < 0 homogeneous part
    if (il >= ilmin) then
       do ig = ntgr-1, ntgl, -1
          g1(ig,2) = -g1(ig+1,2)*r(ig,iglo) + g2(ig,2)
       end do
    end if

    if (nlambda > ng2 .and. il >= ng2+2 .and. il <= lmax) then
       ! match boundary conditions at lower bounce point
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

! case of linked .and. WFB has problem, particularly when nu > 0, dt << 1
! wfb solution below fails with real links.  need to balance this properly
! this has been fixed I believe (should double check)

    if (boundary_option_switch == boundary_option_linked) then
       if (speriod_flag .and. il <= ng2+1) then
          call self_periodic
       else 
       ! save homogeneous solution as necessary
          if (save_h (1, iglo)) g_h(:,1,iglo) = g1(:,1)
          if (save_h (2, iglo)) g_h(:,2,iglo) = g1(:,2)
       end if

! wfb (isolated)
       if (il == ng2+1 .and. .not. connections(iglo)%neighbor) &
            call self_periodic

    else       
       ! add correct amount of homogeneous solution now
       if (kperiod_flag .and. il <= ng2+1) then
          call self_periodic
!!!
       else if (il == ng2 + 1) then
          call self_periodic
       end if

    end if

    ! add correct amount of homogeneous solution for trapped particles
    if (il >= ng2+2 .and. il <= lmax) then
       beta1 = 0.0
       do ig = ntgr-1, ntgl, -1
          if (ittp(ig) == il) cycle
          if (forbid(ig,il)) then
             beta1 = 0.0
          else if (forbid(ig+1,il)) then
             beta1 = (gnew(ig,1,iglo) - gnew(ig,2,iglo))/(1.0 - g1(ig,1))
          end if
          gnew(ig,1,iglo) = gnew(ig,1,iglo) + beta1*g1(ig,1)
          gnew(ig,2,iglo) = gnew(ig,2,iglo) + beta1*g1(ig,2)
       end do
    end if

    if (def_parity) then
       if (even) then
          gnew(ntgl:-1,1,iglo) = gnew( ntgr:1:-1,2,iglo)
          gnew(1:ntgr, 1,iglo) = gnew(-1:ntgl:-1,2,iglo)
       else
          gnew(1:ntgr, 1,iglo) = -gnew(-1:ntgl:-1,2,iglo)
          gnew(ntgl:-1,1,iglo) = -gnew( ntgr:1:-1,2,iglo)
       end if
    end if

    ! zero out spurious gnew outside trapped boundary
    where (forbid(:,il))
       gnew(:,1,iglo) = 0.0
       gnew(:,2,iglo) = 0.0
    end where

!    if (istep == 1) then
!       write(*,*) kperiod_flag
!       if (ik == 2 .and. it == 1) then
!          do ig = -ntgrid, ntgrid
!             write(*,fmt="(5(1x,e10.4),4(1x,i5))") theta(ig), &
!                  gnew(ig,1,iglo), gnew(ig,2,iglo), il, ie, ik, it
!          end do
!          write(*,*) 
!       end if
!    end if

    call prof_leaving ("invert_rhs_1", "dist_fn")

  contains

    subroutine self_periodic

      if (g1(ntgr,1) /= 1.) then
         beta1 = (gnew(ntgr,1,iglo) - gnew(ntgl,1,iglo))/(1.0 - g1(ntgr,1))
         gnew(:,1,iglo) = gnew(:,1,iglo) + beta1*g1(:,1)
      end if

      if (g1(ntgl,2) /= 1.) then
         beta1 = (gnew(ntgl,2,iglo) - gnew(ntgr,2,iglo))/(1.0 - g1(ntgl,2))
         gnew(:,2,iglo) = gnew(:,2,iglo) + beta1*g1(:,2)
      end if
      
    end subroutine self_periodic

  end subroutine invert_rhs_1

  subroutine invert_rhs_linked &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, ng2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx
    use redistribute, only: fill
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    complex, intent (in) :: sourcefac

    complex :: b0, fac, facd, b1
    integer :: il, ik, it, n, i, j, it_star
    integer :: iglo, ncell

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       call invert_rhs_1 (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, iglo, sourcefac)
    end do

!!    if (no_comm .or. istep == 0) then
    if (no_comm) then
       ! nothing
    else       
       call fill (links_p, gnew, g_adj)
       call fill (links_h, g_h, g_adj)
       call fill (wfb_p, gnew, g_adj)
       call fill (wfb_h, g_h, g_adj)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)

          ncell = r_links(ik, it) + l_links(ik, it) + 1
          if (ncell == 1) cycle
! wfb
          if (nlambda > ng2 .and. il == ng2+1) then
             if (save_h(1, iglo)) then

                facd = 1.0
                do j = 1, ncell
                   facd = facd * g_adj(ncell+j,1,iglo)
                end do
                facd = 1./(1.-facd)
                
                b0 = 0.
                do i = 1, ncell-1
                   fac = 1.0
                   do j = i+1, ncell
                      fac = fac * g_adj(ncell+j,1,iglo)
                   end do
                   b0 = b0 + fac * g_adj(ncell+1-i,1,iglo)
                end do
                b0 = (b0 + g_adj(1,1,iglo))*facd

                do i = 1, l_links(ik, it)
                   b0 = b0 * g_adj(ncell+i,1,iglo) + g_adj(ncell+1-i,1,iglo)
                end do
                
                gnew(:,1,iglo) = gnew(:,1,iglo) + b0*g_h(:,1,iglo)
             endif

             if (save_h(2, iglo)) then

                facd = 1.0
                do j = 1, ncell
                   facd = facd * g_adj(ncell+j,2,iglo)
                end do
                facd = 1./(1.-facd)
                
                b0 = 0.
                do i = 1, ncell-1
                   fac = 1.0
                   do j = i+1, ncell
                      fac = fac * g_adj(ncell+j,2,iglo)
                   end do
                   b0 = b0 + fac * g_adj(ncell+1-i,2,iglo)
                end do
                b0 = (b0 + g_adj(1,2,iglo))*facd

                do i = 1, r_links(ik, it)
                   b0 = b0 * g_adj(ncell+i,2,iglo) + g_adj(ncell+1-i,2,iglo)
                end do
                
                gnew(:,2,iglo) = gnew(:,2,iglo) + b0*g_h(:,2,iglo)
             end if
          else
!
! n_links is the number of complex numbers required to fix the boundary 
! conditions in each cell that is a member of a supercell with at least two
! cells, and for which the bounce point is not at theta=pi.
!
             if (save_h(1, iglo)) then
                n = n_links(1, ik, it)
                b0 = 0.0
                do i = 1, l_links(ik, it)
                   fac = 1.0
                   do j = 1, i-1
                      fac = fac * g_adj(n+1-j, 1, iglo)
                   end do
                   b0 = b0 + g_adj(i,1,iglo) * fac
                end do
                
                gnew(:,1,iglo) = gnew(:,1,iglo) + b0*g_h(:,1,iglo)
             end if

             if (save_h(2, iglo)) then
                n = n_links(2, ik, it)
                b0 = 0.0
                do i = 1, r_links(ik, it)
                   fac = 1.0
                   do j = 1, i-1
                      fac = fac * g_adj(n+1-j, 2, iglo)
                   end do
                   b0 = b0 + g_adj(i,2,iglo) * fac
                end do
                
                gnew(:,2,iglo) = gnew(:,2,iglo) + b0*g_h(:,2,iglo)
             end if
          end if

       end do
       
    end if

  end subroutine invert_rhs_linked

  subroutine invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use gs2_time, only: code_time
    use constants
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep

    integer :: iglo

    real :: time, timep
    complex :: sourcefac

    call prof_entering ("invert_rhs", "dist_fn")

    time = code_time
    if (time > t0) then
       sourcefac = source0*exp(-zi*omega0*time+gamma0*time)
    else
       sourcefac = (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time+gamma0*time)
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)
       call invert_rhs_linked &
            (phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac) 
    case default
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          call invert_rhs_1 (phi, apar, bpar, phinew, aparnew, bparnew, &
               istep, iglo, sourcefac)
       end do
    end select

    call prof_leaving ("invert_rhs", "dist_fn")
  end subroutine invert_rhs

  subroutine getan (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew
    use dist_fn_arrays, only: kperp2
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species
    use run_parameters, only: beta, fphi, fapar, fbpar
    use prof, only: prof_entering, prof_leaving
    use gs2_layouts, only: g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt

    integer :: isgn, iglo, ig

    call prof_entering ("getan", "dist_fn")

    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do

       wgt = spec%z*spec%dens
       call integrate_species (g0, wgt, antot)

!    if (kfilter > epsilon(0.0)) call par_filter(antot)
       if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)

    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(ig,iglo)*vpa(ig,isgn,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       
       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g0, wgt, antota)
!    if (kfilter > epsilon(0.0)) call par_filter(antota)

    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj1(ig,iglo)*vperp2(ig,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       wgt = spec%temp*spec%dens
       call integrate_species (g0, wgt, antotp)
!    if (kfilter > epsilon(0.0)) call par_filter(antotp)

    end if
    call prof_leaving ("getan", "dist_fn")
  end subroutine getan

  subroutine getmoms (phi, ntot, density, upar, tpar, tperp)
    use dist_fn_arrays, only: vpa, vperp2, aj0, gnew
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, anon
    use prof, only: prof_entering, prof_leaving
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density, &
         upar, tpar, tperp, ntot

    integer :: ik, it, isgn, ie, is, iglo, ig

! returns moment integrals to PE 0
    call prof_entering ("getmoms", "dist_fn")

! total density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)

       do isgn = 1, 2
          g0(:,isgn,iglo) = (aj0(:,iglo)**2-1.0)*anon(ie,is) &
               *phi(:,it,ik)*spec(is)%zt*spec(is)%dens
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo) + g0(ig,isgn,iglo)
          end do
       end do
    end do
    call integrate_moment (g0, ntot)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
!             g0(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
             g0(ig,isgn,iglo) = gnew(ig,isgn,iglo)
          end do
       end do
    end do

! guiding center density
    call integrate_moment (gnew, density)

! guiding center upar
    g0 = vpa*gnew

    call integrate_moment (g0, upar)

! guiding center tpar
    g0 = 2.*vpa*g0

    call integrate_moment (g0, tpar)
    tpar = tpar - density

! guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
!             g0(ig,isgn,iglo) = vperp2(ig,iglo)*gnew(ig,isgn,iglo)*aj0(ig,iglo)
             g0(ig,isgn,iglo) = vperp2(ig,iglo)*gnew(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (g0, tperp)
    tperp = tperp - density

    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       density(:,:,:,is)=density(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
    end do

    call prof_leaving ("getmoms", "dist_fn")
  end subroutine getmoms

  subroutine init_fieldeq
    use dist_fn_arrays, only: aj0, aj1, vperp2, kperp2
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky
    use le_grids, only: anon, integrate_species
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use run_parameters, only: tite
    implicit none
    integer :: iglo, isgn
    integer :: ik, it, ie, is
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: tot
    real, dimension (nspec) :: wgt
    logical :: done = .false.

    if (done) return
    done = .true.

    allocate (gridfac1(-ntgrid:ntgrid,ntheta0,naky))
    gridfac1 = 1.0
    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       ! nothing
    case (boundary_option_linked)
       do it = 1, ntheta0
          do ik = 1, naky
             if (aky(ik) == 0.0) cycle
             if (itleft(ik,it) < 0) gridfac1(-ntgrid,it,ik) = gridfac
             if (itright(ik,it) < 0) gridfac1(ntgrid,it,ik) = gridfac
          end do
       end do
    case default
       do ik = 1, naky
          if (aky(ik) == 0.0) cycle
          gridfac1(-ntgrid,:,ik) = gridfac
          gridfac1(ntgrid,:,ik) = gridfac
       end do
    end select

    allocate (gamtot(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot1(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot2(-ntgrid:ntgrid,ntheta0,naky))
    if (adiabatic_option_switch == adiabatic_option_fieldlineavg .or. &	
         adiabatic_option_switch == adiabatic_option_noJ) then	
       allocate (gamtot3(-ntgrid:ntgrid,ntheta0,naky))
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
             if (aky(ik) > epsilon(0.0)) gamtot(:,:,ik) = gamtot(:,:,ik) + tite
          end do
       elseif (adiabatic_option_switch == adiabatic_option_fieldlineavg .or. &
            adiabatic_option_switch == adiabatic_option_noJ) then
          gamtot  = gamtot + tite
          gamtot3 = (gamtot-tite) / gamtot
          where (gamtot3 < 2.*epsilon(0.0)) gamtot3 = 1.0
       else
          gamtot = gamtot + tite 
       endif
    endif
  end subroutine init_fieldeq

  subroutine getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use dist_fn_arrays, only: kperp2
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    real, allocatable, dimension(:,:), save :: fl_avg, awgt
    integer :: ik, it
    logical :: first = .true.
    
    if (first) allocate (fl_avg(ntheta0, naky))
    fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_noJ) then
          
          if (first) then 
             allocate (awgt(ntheta0, naky))
             awgt = 0.
             do ik = 1, naky
                do it = 1, ntheta0
                   if (aky(ik) > epsilon(0.0)) cycle
                   awgt(it,ik) = 1.0/sum(delthet*gamtot3(:,it,ik))
                end do
             end do
          endif
          
          do ik = 1, naky
             do it = 1, ntheta0
                fl_avg(it,ik) = tite*sum(delthet*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
             end do
          end do
          
       end if

       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          
          if (first) then 
             allocate (awgt(ntheta0, naky))
             awgt = 0.
             do ik = 1, naky
                do it = 1, ntheta0
                   if (aky(ik) > epsilon(0.0)) cycle
                   awgt(it,ik) = 1.0/sum(delthet*jacob*gamtot3(:,it,ik))
                end do
             end do
          endif
           
          do ik = 1, naky
             do it = 1, ntheta0
                fl_avg(it,ik) = tite*sum(delthet*jacob*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
             end do
          end do

!          do it=2,ntheta0
!             fl_avg(it, 1) = fl_avg(it, 1) - 1.0
!          end do
       end if

    end if

    if (fphi > epsilon(0.0)) then
       fieldeq = antot + bpar*gamtot1 - gamtot*gridfac1*phi 

       if (.not. has_electron_species(spec)) then
          do ik = 1, naky
             do it = 1, ntheta0
                fieldeq(:,it,ik) = fieldeq(:,it,ik) + fl_avg(it,ik)
             end do
          end do
       end if
    end if

    if (fapar > epsilon(0.0)) then
       fieldeqa = antota - kperp2*gridfac1*apar
    end if

    if (fbpar > epsilon(0.0)) then
       fieldeqp = (antotp+bpar*gamtot2+0.5*phi*gamtot1)*beta*apfac
       do ik = 1, naky
          do it = 1, ntheta0
             fieldeqp(:,it,ik) = fieldeqp(:,it,ik)/bmag(:)**2
          end do
       end do
       fieldeqp = fieldeqp + bpar*gridfac1
    end if

    first = .false.

  end subroutine getfieldeq1

  subroutine getfieldeq (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (:,:,:), allocatable :: antot, antota, antotp

    allocate (antot (-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))

    call getan (antot, antota, antotp)
    call getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)

    deallocate (antot, antota, antotp)
  end subroutine getfieldeq
  
  subroutine getfieldexp (phi, apar, bpar)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
    complex, dimension (:,:,:), allocatable :: antot, antota, antotp

    allocate (antot (-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))

    call getan (antot, antota, antotp)
    call getfieldeq2 (phi, apar, bpar, antot, antota, antotp)

    deallocate (antot, antota, antotp)
  end subroutine getfieldexp
  
  subroutine getfieldeq2 (phi, apar, bpar, antot, antota, antotp)
    use dist_fn_arrays, only: kperp2
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    real, allocatable, dimension(:,:,:), save :: f1, f2, f3, f4, kp2
    real, allocatable, dimension(:,:), save :: fl_avg, awgt
    real, allocatable, dimension(:), save :: bfac
    real :: f5
    integer :: ig, ik, it
    logical :: first = .true.
    
    
    if (first) then
! prepare for field-line-average term: 
       allocate (fl_avg(ntheta0, naky))
       fl_avg = 0.
! prepare for field solves: 
       allocate (bfac(-ntgrid:ntgrid))
       allocate (f1(-ntgrid:ntgrid,ntheta0,naky))
       allocate (f2(-ntgrid:ntgrid,ntheta0,naky))
       allocate (f3(-ntgrid:ntgrid,ntheta0,naky))
       allocate (f4(-ntgrid:ntgrid,ntheta0,naky))
       allocate (kp2(-ntgrid:ntgrid,ntheta0,naky))
       bfac = beta*apfac/bmag**2
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                f5 = 1. + 0.5*bfac(ig)*gamtot1(ig,it,ik)**2 &
                     / (1.+bfac(ig)*gamtot2(ig,it,ik))
                f1(ig,it,ik) = 1.0 / gamtot(ig,it,ik) / f5
                f3(ig,it,ik) = -bfac(ig) / (1.0 + bfac(ig)*gamtot2(ig,it,ik))
                f2(ig,it,ik) = gamtot1(ig,it,ik)*f3(ig,it,ik) / f5
                f4(ig,it,ik) = gamtot1(ig,it,ik)*f3(ig,it,ik) * 0.5
             end do
          end do
       end do

! Avoid dividing by zero for the kx=0, ky=0 mode:
       kp2 = kperp2
       where (kp2 == 0.) 
          kp2 = 1.0
       end where

    endif

    if (.not. has_electron_species(spec)) then
       
       if (adiabatic_option_switch == adiabatic_option_noJ) then

          if (first) then 
             allocate (awgt(ntheta0, naky))
             awgt = 0.
             do ik = 1, naky
                do it = 1, ntheta0
                   if (aky(ik) > epsilon(0.0)) cycle
                   awgt(it,ik) = 1.0/sum(delthet*gamtot3(:,it,ik))
                end do
             end do
          endif
          
          do ik = 1, naky
             do it = 1, ntheta0
                fl_avg(it,ik) = tite*sum(delthet*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
             end do
          end do
       end if

       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then

          if (first) then 
             allocate (awgt(ntheta0, naky))
             awgt = 0.
             do ik = 1, naky
                do it = 1, ntheta0
                   if (aky(ik) > epsilon(0.0)) cycle
                   awgt(it,ik) = 1.0/sum(delthet*jacob*gamtot3(:,it,ik))
                end do
             end do
          endif
          
          do ik = 1, naky
             do it = 1, ntheta0
                fl_avg(it,ik) = tite*sum(delthet*jacob*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
             end do
          end do
       end if

    end if

! main loop: 
    do ik = 1, naky
       do it = 1, ntheta0
          phi(:,it,ik) = fphi*(antot(:,it,ik)*f1(:,it,ik) + fl_avg(it,ik)*f1(:,it,ik) &
               + antotp(:,it,ik)*f2(:,it,ik))
          bpar(:,it,ik) = fbpar*antotp(:,it,ik)*f3(:,it,ik) + phi(:,it,ik)*f4(:,it,ik)
          apar(:,it,ik) = fapar*antota(:,it,ik)/kp2(:,it,ik)
       end do
    end do

    first = .false.

  end subroutine getfieldeq2

  elemental function j0 (x)
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

  elemental function j1 (x)
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

  subroutine lambda_flux (phi, lamflux)

    use le_grids, only: pintegrate, e
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0, g, vpa, vperp2
    use run_parameters, only: fphi
    use gs2_layouts, only: is_idx, ie_idx, g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:), intent (out) :: lamflux
    real :: etmp
    integer :: isgn, iglo

    if (fphi > epsilon(0.)) then
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = g(:,isgn,iglo)*aj0(:,iglo)
          end do
       end do
       call pintegrate (g0, phi, lamflux(:,:,1))
       
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          etmp = e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
          g0(:,1,iglo) = g0(:,1,iglo)*etmp
          g0(:,2,iglo) = g0(:,2,iglo)*etmp
       end do
       call pintegrate (g0, phi, lamflux(:,:,2))
       
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = g(:,isgn,iglo)*2.*vpa(:,isgn,iglo)**2*aj0(:,iglo)
          end do
       end do
       call pintegrate (g0, phi, lamflux(:,:,3))
       
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = g(:,isgn,iglo)*vperp2(:,iglo)*aj0(:,iglo)
          end do
       end do
       call pintegrate (g0, phi, lamflux(:,:,4))
    else
       lamflux = 0.
    end if

  end subroutine lambda_flux
      
  subroutine e_flux (phi, enflux)

    use le_grids, only: pe_integrate, e
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0, g, vpa, vperp2
    use run_parameters, only: fphi
    use gs2_layouts, only: is_idx, ie_idx, g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:) :: enflux
    integer :: isgn, iglo

    if (fphi > epsilon(0.)) then
       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*aj0
       end do
       call pe_integrate (g0, phi, enflux(:,:,1))
       
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call pe_integrate (g0, phi, enflux(:,:,2))
       
       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*2.*vpa(:,isgn,:)**2*aj0
       end do
       call pe_integrate (g0, phi, enflux(:,:,3))
       
       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*vperp2*aj0
       end do
       call pe_integrate (g0, phi, enflux(:,:,4))
    else
       enflux = 0.
    end if

  end subroutine e_flux
      
!  subroutine flux (phi, apar, bpar, &
!       pflux, qflux, qflux_par, qflux_perp, vflux, &
!       pmflux, qmflux, qmflux_par, qmflux_perp, vmflux, &
!       pbflux, qbflux, qbflux_par, qbflux_perp, vbflux, anorm)
  subroutine flux (phi, apar, bpar, &
        pflux,  qflux,  vflux, &
       pmflux, qmflux, vmflux, &
       pbflux, qbflux, vbflux, &
       theta_pflux, theta_vflux, theta_qflux, &
       theta_pmflux, theta_vmflux, theta_qmflux, & 
       theta_pbflux, theta_vbflux, theta_qbflux, &
       anorm)
    use species, only: spec
    use theta_grid, only: ntgrid, bmag, gradpar, grho, delthet
    use kt_grids, only: naky, ntheta0
    use le_grids, only: e
    use dist_fn_arrays, only: g, aj0, vpac, vpa, aj1, vperp2
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use mp, only: proc0
    use run_parameters, only: woutunits, fphi, fapar, fbpar
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:,:), intent (out) :: pflux, pmflux, pbflux
    real, dimension (:,:,:), intent (out) :: vflux, vmflux, vbflux
    real, dimension (:,:,:,:), intent (out) :: qflux, qmflux, qbflux
    real, dimension (-ntgrid:,:), intent (out) :: theta_pflux, theta_pmflux, theta_pbflux
    real, dimension (-ntgrid:,:), intent (out) :: theta_vflux, theta_vmflux, theta_vbflux
    real, dimension (-ntgrid:,:,:), intent (out) :: theta_qflux, theta_qmflux, theta_qbflux
!    real, dimension (:,:,:), intent (out) :: qflux_par, qmflux_par, qbflux_par
!    real, dimension (:,:,:), intent (out) :: qflux_perp, qmflux_perp, qbflux_perp
    real, intent (out) :: anorm
    real, dimension (:,:,:), allocatable :: dnorm
    integer :: ig, it, ik, is, isgn
    integer :: iglo

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))

    if (proc0) then
       pflux = 0.0;   qflux = 0.0;   vflux = 0.0
       pmflux = 0.0;  qmflux = 0.0;  vmflux = 0.0
       pbflux = 0.0;  qbflux = 0.0;  vbflux = 0.0
       theta_pflux = 0.0  ; theta_pmflux = 0.0  ; theta_pbflux = 0.0
       theta_vflux = 0.0  ; theta_vmflux = 0.0  ; theta_vbflux = 0.0
       theta_qflux = 0.0  ; theta_qmflux = 0.0  ; theta_qbflux = 0.0
    end if

    if (adiabatic_option_switch == adiabatic_option_noJ) then
       dnorm = 1.
    else
       do ik = 1, naky
          do it = 1, ntheta0
! BD bug fix.  actually done in summer sometime, but in main source 2.2.04
!             dnorm(:,it,ik) = grho/bmag/gradpar
             dnorm(:,it,ik) = 1./bmag/gradpar
          end do
       end do
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          dnorm(:,it,ik) = dnorm(:,it,ik)*delthet*woutunits(ik)/sqrt(2.0)
       end do
    end do
    
! anorm is only actually used for QL estimates.  It is not up to date.

! weird definition was here.  changed 7.13.01
    anorm = sum(dnorm*(abs(phi)**2 + abs(apar)**2))

    if (fphi > epsilon(0.0)) then
       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*aj0
       end do
       call get_flux (phi, pflux, theta_pflux, anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call get_flux (phi, qflux(:,:,:,1), theta_qflux(:,:,1), anorm, dnorm)

       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*2.*vpa(:,isgn,:)**2*aj0
       end do
       call get_flux (phi, qflux(:,:,:,2), theta_qflux(:,:,2), anorm, dnorm)

       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*vperp2*aj0
       end do
       call get_flux (phi, qflux(:,:,:,3), theta_qflux(:,:,3), anorm, dnorm)

       do isgn = 1, 2
          g0(:,isgn,:) = g(:,isgn,:)*aj0*vpac(:,isgn,:)
       end do
       call get_flux (phi, vflux, theta_vflux, anorm, dnorm)

    else
       pflux = 0.
       qflux = 0.
!       qflux_par = 0.
!       qflux_perp = 0.
       vflux = 0.
    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo)
          end do
       end do
       call get_flux (apar, pmflux, theta_pmflux, anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call get_flux (apar, qmflux(:,:,:,1), theta_qmflux(:,:,1), anorm, dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo) &
                  *2.*vpa(:,isgn,iglo)**2
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,2), theta_qmflux(:,:,2), anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo) &
                  *vperp2(:,iglo)
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,3), theta_qmflux(:,:,3), anorm, dnorm)
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm &
                  *vpa(:,isgn,iglo)*vpac(:,isgn,iglo)
          end do
       end do
       call get_flux (apar, vmflux, theta_vmflux, anorm, dnorm)
    else
       pmflux = 0.
       qmflux = 0.
!       qmflux_par = 0.
!       qmflux_perp = 0.
       vmflux = 0.
    end if

    if (fbpar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = g(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
          end do
       end do
       call get_flux (bpar, pbflux, theta_pbflux, anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call get_flux (bpar, qbflux(:,:,:,1), theta_qbflux(:,:,1), anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = g(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz &
                    *2.*vpa(:,isgn,iglo)**2
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,2), theta_qbflux(:,:,2), anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = g(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz &
                    *vperp2(:,iglo)
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,3), theta_qbflux(:,:,3), anorm, dnorm)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = g(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo) &
                  *spec(is)%tz*vpac(:,isgn,iglo)
          end do
       end do
       call get_flux (bpar, vbflux, theta_vbflux, anorm, dnorm)
    else
       pbflux = 0.
       qbflux = 0.
!       qbflux_par = 0.
!       qbflux_perp = 0.
       vbflux = 0.
    end if

    deallocate (dnorm)
  end subroutine flux

  subroutine get_flux (fld, flx, theta_flx, anorm, dnorm)
    use theta_grid, only: ntgrid, delthet
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use theta_grid, only: grho
    use species, only: nspec
    use mp, only: proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (in out) :: flx
    real, dimension (-ntgrid:,:), intent (in out) :: theta_flx
    real :: anorm
    real, dimension (-ntgrid:,:,:) :: dnorm
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is, ig

    if (abs(anorm) < 2.0*epsilon(0.0)) then
       if (proc0) flx = 0.0
       return
    end if

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (g0, total)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                wgt = sum(dnorm(:,it,ik)*grho)
                flx(it,ik,is) = sum(aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     *dnorm(:,it,ik)*aky(ik))/wgt/anorm
             end do
          end do
       end do

       ! factors of 0.5, kxfac, included 8.16.00
       ! flx = flx*0.5*kxfac 
       ! factor of kxfac was wrong!  Fixed 2.2.04 BD
       flx = flx*0.5

       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                do ig = -ntgrid, ntgrid
                   wgt = sum(dnorm(:,it,ik)*grho)*anorm*delthet(ig)
                   theta_flx(ig,is) = theta_flx(ig,is) + &
                        aimag(total(ig,it,ik,is)*conjg(fld(ig,it,ik)) &
                        *dnorm(ig,it,ik)*aky(ik))/wgt
                end do
             end do
          end do
       end do

       theta_flx = theta_flx*0.5
    end if

    deallocate (total)
  end subroutine get_flux
!>GGH
!=============================================================================
! Density: Calculate Density perturbations
!=============================================================================
  subroutine get_dens_vel (dv, dvk, phi, apar, bpar, phinew, aparnew, bparnew)
    use constants, only: pi
    use dist_fn_arrays, only: aj0, aj1, vpar, vperp2, g, gnew
    use gs2_heating, only: dens_vel_diagnostics
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use kt_grids, only: ntheta0, naky,aky
    use le_grids, only: integrate_moment
    use mp, only: proc0
    use nonlinear_terms, only: nonlin
    use run_parameters, only: fphi, fbpar
    use species, only: spec,nspec
    use theta_grid, only: jacob, delthet, ntgrid
    implicit none
    !Passed
    type (dens_vel_diagnostics) :: dv
    type (dens_vel_diagnostics), dimension(:,:) :: dvk
    complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar, phinew, aparnew, bparnew
    !Local 
    integer :: isgn, iglo, ig, is, ik, it            !Indices
    complex :: j0phi                                 !J_0 q phi/T 
    complex :: j1bpar                               !2 J_1 v_perp^2 A_perp
    complex :: hvpar                                 !Dist Func h v_par
    complex :: hvperp                                !Dist Func h v_perp
    complex :: hh                                    !Dist Func h 
    complex, dimension(:,:,:,:), allocatable :: tot  !Integrated DF
    real :: fac2                                     !Factor
    real, dimension (:), allocatable :: wgt

    !Allocate tot variable for each calculation
    allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))

    !Set up weighting factors for z-sums on proc0
    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
          wgt(ig) = delthet(ig)*jacob(ig)
       end do
       wgt = wgt/sum(wgt)         
    endif

    !Parallel velocity perturbation--------------------------------------------
    !Loop through and constuct time and z-centered integrands
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid, ntgrid-1
             !Construct time and space centered v_par h
             !NOTE: h=g + J_0 q phi/T + 2 J_1 v_perp^2  A_perp
             j0phi   = aj0(ig,iglo)*phi(ig,it,ik)*spec(is)%zt
             j1bpar = aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpar(ig,it,ik)
             
             hvpar = vpar(ig,isgn,iglo)*( &
                  g(ig,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             j0phi   = aj0(ig+1,iglo)*phi(ig+1,it,ik)*spec(is)%zt
             j1bpar = aj1(ig+1,iglo)*2.0*vperp2(ig+1,iglo)*bpar(ig+1,it,ik)
             
             hvpar = hvpar +  vpar(ig+1,isgn,iglo)*( &
                  g(ig+1,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             j0phi   = aj0(ig,iglo)*phinew(ig,it,ik)*spec(is)%zt
             j1bpar = aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bparnew(ig,it,ik)
             
             hvpar = hvpar + vpar(ig,isgn,iglo)*( &
                  gnew(ig,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             j0phi   = aj0(ig+1,iglo)*phinew(ig+1,it,ik)*spec(is)%zt
             j1bpar = aj1(ig+1,iglo)*2.0*vperp2(ig+1,iglo)*bparnew(ig+1,it,ik)
             
             hvpar = hvpar + vpar(ig+1,isgn,iglo)*( &
                  gnew(ig+1,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             hvpar = 0.25 * hvpar 

             g0(ig,isgn,iglo) = spec(is)%dens*hvpar

          end do
       end do
    end do

    !Integrate over velocity
    call integrate_moment (g0, tot)
    
    !Sum over theta (z-direction) to get perpendicular Fourier components only
    if (proc0) then
       do ik = 1, naky
          !GGH QUES: What is fac2 all about?
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             !Skip mean value (k=0) components
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   dvk(it,ik)%dvpar(is) = dvk(it,ik)%dvpar(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                dv % dvpar(is) = dv%dvpar(is) + dvk(it,ik)%dvpar(is)
             end do
          end do
       end do
    end if

    !Perpendicular velocity perturbation---------------------------------------
    !Loop through and constuct time and z-centered integrands
    !Non-adiabatic part
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid, ntgrid-1
             !Construct time and space centered v_par h
             !NOTE: h=g + J_0 q phi/T + 2 J_1 v_perp^2  A_perp
             j0phi   = aj0(ig,iglo)*phi(ig,it,ik)*spec(is)%zt
             j1bpar = aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpar(ig,it,ik)
             
             hvperp = sqrt(vperp2(ig,iglo))*( &
                  g(ig,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             j0phi   = aj0(ig+1,iglo)*phi(ig+1,it,ik)*spec(is)%zt
             j1bpar = aj1(ig+1,iglo)*2.0*vperp2(ig+1,iglo)*bpar(ig+1,it,ik)
             
             hvperp = hvperp +  sqrt(vperp2(ig+1,iglo))*( &
                  g(ig+1,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             j0phi   = aj0(ig,iglo)*phinew(ig,it,ik)*spec(is)%zt
             j1bpar = aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bparnew(ig,it,ik)
             
             hvperp = hvperp + sqrt(vperp2(ig,iglo))*( &
                  gnew(ig,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             j0phi   = aj0(ig+1,iglo)*phinew(ig+1,it,ik)*spec(is)%zt
             j1bpar = aj1(ig+1,iglo)*2.0*vperp2(ig+1,iglo)*bparnew(ig+1,it,ik)
             
             hvperp = hvperp + sqrt(vperp2(ig+1,iglo))*( &
                  gnew(ig+1,isgn,iglo) + j0phi*fphi + j1bpar*fbpar)
             
             hvperp = 0.25 * hvperp 

             g0(ig,isgn,iglo) = spec(is)%dens*hvperp

          end do
       end do
    end do

    !Integrate over velocity
    call integrate_moment (g0, tot)
    
    !Add in adiabatic part 
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do ig= -ntgrid, ntgrid-1
          !Time and z-center J_0 Phi
          j0phi= 0.25*spec(is)%zt* &
               ( aj0(ig,iglo)  *(phi(ig,it,ik)  + phinew(ig,it,ik)  ) + &
                 aj0(ig+1,iglo)*(phi(ig+1,it,ik)+ phinew(ig+1,it,ik)) )

          tot(ig,it,ik,is)=tot(ig,it,ik,is)+ sqrt(pi/2.)*spec(is)%stm* &
               spec(is)%dens*(1.-j0phi)               
       enddo
    enddo

    !Sum over theta (z-direction) to get perpendicular Fourier components only
    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
           do it = 1, ntheta0
             !Skip mean value (k=0) components
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   dvk(it,ik)%dvperp(is) = dvk(it,ik)%dvperp(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                dv % dvperp(is) = dv%dvperp(is) + dvk(it,ik)%dvperp(is)
             end do
          end do
       end do
    end if
    !Density perturbation------------------------------------------------------
    !Loop through and constuct time and z-centered integrands
    !Non-adiabatic part
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid, ntgrid-1
             !Construct time and space centered v_par h
             !NOTE: h=g + J_0 q phi/T + 2 J_1 v_perp^2  A_perp
             j0phi   = aj0(ig,iglo)*phi(ig,it,ik)*spec(is)%zt
             j1bpar = aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpar(ig,it,ik)
             
             hh = g(ig,isgn,iglo) + j0phi*fphi + j1bpar*fbpar
             
             j0phi   = aj0(ig+1,iglo)*phi(ig+1,it,ik)*spec(is)%zt
             j1bpar = aj1(ig+1,iglo)*2.0*vperp2(ig+1,iglo)*bpar(ig+1,it,ik)
             
             hh = hh +  &
                  g(ig+1,isgn,iglo) + j0phi*fphi + j1bpar*fbpar
             
             j0phi   = aj0(ig,iglo)*phinew(ig,it,ik)*spec(is)%zt
             j1bpar = aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bparnew(ig,it,ik)
             
             hh = hh + &
                  gnew(ig,isgn,iglo) + j0phi*fphi + j1bpar*fbpar
             
             j0phi   = aj0(ig+1,iglo)*phinew(ig+1,it,ik)*spec(is)%zt
             j1bpar = aj1(ig+1,iglo)*2.0*vperp2(ig+1,iglo)*bparnew(ig+1,it,ik)
             
             hh = hh +  &
                  gnew(ig+1,isgn,iglo) + j0phi*fphi + j1bpar*fbpar
             
             hh = 0.25 * hh 

             g0(ig,isgn,iglo) = spec(is)%dens*hh

          end do
       end do
    end do

    !Integrate over velocity
    call integrate_moment (g0, tot)
    
    !Add in adiabatic part 
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do ig= -ntgrid, ntgrid-1
          !Time and z-center J_0 Phi
          j0phi= 0.25*spec(is)%zt* &
               ( aj0(ig,iglo)  *(phi(ig,it,ik)  + phinew(ig,it,ik)  ) + &
                 aj0(ig+1,iglo)*(phi(ig+1,it,ik)+ phinew(ig+1,it,ik)) )

          tot(ig,it,ik,is)=tot(ig,it,ik,is) - spec(is)%dens*j0phi              
       enddo
    enddo

    !Sum over theta (z-direction) to get perpendicular Fourier components only
    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             !Skip mean value (k=0) components
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   dvk(it,ik)%dn(is) = dvk(it,ik)%dn(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                dv % dn(is) = dv%dn(is) + dvk(it,ik)%dn(is)
             end do
          end do
       end do
    end if

    deallocate (tot)
  end subroutine get_dens_vel
!=============================================================================
! Density: Calculate Density perturbations
!=============================================================================
  subroutine get_jext(j_ext)
    use dist_fn_arrays, only: kperp2
    use kt_grids, only: ntheta0, naky,aky
    use mp, only: proc0
    use theta_grid, only: jacob, delthet, ntgrid
    use antenna, only: antenna_apar
    implicit none
    !Passed
    real, dimension(:,:) ::  j_ext
    !Local 
    complex, dimension(:,:,:), allocatable :: j_extz
    integer :: ig,ik, it                             !Indices
    real :: fac2                                     !Factor
    real, dimension (:), allocatable :: wgt
    

    !Get z-centered j_ext at current time
    allocate (j_extz(-ntgrid:ntgrid, ntheta0, naky)) ; j_extz = 0.
    call antenna_apar (kperp2, j_extz)       

    !Set weighting factor for z-averages
!    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       !GGH NOTE: Here wgt is 1/(2*ntgrid)
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
          wgt(ig) = delthet(ig)*jacob(ig)
       end do
       wgt = wgt/sum(wgt)         
!    endif

    !Take average over z
    do ig=-ntgrid, ntgrid-1
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             j_ext(it,ik)=j_ext(it,ik)+real(j_extz(ig,it,ik))*wgt(ig)*fac2
          end do
       end do
    enddo

    deallocate (wgt,j_extz)

  end subroutine get_jext
!<GGH
!=============================================================================
  subroutine get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)
    use mp, only: proc0, iproc
    use constants, only: pi, zi
    use kt_grids, only: ntheta0, naky, aky, akx
    use dist_fn_arrays, only: vpa, vpac, aj0, aj1, vperp2, g, gnew, kperp2
    use gs2_heating, only: heating_diagnostics
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
    use le_grids, only: integrate_moment
    use species, only: spec, nspec,has_electron_species
    use theta_grid, only: jacob, delthet, ntgrid
    use run_parameters, only: fphi, fapar, fbpar, tunits, beta, tnorm, tite
    use gs2_time, only: code_dt
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_apar, a_ext_data
    use hyper, only: D_v, D_eta, nexp
    implicit none
    type (heating_diagnostics) :: h
    type (heating_diagnostics), dimension(:,:) :: hk
!    complex, dimension (-ntgrid:,:,:), pointer :: hh, hnew
    complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar, phinew, aparnew, bparnew
    complex, dimension(:,:,:,:), allocatable :: tot
    complex, dimension(:,:,:), allocatable :: epar, bpardot, apardot, phidot, j_ext, a_ext_old, a_ext_new
    complex :: fac, hfac, pfac, jpold, afac, chi, havg, gavg
    complex :: bfac                                   !B_|| intermediate calculation
    complex :: chidot, eparavg, j0phidot, j0phiavg, j1bparavg, j0aparavg
    complex :: phi_m, apar_m, bpar_m, pstar, pstardot, gdot, hdot, hold, adot
    complex :: phi_avg, apar_avg, bpar_avg, bperp_m, bperp_avg
    complex :: de, denew, chinew,chi1,chi2, havg2,havg2new,h2dot
    real, dimension (:), allocatable :: wgt
    real :: fac2, dtinv, akperp4 , akperp6, akperp2
    integer :: isgn, iglo, ig, is, ik, it, ie

    g0(ntgrid,:,:) = 0.

! ==========================================================================
! Ion/Electron heating------------------------------------------------------
! ==========================================================================

    allocate ( phidot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (apardot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (bpardot(-ntgrid:ntgrid, ntheta0, naky))

    call dot ( phi,  phinew,  phidot, fphi)
    call dot (apar, aparnew, apardot, fapar)
    call dot (bpar, bparnew, bpardot, fbpar)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Next two calls make the variables g, gnew = h, hnew 
!!! until the end of this procedure!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call g_adjust (g,    phi,    bpar,    fphi, fbpar)
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2
          
          do ig=-ntgrid, ntgrid-1
             
             chidot = aj0(ig,iglo)*(phidot(ig,it,ik) &
                  - vpac(ig,isgn,iglo) * spec(is)%stm * apardot(ig,it,ik)) &
                  + aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpardot(ig,it,ik)*spec(is)%tz
             
             hdot = fdot (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo), dtinv)
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))
             
! First term on RHS and LHS of Eq B-10 of H1:

             g0(ig,isgn,iglo) = spec(is)%dens*conjg(havg)* &
                  (chidot*spec(is)%z-hdot*spec(is)%temp)

          end do
       end do
    end do

    deallocate (phidot, apardot, bpardot)

    allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))

    call integrate_moment (g0, tot)

    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
          wgt(ig) = delthet(ig)*jacob(ig)
       end do
       wgt = wgt/sum(wgt)         

       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                    hk(it,ik) % heating(is) = hk(it,ik) % heating(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2 
                end do
                h % heating(is) = h % heating(is) + hk(it,ik) % heating(is)
             end do
          end do
       end do
    end if

! ==========================================================================
! Antenna contribution to entropy equations for each species
! ==========================================================================

    allocate (a_ext_old(-ntgrid:ntgrid, ntheta0, naky))
    allocate (a_ext_new(-ntgrid:ntgrid, ntheta0, naky))

    call a_ext_data (a_ext_old, a_ext_new)

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2
          do ig=-ntgrid, ntgrid-1
             
             adot = fdot (a_ext_old(ig  ,it,ik), &
                          a_ext_old(ig+1,it,ik), &
                          a_ext_new(ig  ,it,ik), &
                          a_ext_new(ig+1,it,ik), dtinv)
             
             havg = favg (g   (ig  ,isgn,iglo)*aj0(ig  ,iglo)*vpa(ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo)*aj0(ig+1,iglo)*vpa(ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo)*aj0(ig  ,iglo)*vpa(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo)*aj0(ig+1,iglo)*vpa(ig+1,isgn,iglo))
             
             g0(ig,isgn,iglo) = conjg(havg) * spec(is)%z * spec(is)%stm * spec(is)%dens * adot

          end do
       end do
    end do

    deallocate (a_ext_old, a_ext_new)

    call integrate_moment (g0, tot)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                    hk(it,ik) % S_ext(is) = hk(it,ik) % S_ext(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                h % S_ext(is) = h % S_ext(is) + hk(it,ik) % S_ext(is)
             end do
          end do
       end do
    end if

! ==========================================================================
! Antenna Power and B-field contribution to E and E_dot---------------------
! ==========================================================================
    if (proc0) then
!       allocate (epar (-ntgrid:ntgrid, ntheta0, naky)) ; epar = 0.
       allocate (j_ext(-ntgrid:ntgrid, ntheta0, naky)) ; j_ext = 0.
       call antenna_apar (kperp2, j_ext)       
!       call get_epar (phi, apar, phinew, aparnew, epar)
       
       if (beta > epsilon(0.)) then
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0
                
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                
                do ig=-ntgrid, ntgrid-1
                   
                   !GGH Time and space averaged estimate of d/dt(apar)
                   apar_m = fdot (apar   (ig  ,it,ik), &
                                  apar   (ig+1,it,ik), &
                                  aparnew(ig  ,it,ik), &
                                  aparnew(ig+1,it,ik), dtinv)*fapar

! J_ext.E when driving antenna only includes A_parallel:

                   hk(it,ik) % antenna = hk(it, ik) % antenna + real(conjg(j_ext(ig,it,ik))*apar_m)*wgt(ig)*fac2

                   !GGH Time and space averaged estimate of d/dt(bperp)
                   bperp_m = fdot (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                   aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), dtinv) * fapar

                   !GGH Time and space averaged estimate of d/dt(bpar)
                   bpar_m = fdot (bpar   (ig  ,it,ik), &
                                  bpar   (ig+1,it,ik), &
                                  bparnew(ig  ,it,ik), &
                                  bparnew(ig+1,it,ik), dtinv)*fbpar

                   !GGH Time and space averaged estimate of bperp
                   bperp_avg = favg (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                     aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik))) * fapar

                   !GGH Time and space averaged estimate of bpar
                   bpar_avg = favg (bpar   (ig  ,it,ik), &
                                    bpar   (ig+1,it,ik), &
                                    bparnew(ig  ,it,ik), &
                                    bparnew(ig+1,it,ik)) * fbpar

! 1/2 * d/dt B**2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        real(0.25 * conjg(bperp_m)*bperp_avg + conjg(bpar_m)*bpar_avg) &
                        * wgt(ig)*fac2*(2.0/beta)

! B**2/2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
                   hk(it,ik) % energy = hk(it,ik) % energy &
                        + 0.5*real((0.25*conjg(bperp_avg)*bperp_avg + conjg(bpar_avg)*bpar_avg)) &
                        * wgt(ig)*fac2*(2.0/beta)

                   !Eapar = int k_perp^2 A_par^2/(8 pi)                   
                   hk(it,ik) % eapar = hk(it,ik) % eapar &
                        + 0.5*real(0.25*conjg(bperp_avg)*bperp_avg) * wgt(ig)*fac2*(2.0/beta)

                   !Ebpar = int B_par^2/(8 pi)
                   hk(it,ik) % ebpar = hk(it,ik) % ebpar &
                        + 0.5*real(conjg(bpar_avg)*bpar_avg) * wgt(ig)*fac2*(2.0/beta)
                   
                end do
                h % antenna = h % antenna + hk(it, ik) % antenna
                h % eapar = h % eapar + hk(it, ik) % eapar
                h % ebpar = h % ebpar + hk(it, ik) % ebpar
             end do
          end do
       else
          hk % antenna = 0.
          h  % antenna = 0.
          hk % energy_dot = 0.
          hk % energy = 0.
          hk % eapar = 0.
          h  % eapar = 0.
          hk % ebpar = 0.
          h  % ebpar = 0.
       end if
       deallocate (j_ext)
    end if

! ==========================================================================
! Finish E_dot--------------------------------------------------------------
! ==========================================================================

!GGH Include response of Boltzmann species for single-species runs

    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   phi_m   = fdot (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik), dtinv)

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        fphi * real(conjg(phi_avg)*phi_m) &
                        * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to E_dot for single species runs---------------------
 
!GGH New E_dot calc
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             !Calculate old fluctuating energy de
             havg = favg_x (g(ig  ,isgn,iglo), &
                            g(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phi(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phi(ig  ,it,ik), &
                               phi(ig+1,it,ik)) * fphi * spec(is)%zt

             de=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)

            !Calculate new fluctuating energy denew
             havg = favg_x (gnew(ig  ,isgn,iglo), &
                            gnew(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phinew(ig  ,it,ik), &
                               phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             denew=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg) 

             !Set g0 as the change of energy (denew-de)/dt
             g0(ig,isgn,iglo) = fdot_t(de,denew,dtinv)

          end do
       end do
    end do
    !GGH -END new e_dot calc

   call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy_dot = h % energy_dot + hk(it,ik) % energy_dot
          end do
       end do
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ==========================================================================
! Gradient Contributions to Heating-----------------------------------------
! ==========================================================================

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle

       do isgn=1,2
          do ig=-ntgrid, ntgrid-1

             chi = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phi   (ig+1,it,ik),  &
                         aj0(ig  ,iglo)*phinew(ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phinew(ig+1,it,ik)) &
                         *fphi

!!GGH Bug fix: The apar part should be subtracted (because chi= phi - v|| A|| + B||)
             chi = chi - &
                  favg (aj0(ig  ,iglo)*apar   (ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*apar   (ig+1,it,ik)*vpac(ig+1,isgn,iglo),  &
                        aj0(ig  ,iglo)*aparnew(ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*aparnew(ig+1,it,ik)*vpac(ig+1,isgn,iglo)) &
                        *spec(is)%stm*fapar
                
             chi = chi + &
                  favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                        aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                        *spec(is)%tz*fbpar

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             g0(ig,isgn,iglo) = zi * wstar(ik,ie,is)/code_dt * conjg(havg)*chi &
                  * spec(is)%dens
            
          end do
       end do
    end do

    call integrate_moment (g0, tot)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % gradients(is) = hk(it,ik) % gradients(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                h % gradients(is) = h % gradients(is) + hk(it,ik) % gradients(is)
             end do
          end do
       end do
    end if
! ==========================================================================
! Hyperviscosity------------------------------------------------------------
! ==========================================================================

    if (D_v > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1
                
                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo)) 

                j0phiavg = favg (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phi(ig+1,it,ik), &
                                 aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                j1bparavg= favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                                 aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                                 *fbpar
!Set g0 for hyperviscous heating
                g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_v*akperp4* &
                     ( conjg(havg)*havg - conjg(havg)*j0phiavg - &
                     conjg(havg)*j1bparavg)

             end do
          end do
       end do

       call integrate_moment (g0, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hypervisc(is) = hk(it,ik) % hypervisc(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hypervisc(is) = h % hypervisc(is) + hk(it,ik) % hypervisc(is)
                end do
             end do
          end do
       end if

    end if !End Hyperviscous Heating Calculation


! ==========================================================================
! Hyperresistivity------------------------------------------------------------
! ==========================================================================
 
    if (D_eta > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1
                
                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo)) 

                j0aparavg = favg (aj0(ig  ,iglo) * apar(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * apar(ig+1,it,ik), &
                                 aj0(ig  ,iglo)  * aparnew(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * aparnew(ig+1,it,ik)) & 
                                 * fapar * spec(is)%zstm * vpac(ig,isgn,iglo)

!Set g0 for hyperresistive heating
                g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_eta*akperp4* &
                     conjg(havg)*j0aparavg

             end do
          end do
       end do

       call integrate_moment (g0, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hyperres(is) = hk(it,ik) % hyperres(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hyperres(is) = h % hyperres(is) + hk(it,ik) % hyperres(is)
                end do
             end do
          end do
       end if

    end if !End Hyperresistivity Heating Calculation

!==========================================================================
!Finish Energy-------------------------------------------------------------
!==========================================================================

!GGH Calculate hs2-------------------------------------------------------------
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             g0(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg)
          end do
       end do
    end do

    call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec             
                do ig = -ntgrid, ntgrid-1

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % hs2(is) = hk(it,ik) % hs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                end do
             end do
             h % hs2(:) = h % hs2(:) + hk(it,ik) % hs2(:)
          end do
       end do
    end if

!Calculate phis2-------------------------------------------------------------
    if (proc0) then
       do ik=1,naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1,ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do ig=-ntgrid, ntgrid-1
                do is = 1, nspec
                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % phis2(is) = hk(it,ik) % phis2(is)  &
                        +0.5*spec(is)%temp*spec(is)%dens*real(conjg(phi_avg)*phi_avg) &
                        * wgt(ig) * fac2
                enddo
             end do
             h % phis2(:) = h % phis2(:) + hk(it,ik) % phis2(:)
          end do
       end do
    endif

! Calculate delfs2 (rest of energy)-----------------------------------------------

!GGH  Include response of Boltzmann species for single species runs
    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy = hk(it,ik) % energy + &
                        fphi * real(conjg(phi_avg)*phi_avg) &
                        * 0.5 * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to energy for single species runs---------------------

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             j0phiavg = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik), &
                              aj0(ig+1,iglo)*phi   (ig+1,it,ik), &  
                              aj0(ig  ,iglo)*phinew(ig  ,it,ik), &  
                              aj0(ig+1,iglo)*phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg (phi   (ig  ,it,ik), &
                             phi   (ig+1,it,ik), &
                             phinew(ig  ,it,ik), &
                             phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             g0(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)
          end do
       end do
    end do

    call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec             
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy = hk(it,ik) % energy  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                   !Delfs2 = int_r int_v T/F0 dfs^2/2
                   hk(it,ik) % delfs2(is) = hk(it,ik) % delfs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy = h % energy + hk(it,ik) % energy
             h % delfs2(:) = h % delfs2(:) + hk(it,ik) % delfs2(:)
          end do
       end do
       deallocate (wgt)
    end if

    deallocate (tot)

!!
!! Put g, gnew back to their usual meanings
!!
    call g_adjust (g,    phi,    bpar,    -fphi, -fbpar)
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar)

  end subroutine get_heat
!==============================================================================

  subroutine get_stress (rstress, ustress)

    use mp, only: proc0
    use kt_grids, only: ntheta0, akx
    use dist_fn_arrays, only: vpa
    use gs2_layouts, only: g_lo, is_idx, il_idx, ie_idx
    use le_grids, only: integrate, orbit_avg, integrate_stress, ng2
    use species, only: spec, nspec
    use theta_grid, only: drhodpsi, bmag, ntgrid
    use constants, only: zi
    implicit none
    complex, dimension (:,:) :: rstress, ustress
    real, dimension(nspec) :: wgt
    real :: fac1, fac2
    integer :: isgn, ie, is, il, iglo, ig, it

    wgt = spec%z*spec%dens
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is=is_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid
             g0(ig,isgn,iglo) = gnl_1(ig,isgn,iglo)*wgt(is)
          end do
       end do
    end do

    call integrate_stress (g0, rstress)

    if (proc0) then
       do is=1,nspec
          do it = 2, ntheta0
             rstress(it, is) = rstress (it, is) / (zi*akx(it))
          end do
       end do
    end if

    wgt = spec%mass*spec%dens
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il=il_idx(g_lo, iglo)
       if (il > ng2) cycle
       is=is_idx(g_lo, iglo)
       ie=ie_idx(g_lo, iglo)
       
       do ig=-ntgrid,ntgrid
          fac1=(vpa(ig,1,iglo)/bmag(ig)-orbit_avg(il,ie))*drhodpsi
          fac2=(vpa(ig,2,iglo)/bmag(ig)+orbit_avg(il,ie))*drhodpsi
          g0(ig,1,iglo) = gnl_1(ig,1,iglo)*fac1
          g0(ig,2,iglo) = gnl_1(ig,2,iglo)*fac2
       end do
    end do

    call integrate_stress (g0, ustress)

  end subroutine get_stress

  subroutine neoclassical_flux (pflux, qflux)
    use species, only: nspec
    use theta_grid, only: ntgrid, gradpar, bmag, delthet
    use kt_grids, only: naky, ntheta0
    use le_grids, only: e, integrate
    use dist_fn_arrays, only: g
    use gs2_layouts, only: g_lo, geint_lo, ik_idx, it_idx, ie_idx, is_idx
    use gs2_layouts, only: idx, proc_id, idx_local
    use mp, only: proc0, send, receive, broadcast
    use run_parameters, only: woutunits
    use gs2_time, only: code_dt
    use constants
    implicit none
    complex, dimension (:,:,:), intent (out) :: pflux, qflux
    real, dimension (:,:,:), allocatable :: dnorm
    complex, dimension (:,:,:), allocatable :: anorm
    complex, dimension (:,:), allocatable :: geint
    integer :: ig, ik, it, ie, is, ng
    integer :: iglo, igeint
    complex :: x

    allocate (dnorm (-ntgrid:ntgrid,ntheta0,naky))
    allocate (anorm (ntheta0,naky,nspec))
    allocate (geint (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_alloc))
    ng = ntgrid

    if (proc0) then
       pflux = 0.0
       qflux = 0.0
    end if

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             dnorm(ig,it,ik) = 1.0/gradpar(ig)/bmag(ig)*woutunits(ik)/sqrt(2.0)
          end do
       end do
    end do
    dnorm(-ntgrid,:,:) = 0.5*dnorm(-ntgrid,:,:)
    dnorm(ntgrid,:,:)  = 0.5*dnorm(ntgrid,:,:)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       g0(:,1,iglo) = dnorm(:,it,ik)
       g0(:,2,iglo) = dnorm(:,it,ik)
    end do
    call integrate (g0, geint)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             igeint = idx(geint_lo,ik,it,is)
             if (idx_local(geint_lo,igeint)) then
                anorm(it,ik,is) &
                     = sum(geint(-ng:ng-1,igeint)*delthet(-ng:ng-1))
             end if
             call broadcast (anorm(it,ik,is), proc_id(geint_lo,igeint))
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       g0(:,1,iglo) = dnorm(:,it,ik)*zi*wdrift(:,iglo)/code_dt*g(:,1,iglo)
       g0(:,2,iglo) = dnorm(:,it,ik)*zi*wdrift(:,iglo)/code_dt*g(:,2,iglo)
    end do
    call integrate (g0, geint)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             igeint = idx(geint_lo,ik,it,is)
             if (proc0) then
                if (idx_local(geint_lo,igeint)) then
                   pflux(it,ik,is) &
                        = sum(geint(-ng:ng-1,igeint) &
                              *delthet(-ng:ng-1))/anorm(it,ik,is)
                else
                   call receive (pflux(it,ik,is), proc_id(geint_lo,igeint))
                end if
             else if (idx_local(geint_lo,igeint)) then
                x = sum(geint(-ng:ng-1,igeint)*delthet(-ng:ng-1)) &
                     /anorm(it,ik,is)
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
       do ik = 1, naky
          do it = 1, ntheta0
             igeint = idx(geint_lo,ik,it,is)
             if (proc0) then
                if (idx_local(geint_lo,igeint)) then
                   qflux(it,ik,is) &
                        = sum(geint(-ng:ng-1,igeint) &
                              *delthet(-ng:ng-1))/anorm(it,ik,is)
                else
                   call receive (qflux(it,ik,is), proc_id(geint_lo,igeint))
                end if
             else if (idx_local(geint_lo,igeint)) then
                x = sum(geint(-ng:ng-1,igeint)*delthet(-ng:ng-1)) &
                     /anorm(it,ik,is)
                call send (x, 0)
             end if
          end do
       end do
    end do

    deallocate (dnorm, anorm, geint)
  end subroutine neoclassical_flux

  subroutine init_fieldcheck
  end subroutine init_fieldcheck

  subroutine finish_fieldcheck
  end subroutine finish_fieldcheck

  subroutine fieldcheck (phi, apar, bpar)
    use file_utils, only: open_output_file, close_output_file
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, aky, ntheta0, theta0
    use dist_fn_arrays, only: kperp2
    use mp, only: proc0
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (:,:,:), allocatable :: fieldeq, fieldeqa, fieldeqp
    integer :: ig, ik, it, unit

    allocate (fieldeq (-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqa(-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqp(-ntgrid:ntgrid,ntheta0,naky))
    call getfieldeq (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)

    if (proc0) then
       call open_output_file (unit, ".fieldcheck")
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
!                write (unit,*) "j,theta,ik,aky,it,theta0 ", &
!                     ig, theta(ig), ik, aky(ik), it, theta0(it,ik)
!                write (unit,"(8(1pe12.5,1x))") phi(ig,it,ik),fieldeq(ig,it,ik),&
!                     gamtot(ig,it,ik),kperp2(ig,it,ik)
!                write (unit,"(8(1pe12.5,1x))") &
!                     kperp2(ig,it,ik)*gridfac1(ig,it,ik)*apar(ig,it,ik), &
!                     fieldeqa(ig,it,ik), gamtot1(ig,it,ik), kperp2(ig,it,ik)
!                write (unit,"(8(1pe12.5,1x))") &
!                     gridfac1(ig,it,ik)*bpar(ig,it,ik), fieldeqp(ig,it,ik), &
!                     gamtot2(ig,it,ik), kperp2(ig,it,ik)
                write (unit,"(8(1pe12.5,1x))") &
                     gamtot(ig,it,ik), gamtot1(ig,it,ik), gamtot2(ig,it,ik), &
                     kperp2(ig,it,ik)
             end do
          end do
       end do
       call close_output_file (unit)
     end if
     deallocate (fieldeq, fieldeqa, fieldeqp)
  end subroutine fieldcheck

  subroutine init_vortcheck
  end subroutine init_vortcheck

  subroutine finish_vortcheck
  end subroutine finish_vortcheck

  subroutine vortcheck (phi, bpar)
    use file_utils, only: open_output_file, close_output_file
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag, gbdrift, gbdrift0, cvdrift, cvdrift0
    use le_grids, only: e, al, anon, integrate_species
    use kt_grids, only: naky, aky, ntheta0, theta0
    use dist_fn_arrays, only: aj0, aj1, vperp2, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use run_parameters, only: fphi, fbpar
    use gs2_time, only: code_dt
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, dimension (nspec) :: wgt
    complex, dimension (:,:,:,:), allocatable :: apchk
    integer :: iglo, ig, ik, it, il, ie, is
    real :: temp, z
    real, dimension (-ntgrid:ntgrid) :: vplus, cvtot, gbtot, delwd
    integer :: unit
    
    allocate (apchk (-ntgrid:ntgrid,ntheta0,naky,2))
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
       cvtot = cvdrift + theta0(it,ik)*cvdrift0
       gbtot = gbdrift + theta0(it,ik)*gbdrift0
       delwd = temp/z*code_dt*(gbtot-cvtot)*vplus/2.0

       g0(:,1,iglo) &
         = (gnew(:,1,iglo) &
            +anon(ie,is)*2.0*vperp2(:,iglo)*aj1(:,iglo)*fbpar*bpar(:,it,ik) &
            +fphi*z*anon(ie,is)*phi(:,it,ik)*aj0(:,iglo)/temp) &
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

       g0(:,1,iglo) = bpar(:,it,ik)*fbpar*vperp2(:,iglo)*temp/z
    end do
    wgt = spec%dens*spec%z
    call integrate_species (g0, wgt, apchk(:,:,:,2))

    call open_output_file (unit, ".vortcheck")
    do ik = 1, naky
       do it = 1, ntheta0
          write (unit,*) 'aky=',aky(ik), ' theta0=',theta0(it,ik)
          do ig = -ntgrid, ntgrid
             write (unit,*) apchk(ig,it,ik,1), apchk(ig,it,ik,2)
          end do
       end do
    end do
    call close_output_file (unit)
    deallocate (apchk)
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

    use dist_fn_arrays, only: gnew, g
    initializing  = .true.
    initialized = .false.
    
    wdrift = 0.
    wdriftttp = 0.
    a = 0.
    b = 0.
    r = 0.
    ainv = 0.
    gnew = 0.
    g0 = 0.
    g = 0.

  end subroutine reset_init

  subroutine reset_physics

    call init_wstar

  end subroutine reset_physics

  subroutine write_f (last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx, ie_idx
    use gs2_layouts, only: idx_local, proc_id
    use le_grids, only: al, e, forbid
    use theta_grid, only: bmag
    use gs2_time, only: user_time

    integer :: iglo, ik, it, is
    integer :: ie, il, ig, ie_last
    integer, save :: unit
    real :: vpa, vpe
    complex, dimension(2) :: gtmp
    logical :: first = .true.
    logical, intent(in)  :: last 

    if (proc0) then
       if (first) then 
          call get_unused_unit (unit)
          call open_output_file (unit, ".dist")
          first = .false.
       end if
    endif
    ie_last = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo, iglo) ; if (ik /= 1) cycle
       it = it_idx(g_lo, iglo) ; if (it /= 1) cycle
       is = is_idx(g_lo, iglo) !; if (is /= 1) cycle
       ie = ie_idx(g_lo, iglo) 
       ig = 0
       il = il_idx(g_lo, iglo)
       if (idx_local (g_lo, ik, it, il, ie, is)) then
          if (proc0) then 
             gtmp = g0(ig,:,iglo)
          else
             call send (g0(ig,:,iglo), 0)
          endif
       else if (proc0) then
          call receive (gtmp, proc_id(g_lo, iglo))
       endif
       if (proc0) then
          if (ie /= ie_last) then
             write (unit, *)
             write (unit, fmt="('# Energy = ',e16.10)") e(ie,is)
          end if
          if (.not. forbid(ig,il)) then
             vpa = sqrt(e(ie,is)*max(0.0, 1.0-al(il)*bmag(ig)))
             vpe = sqrt(e(ie,is)*al(il)*bmag(ig))
!             write (unit, "(8(1x,e12.6))") vpe,vpa,gtmp(1),-vpa,gtmp(2)             
             gtmp = gtmp !* exp(-e(ie,is))
             write (unit, "(8(1x,e12.6),1x,i2)") vpe,vpa,gtmp(1),-vpa,gtmp(2),user_time,is
          end if
       end if
       ie_last = ie
    end do
!    if (proc0) write (unit, *)
    if (proc0) write (unit, *)
    if (last .and. proc0) call close_output_file (unit)
    
  end subroutine write_f

  subroutine boundary(linked)

    logical :: linked

    call init_dist_fn
    linked = boundary_option_switch == boundary_option_linked

  end subroutine boundary

! This subroutine only returns epar correctly for linear runs.
  subroutine get_epar (phi, apar, phinew, aparnew, epar)

    use theta_grid, only: ntgrid, delthet, gradpar
    use run_parameters, only: tunits, fphi, fapar
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    complex, dimension(-ntgrid:,:,:) :: phi, apar, phinew, aparnew, epar
    complex :: phi_m, apar_m

    integer :: ig, ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid-1
             ! ignoring decentering in time and space for now
             phi_m = 0.5*(phi(ig+1,it,ik)-phi(ig,it,ik) + &
                  phinew(ig+1,it,ik)-phinew(ig,it,ik))*fphi
             apar_m = 0.5*(aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
                  -apar(ig+1,it,ik)-apar(ig,it,ik))*fapar
             
             epar(ig,it,ik) = -phi_m/delthet(ig)* &
                  0.5*(abs(gradpar(ig)) + abs(gradpar(ig+1))) &
                  -apar_m/tunits(ik)/code_dt
          end do
       end do
    end do    

  end subroutine get_epar

  subroutine find_leftmost_link (iglo, iglo_left, ipleft)

    integer, intent (in) :: iglo
    integer, intent (in out) :: iglo_left, ipleft
    integer :: iglo_star, iglo_left_star, ipleft_star

    iglo_star = iglo
    do 
       call get_left_connection (iglo_star, iglo_left_star, ipleft_star)
    
       if (ipleft_star == -1) exit
       iglo_star = iglo_left_star
       iglo_left = iglo_left_star
       ipleft = ipleft_star
    end do

  end subroutine find_leftmost_link

  subroutine find_rightmost_link (iglo, iglo_right, ipright)

    integer, intent (in) :: iglo
    integer, intent (in out) :: iglo_right, ipright
    integer :: iglo_star, iglo_right_star, ipright_star

    iglo_star = iglo
    do 
       call get_right_connection (iglo_star, iglo_right_star, ipright_star)
    
       if (ipright_star == -1) exit
       iglo_star = iglo_right_star
       iglo_right = iglo_right_star
       ipright = ipright_star
    end do

  end subroutine find_rightmost_link

  subroutine timer (i, place)
    
    character (len=10) :: zdate, ztime, zzone
    character (*) :: place
    integer, intent (in) :: i
    integer, dimension(8) :: ival
    real, save :: told=0., tnew=0.
    
    call date_and_time (zdate, ztime, zzone, ival)
    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
    if (i == 0) told = -1.
    if (told > 0.) then
       print *, ': Elapsed time = ',tnew-told,' seconds in '//trim(place)
    end if
    told = tnew
  end subroutine timer

  subroutine dot (a, anew, adot, fac)

! Get a theta-centered and time-centered estimate of the time derivative 
! of a field.
! 
! tunits(ky) == 1. unless the "omega_*" units are chosen.
! omega_* units normalize time by an additional factor of ky.
!    

    use run_parameters, only: tunits
    use gs2_time, only: code_dt
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid

    implicit none
    complex, intent (in), dimension (-ntgrid:,:,:) :: a, anew
    complex, intent (out), dimension (-ntgrid:,:,:) :: adot
    real, intent (in) :: fac
    real :: dtinv
    integer :: ig, it, ik

    do ik=1,naky
       dtinv = 1./(code_dt*tunits(ik))
       do it=1,ntheta0
          do ig=-ntgrid,ntgrid-1
             adot(ig,it,ik) = 0.5*fac*(anew(ig+1,it,ik)+anew(ig,it,ik) - &
                  (a(ig+1,it,ik)+a(ig,it,ik)))*dtinv
          end do
       end do
    end do
    
  end subroutine dot

  complex function fdot (fl, fr, fnewl, fnewr, dtinv)

    complex, intent (in) :: fl, fr, fnewl, fnewr
    real, intent (in) :: dtinv
    
    fdot = 0.5*(fnewl+fnewr-fl-fr)*dtinv

  end function fdot

! construct time and space-centered quantities
! (should use actual bakdif and fexpr values?)
!
  complex function favg (fl, fr, fnewl, fnewr)

    complex, intent (in) :: fl, fr, fnewl, fnewr
    
    favg = 0.25*(fnewl+fnewr+fl+fr)

  end function favg
 
   complex function fdot_t (f,fnew, dtinv)

    complex, intent (in) :: f, fnew
    real, intent (in) :: dtinv
    
    fdot_t = (fnew-f)*dtinv

  end function fdot_t

 complex function favg_x (fl, fr)

    complex, intent (in) :: fl, fr
    
    favg_x = 0.5*(fl+fr)

  end function favg_x

end module dist_fn


!             chidot = 0.5*(aj0(ig,iglo)+aj0(ig+1,iglo))* &
!                  (phidot(ig,it,ik)-spec(is)%stm*vpa(ig,isgn,iglo)*apardot(ig,it,ik)) &
!                  +2.0*vperp2(ig,iglo)*spec(is)%tz*bpardot(ig,it,ik)*aj1(ig,iglo)

!             chidot = aj0(ig,iglo)* &
!                  (phidot(ig,it,ik)-spec(is)%stm*vpa(ig,isgn,iglo)*apardot(ig,it,ik)) &
!                  +2.0*vperp2(ig,iglo)*spec(is)%tz*bpardot(ig,it,ik)*aj1(ig,iglo)

!             chidot = aj0(ig,iglo)*spec(is)%stm*vpa(ig,isgn,iglo)*apardot(ig,it,ik) &
!                   +aj1(ig,iglo)*2.0*vperp2(ig,iglo)*spec(is)%tz*bpardot(ig,it,ik)


!    case (heating_option_qdh)
!       
!       allocate ( apardot(-ntgrid:ntgrid, ntheta0, naky))
!
!       do ik=1,naky
!          !GGH: Get correct time interval for the chosen ky value
!          dtinv = 1./(code_dt*tunits(ik))
!          do it=1,ntheta0
!             do ig=-ntgrid,ntgrid-1
!                apardot(ig,it,ik) = 0.5*fapar*(aparnew(ig+1,it,ik)+aparnew(ig,it,ik) - &
!                     (apar(ig+1,it,ik)+apar(ig,it,ik)))*dtinv
!             end do
!          end do
!       end do
!
!! GGH Loop through s,x,y,sgn,z and determine heating calculation
!       do iglo=g_lo%llim_proc, g_lo%ulim_proc
!          is = is_idx(g_lo, iglo)
!          it = it_idx(g_lo, iglo)
!          ik = ik_idx(g_lo, iglo)
!! GGH If nonlinear run, do not do this for mean values (kx=ky=0)
!          if (nonlin .and. it == 1 .and. ik == 1) cycle
!
!          dtinv = 1./(code_dt*tunits(ik))
!          do isgn=1,2
!             do ig=-ntgrid, ntgrid-1
!
!                !GGH NOTE: Why is the following not averaging vperp2 for ig, ig+1
!                !      and current and new times?
!                !BD: vperp2 is a coordinate, constant in time.  spatial variation, 
!                !      we could look into
!                pstar = 0.25*(fphi*aj0(ig,iglo)*(phi(ig,it,ik)+phi(ig+1,it,ik) &
!                     + phinew(ig,it,ik)+phinew(ig+1,it,ik)) &
!                     + fbpar*aj1(ig,iglo)*2.0*vperp2(ig,iglo) &
!                     * (bpar(ig,it,ik)+bpar(ig+1,it,ik) &
!                     + bparnew(ig,it,ik)+bparnew(ig+1,it,ik))*spec(is)%tz)
!
!                gdot = 0.5*(gnew(ig+1,isgn,iglo)+gnew(ig,isgn,iglo) - &
!                           (g   (ig+1,isgn,iglo)+g   (ig,isgn,iglo)))*dtinv
!             
!                fac = 0.25*(g   (ig,isgn,iglo) + g   (ig+1,isgn,iglo) &
!                          + gnew(ig,isgn,iglo) + gnew(ig+1,isgn,iglo))
!
!                g0(ig,isgn,iglo) = - spec(is)%z*(conjg(fac)* &
!                     aj0(ig,iglo)*vpac(ig,isgn,iglo)*spec(is)%stm*apardot(ig,it,ik)*fapar &
!                     + conjg(pstar)*gdot) - conjg(fac)*gdot*spec(is)%temp
!!
!! The above is equivalent to what we actually derived except for the g*gdot term:
!!                g0(ig,isgn,iglo) = spec(is)%z*(conjg(fac)*chidot &
!!                     - conjg(pstar)*gdot - conjg(fac)*pstardot)
!!
!             end do
!          end do
!       end do
!
!       deallocate (apardot)
!
!
