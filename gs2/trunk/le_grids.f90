module le_grids
  implicit none

  public :: init_le_grids
  public :: integrate, lintegrate, integrate_species
  public :: e, anon, al, delal, jend, forbid
  public :: negrid, nlambda, ng2, lmax
  public :: geint2g, gint2g
  public :: intcheck
  public :: fcheck

  private

  real, dimension (:,:), allocatable :: e, w, anon ! (negrid,nspec)

  real, dimension (:), allocatable :: al, delal ! (nlambda)
  real, dimension (:,:), allocatable :: wl ! (-ntgrid:ntgrid,nlambda)
  integer, dimension (:), allocatable :: jend ! (-ntgrid:ntgrid)
  logical, dimension (:,:), allocatable :: forbid ! (-ntgrid:ntgrid,nlambda)

  integer :: lint_lo, lint_hi, eint_lo, eint_hi
  integer :: geint2g_lo, geint2g_hi, gint2g_lo, gint2g_hi
  complex, dimension (:,:), allocatable :: integration_work
  ! (-ntgrid:ntgrid, -*- processor-dependent -*-)

 ! knobs
  integer :: ngauss, negrid
  real :: ecut, bouncefuzz

  integer :: nlambda, ng2, lmax

contains

  subroutine init_le_grids
    use mp, only: proc0
    use species, only: init_species
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_species
    call init_theta_grid
    call init_kt_grids

    if (proc0) then
       call read_parameters
       call set_grids
    end if
    call broadcast_results
    call init_integrations
  end subroutine init_le_grids

  subroutine broadcast_results
    use mp, only: proc0, broadcast
    use species, only: nspec
    use theta_grid, only: ntgrid
    implicit none
    integer :: il, is

    call broadcast (ngauss)
    call broadcast (negrid)
    call broadcast (ecut)
    call broadcast (bouncefuzz)
    call broadcast (nlambda)
    call broadcast (ng2)
    call broadcast (lmax)

    if (.not. proc0) then
       allocate (e(negrid,nspec), w(negrid,nspec), anon(negrid,nspec))
       allocate (al(nlambda), delal(nlambda))
       allocate (wl(-ntgrid:ntgrid,nlambda))
       allocate (jend(-ntgrid:ntgrid))
       allocate (forbid(-ntgrid:ntgrid,nlambda))
    end if

    call broadcast (al)
    call broadcast (delal)
    call broadcast (jend)
    do is = 1, nspec
       call broadcast (e(:,is))
       call broadcast (w(:,is))
       call broadcast (anon(:,is))
    end do
    do il = 1, nlambda
       call broadcast (wl(:,il))
       call broadcast (forbid(:,il))
    end do
  end subroutine broadcast_results

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none
    namelist /le_grids_knobs/ ngauss, negrid, ecut, bouncefuzz

    ngauss = 5
    negrid = 10
    ecut = 2.5
    bouncefuzz = 1e-5
    read (unit=input_unit("le_grids_knobs"), nml=le_grids_knobs)
  end subroutine read_parameters

  subroutine init_integrations
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use gs2_layouts, only: init_dist_fn_layouts
    use gs2_layouts, only: g_lo, gint_lo, geint_lo, idx_local
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx, idx
    implicit none
    integer :: ig, igint, igeint
    integer :: lint_ulim, geint2g_ulim, eint_ulim, ulim

    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec)

    lint_lo = gint_lo%ulim_world
    lint_hi = gint_lo%llim_world
    geint2g_lo = geint_lo%ulim_world
    geint2g_hi = geint_lo%llim_world
    do ig = g_lo%llim_proc, g_lo%ulim_proc
       igint = idx(gint_lo,ik_idx(g_lo,ig),it_idx(g_lo,ig), &
                   ie_idx(g_lo,ig),is_idx(g_lo,ig))
       lint_lo = min(igint,lint_lo)
       lint_hi = max(igint,lint_hi)
       igeint = idx(geint_lo,ik_idx(gint_lo,igint),it_idx(gint_lo,igint), &
                    is_idx(gint_lo,igint))
       geint2g_lo = min(igeint,geint2g_lo)
       geint2g_hi = max(igeint,geint2g_lo)
    end do
    lint_ulim = lint_hi - lint_lo
    geint2g_ulim = geint2g_hi - geint2g_lo

    eint_lo = geint_lo%ulim_world
    eint_hi = geint_lo%llim_world
    do igint = gint_lo%llim_proc, gint_lo%ulim_proc
       igeint = idx(geint_lo,ik_idx(gint_lo,igint),it_idx(gint_lo,igint), &
                    is_idx(gint_lo,igint))
       eint_lo = min(igeint,eint_lo)
       eint_hi = max(igeint,eint_hi)
    end do
    eint_ulim = eint_hi - eint_lo

    ulim = max(eint_ulim,geint2g_ulim,lint_ulim)
    if (ulim >= 0) then
       allocate (integration_work(-ntgrid:ntgrid,0:ulim))
    end if
  end subroutine init_integrations

  subroutine geint2g (geint, g)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, geint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, is_idx, idx
    use mp, only: broadcast
    implicit none
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (in) :: geint
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g

    complex, dimension (-ntgrid:ntgrid) :: collector
    integer :: iglo, igeint, ik, it, is

    ! essentially, spread geint_lo laid out data into g_lo layout,
    ! replicating along the virtual sign,negrid,nlambda dimensions

    do igeint = geint_lo%llim_world, geint_lo%ulim_world
       if (idx_local(geint_lo,igeint)) then
          collector = geint(:,igeint)
       end if
       call broadcast (collector, proc_id(geint_lo,igeint))
       if (igeint >= geint2g_lo .and. igeint <= geint2g_hi) then
          integration_work(:,igeint-geint2g_lo) = collector
       end if
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       igeint = idx(geint_lo,ik,it,is)
       g(:,1,iglo) = integration_work(:,igeint-geint2g_lo)
       g(:,2,iglo) = integration_work(:,igeint-geint2g_lo)
    end do
  end subroutine geint2g

  subroutine gint2g (gint, g)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, gint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, idx
    use mp, only: broadcast
    implicit none
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (in) :: gint
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g

    complex, dimension (-ntgrid:ntgrid) :: collector
    integer :: iglo, igint, ik, it, ie, is

    ! essentially, spread gint_lo laid out data into g_lo layout,
    ! replicating along the virtual sign,nlambda dimensions

    do igint = gint_lo%llim_world, gint_lo%ulim_world
       if (idx_local(gint_lo,igint)) then
          collector = gint(:,igint)
       end if
       call broadcast (collector, proc_id(gint_lo,igint))
       if (igint >= lint_lo .and. igint <= lint_hi) then
          integration_work(:,igint-lint_lo) = collector
       end if
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       igint = idx(gint_lo,ik,it,ie,is)
       g(:,1,iglo) = integration_work(:,igint-lint_lo)
       g(:,2,iglo) = integration_work(:,igint-lint_lo)
    end do
  end subroutine gint2g

  subroutine lintegrate (g1, gint)
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce
    use gs2_layouts, only: g_lo, gint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx, idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g1
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (out) :: gint

    complex, dimension (-ntgrid:ntgrid) :: collector
    integer :: ig, igint
    integer :: il

    if (lint_hi >= lint_lo) then
       integration_work = 0.0
       do ig = g_lo%llim_proc, g_lo%ulim_proc
          igint = idx(gint_lo,ik_idx(g_lo,ig),it_idx(g_lo,ig), &
                      ie_idx(g_lo,ig),is_idx(g_lo,ig))
          il = il_idx(g_lo,ig)
          integration_work(:,igint-lint_lo) &
               = integration_work(:,igint-lint_lo) &
                 + wl(:,il)*(g1(:,1,ig) + g1(:,2,ig))
       end do
    end if

    collector = 0.0
    do igint = gint_lo%llim_world, gint_lo%ulim_world
       if (igint >= lint_lo .and. igint <= lint_hi) then
          call sum_reduce (integration_work(:,igint - lint_lo), &
                           proc_id(gint_lo, igint))
          if (idx_local(gint_lo, igint)) then
             gint(:,igint) = integration_work(:,igint - lint_lo)
          end if
       else
          call sum_reduce (collector, proc_id(gint_lo, igint))
          if (idx_local(gint_lo, igint)) then
             gint(:,igint) = collector
             collector = 0.0
          end if
       end if
    end do
  end subroutine lintegrate

  subroutine eintegrate (gint, geint)
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce
    use gs2_layouts, only: gint_lo, geint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx, idx
    implicit none
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (in) :: gint
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (out) :: geint

    complex, dimension (-ntgrid:ntgrid) :: collector
    integer :: igint, igeint

    integer :: ie, is

    if (eint_hi >= eint_lo) then
       integration_work = 0.0
       do igint = gint_lo%llim_proc, gint_lo%ulim_proc
          ie = ie_idx(gint_lo,igint)
          is = is_idx(gint_lo,igint)
          igeint = idx(geint_lo,ik_idx(gint_lo,igint),it_idx(gint_lo,igint),is)
          integration_work(:,igeint-eint_lo) &
               = integration_work(:,igeint-eint_lo) + w(ie,is)*gint(:,igint)
       end do
    end if

    collector = 0.0
    do igeint = geint_lo%llim_world, geint_lo%ulim_world
       if (igeint >= eint_lo .and. igeint <= eint_hi) then
          call sum_reduce (integration_work(:,igeint - eint_lo), &
                           proc_id(geint_lo, igeint))
          if (idx_local(geint_lo, igeint)) then
             geint(:,igeint) = integration_work(:,igeint - eint_lo)
          end if
       else
          call sum_reduce (collector, proc_id(geint_lo, igeint))
          if (idx_local(geint_lo, igeint)) then
             geint(:,igeint) = collector
             collector = 0.0
          end if
       end if
    end do
  end subroutine eintegrate

  subroutine lintegrate_old (g1, gint)
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce
    use gs2_layouts, only: g_lo, gint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g1
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (out) :: gint

    complex, dimension (-ntgrid:ntgrid) :: collector
    integer :: ig, igint
    integer :: ik, it, ie, is
    integer :: il

    do igint = gint_lo%llim_world, gint_lo%ulim_world
       ik = ik_idx(gint_lo,igint)
       it = it_idx(gint_lo,igint)
       ie = ie_idx(gint_lo,igint)
       is = is_idx(gint_lo,igint)
       collector = 0.0
       do ig = g_lo%llim_proc, g_lo%ulim_proc
          if (       ik == ik_idx(g_lo,ig) &
               .and. it == it_idx(g_lo,ig) &
               .and. ie == ie_idx(g_lo,ig) &
               .and. is == is_idx(g_lo,ig)) &
          then
             il = il_idx(g_lo,ig)
             collector = collector + wl(:,il)*(g1(:,1,ig) + g1(:,2,ig))
          end if
       end do
       call sum_reduce (collector, proc_id(gint_lo, igint))
       if (idx_local(gint_lo, igint)) gint(:,igint) = collector
    end do
  end subroutine lintegrate_old

  subroutine eintegrate_old (gint, geint)
    use theta_grid, only: ntgrid
    use mp, only: sum_reduce
    use gs2_layouts, only: gint_lo, geint_lo
    use gs2_layouts, only: proc_id, idx_local
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,gint_lo%llim_proc:), intent (in) :: gint
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (out) :: geint

    complex, dimension (-ntgrid:ntgrid) :: collector
    integer :: igint, igeint

    integer :: ik, it, is
    integer :: ie

    do igeint = geint_lo%llim_world, geint_lo%ulim_world
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       collector = 0.0
       do igint = gint_lo%llim_proc, gint_lo%ulim_proc
          if (       ik == ik_idx(gint_lo,igint) &
               .and. it == it_idx(gint_lo,igint) &
               .and. is == is_idx(gint_lo,igint)) &
          then
             ie = ie_idx(gint_lo,igint)
             collector = collector + w(ie,is)*gint(:,igint)
          end if
       end do
       call sum_reduce (collector, proc_id(geint_lo, igeint))
       if (idx_local(geint_lo, igeint)) geint(:,igeint) = collector
    end do
  end subroutine eintegrate_old

  subroutine integrate (g1, geint)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, gint_lo, geint_lo
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g1
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (out) :: geint

    complex, dimension (-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_proc) ::&
         gint

    call lintegrate (g1, gint)
    call eintegrate (gint, geint)
  end subroutine integrate

  subroutine sum_species (geint, weights, total)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: geint_lo, ik_idx, it_idx, is_idx
    use mp, only: sum_allreduce
    implicit none
    complex, dimension (-ntgrid:,geint_lo%llim_proc:), intent (in) :: geint
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total

    integer :: ik, it, is
    integer :: igeint

    total = 0.0
    do igeint = geint_lo%llim_proc, geint_lo%ulim_proc
       ik = ik_idx(geint_lo,igeint)
       it = it_idx(geint_lo,igeint)
       is = is_idx(geint_lo,igeint)
       total(:,ik,it) = total(:,ik,it) + weights(is)*geint(:,igeint)
    end do
    do ik = 1, naky
       do it = 1, ntheta0
          call sum_allreduce (total(:,ik,it))
       end do
    end do
  end subroutine sum_species

  subroutine integrate_species (g, weights, total)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: geint_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total

    complex, &
         dimension (-ntgrid:ntgrid,geint_lo%llim_proc:geint_lo%ulim_proc) :: &
         geint

    call integrate (g, geint)
    call sum_species (geint, weights, total)
  end subroutine integrate_species

  subroutine set_grids
    use species, only: init_species, nspec
    use theta_grid, only: init_theta_grid, ntgrid, nbset, bset, eps
    implicit none

    call init_theta_grid
    call init_species

    allocate (e(negrid,nspec), w(negrid,nspec), anon(negrid,nspec))
    call egridset

    ng2 = 2*ngauss
    if (eps > epsilon(0.0)) then
       nlambda = ng2+nbset
       lmax = nlambda-1
    else
       nlambda = ng2
       lmax = nlambda
    end if
    allocate (al(nlambda), delal(nlambda))
    allocate (wl(-ntgrid:ntgrid,nlambda))
    allocate (jend(-ntgrid:ntgrid))
    allocate (forbid(-ntgrid:ntgrid,nlambda))
    al(ng2+1:nlambda) = 1.0/bset
    call lgridset
    delal(1) = 0.0
    delal(2:) = al(2:) - al(:nlambda-1)
  end subroutine set_grids

  subroutine egridset
    use species, only: nspec, spec, specie, slowing_down_species
    use constants
    implicit none
    real, dimension (16, 16), parameter :: gaus = reshape((/ &
         0.50000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.21132,0.78868,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.11270,0.50000,0.88730,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.06943,0.33001,0.66999,0.93057,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.04691,0.23077,0.50000,0.76923,0.95309,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.03377,0.16940,0.38069,0.61931,0.83060,0.96623,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.02545,0.12923,0.29708,0.50000,0.70292,0.87077,0.97455,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.01986,0.10167,0.23723,0.40828,0.59172,0.76277,0.89833,0.98014, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.01592,0.08199,0.19332,0.33788,0.50000,0.66213,0.80669,0.91802, &
         0.98408,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.01305,0.06747,0.16030,0.28330,0.42557,0.57444,0.71670,0.83971, &
         0.93253,0.98696,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00922,0.04794,0.11505,0.20634,0.31609,0.43739,0.56262,0.68392, &
         0.79366,0.88495,0.95206,0.99078,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00530,0.02771,0.06719,0.12230,0.19106,0.27099,0.35920,0.45250, &
         0.54754,0.64080,0.72901,0.80894,0.87770,0.93282,0.97229,0.99470  &
         /), (/ 16, 16 /))
    real, dimension (16, 16), parameter :: w1 = reshape((/ &
         1.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.50000,0.50000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.27778,.444444,.277778,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.17393,0.32607,0.32607,0.17393,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.11846,0.23931,0.28444,0.23931,0.11846,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.08566,0.18038,0.23395,0.23395,0.18038,0.08566,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.06474,0.13985,0.19091,0.20897,0.19091,0.13985,0.06474,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.05061,0.11119,0.15685,0.18134,0.18134,0.15685,0.11119,0.05061, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.04064,0.09033,0.13031,0.15618,0.16512,0.15618,0.13031,0.09033, &
         0.04064,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.03334,0.07473,0.10955,0.13464,0.14776,0.14776,0.13464,0.10955, &
         0.07473,0.03334,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.02359,0.05347,0.08004,0.10159,0.11675,0.12458,0.12458,0.11675, &
         0.10159,0.08004,0.05347,0.02359,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000, &

         0.01358,0.03113,0.04758,0.06232,0.07480,0.08458,0.09130,0.09473, &
         0.09473,0.09130,0.08458,0.07480,0.06232,0.04758,0.03113,0.01358  &
         /), (/ 16, 16 /))
    ! note that gaus and w1 are transposed wrt original code

    integer :: is, ie
    integer :: ng1
    real :: cut

    if (negrid < 1 .or. negrid > 18) call stop_invalid ("negrid", negrid)
    if (any(gaus(:negrid-2,negrid-2) == 0.0)) then
       call stop_invalid ("negrid", negrid)
    end if
    if (any((spec(:)%type == slowing_down_species))) then
       if (negrid > 16) call stop_invalid ("negrid", negrid)
       if (any(gaus(:negrid-1,negrid-1) == 0.0)) then
          call stop_invalid ("negrid", negrid)
       end if
    end if

    do is = 1, nspec
       if (spec(is)%type == slowing_down_species) then
          e(:negrid-1,is) = gaus(:negrid-1,negrid-1)
          w(:negrid-1,is) = w1(:negrid-1,negrid-1)*e(:negrid-1,is)**2
          anon(:negrid-1,is) = -1.5/e(:negrid-1,is)
          e(negrid,is) = 1.0
          w(negrid,is) = 1e-6
          anon(negrid,is) = 1e6
          w(:,is) = w(:,is)*0.75
       else
          ng1 = max(negrid-2,0)
          if (negrid <= 2) then
             cut = 0.0
          else
             cut = ecut
          end if

          e(ng1+1,is) = cut + 0.58578
          w(ng1+1,is) = 0.853553*exp(-cut)*sqrt(e(ng1+1,is))

          if (negrid > 1) then
             e(ng1+2,is) = cut + 3.41421
             w(ng1+2,is) = 0.146446*exp(-cut)*sqrt(e(ng1+2,is))

             e(:ng1,is) = cut*gaus(:ng1,ng1)**(2.0/3.0)
             w(:ng1,is) = (2.0/3.0)*exp(-e(:ng1,is))*w1(:ng1,ng1)*cut**1.5
          end if
          w(:,is) = w(:,is)*0.5/sqrt(pi)
          anon(:,is) = 1.0

          do ie = 1, negrid
             if(e(ie,is) < spec(is)%e_flat) then
                anon(ie,is) = spec(is)%temp_flat
                write(*,*) 
                print *, "At energy = ",e(ie,is)
                print *, "Flattening energy distribution for species ",is
             endif
          enddo

       end if
    end do
  end subroutine egridset

  subroutine lgridset
    use theta_grid, only: ntgrid, bmag, bmax, eps
    use constants
    implicit none
    real, dimension (12,12), parameter :: xgauss = reshape((/ &
         0.57735,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.33998,0.86114,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.23861,0.66120,0.93246,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.18343,0.52553,0.79666,0.96028,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.14887,0.43339,0.67940,0.86506,0.97390, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.12523,0.36783,0.58732,0.76990,0.90411, &
         0.98156,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.09501,0.28160,0.45802,0.61788,0.75540, &
         0.86563,0.94458,0.98940,0.00000,0.00000, &
         0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.07653,0.22779,0.37371,0.51087,0.63605, &
         0.74633,0.83912,0.91223,0.96397,0.99313, &
         0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.06406,0.19112,0.31504,0.43379,0.54542, &
         0.64809,0.74012,0.82000,0.88642,0.93827, &
         0.97473,0.99519  &
         /), (/ 12, 12 /))
    real, dimension (12,12), parameter :: wgauss = reshape((/ &
         1.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.65215,0.34785,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.46791,0.36076,0.17132,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.36268,0.31370,0.22238,0.10122,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.29552,0.26926,0.21908,0.14945,0.06667, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.24915,0.23349,0.20317,0.16008,0.10694, &
         0.04718,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.18945,0.18260,0.16916,0.14960,0.12463, &
         0.09516,0.06225,0.02715,0.00000,0.00000, &
         0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.15275,0.14917,0.14210,0.13169,0.11819, &
         0.10193,0.08328,0.06267,0.04060,0.01761, &
         0.00000,0.00000, &

         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000,0.00000,0.00000,0.00000, &
         0.00000,0.00000, &

         0.12794,0.12584,0.12167,0.11551,0.10744, &
         0.09762,0.08619,0.07335,0.05930,0.04428, &
         0.02853,0.01234  &
         /), (/ 12, 12 /))
    ! note that xgauss and wgauss are transposed wrt original code

    real, dimension (2*ngauss) :: wx, xx
    real :: ww
    integer :: ig, il

    if (ngauss < 1 .or. ngauss > 12) call stop_invalid ("ngauss", ngauss)
    if (any(xgauss(1:ngauss,ngauss) == 0.0)) then
       call stop_invalid ("ngauss", ngauss)
    end if

    wl = 0.0
    xx(:ngauss) = 0.5*(1.0+xgauss(ngauss:1:-1,ngauss))
    xx(ngauss+1:2*ngauss) = 0.5*(1.0-xgauss(1:ngauss,ngauss))
    wx(:ngauss) = 0.5*wgauss(ngauss:1:-1,ngauss)
    wx(ngauss+1:2*ngauss) = 0.5*wgauss(:ngauss,ngauss)

    al(:ng2) = (1.0 - xx(:ng2)**2)/bmax

    do il = 1, ng2
       do ig = -ntgrid, ntgrid
          wl(ig,il) = wx(il)*2.0*sqrt((bmag(ig)/bmax) &
               *((1.0/bmax-al(il))/(1.0/bmag(ig)-al(il))))
       end do
    end do

    jend = 0
    forbid = .false.

    if (eps <= epsilon(0.0)) return

    jend = ng2 + 1
    do il = ng2+1, nlambda-1
       do ig = -ntgrid, ntgrid
          if (1.0-al(il)*bmag(ig) > -bouncefuzz &
               .and. 1.0-al(il+1)*bmag(ig) > -bouncefuzz) &
          then
             jend(ig) = jend(ig) + 1
             ww = sqrt(max(1.0 -   al(il)*bmag(ig),0.0)) - &
                  sqrt(max(1.0 - al(il+1)*bmag(ig),0.0))
             wl(ig,il)   = wl(ig,il)   + ww
             wl(ig,il+1) = wl(ig,il+1) + ww
          end if
       end do
    end do

    do il = ng2+1, nlambda
       do ig = -ntgrid, ntgrid
          forbid(ig,il) = 1.0 - al(il)*bmag(ig) < -bouncefuzz
       end do
    end do
  end subroutine lgridset

  subroutine stop_invalid (name, val)
    use file_utils, only: error_unit
    implicit none
    character(*), intent (in) :: name
    integer, intent (in) :: val
    integer :: ierr

    ierr = error_unit()
    write (unit=ierr, fmt='("Invalid value for ",a,": ",i5)') name, val
    stop
  end subroutine stop_invalid

  subroutine intcheck (g0)
    use file_utils, only: open_output_file, close_output_file
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, bmag, bmax
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, gint_lo, idx_local, proc_id
    use gs2_layouts, only: idx, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use mp, only: proc0, send, receive
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    integer :: unit
    integer :: iglo, igint, ig, ik, it, il, ie, is
    complex, dimension (-ntgrid:ntgrid,gint_lo%llim_proc:gint_lo%ulim_proc) &
         :: gint
    complex, dimension (negrid) :: dumout

    call open_output_file (unit, ".intcheck")

    write (unit,*) "bmax ", bmax

    write (unit,*) "untrapped check:"
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       if (il <= 2*ngauss) then
          g0(:,1,iglo) = bmax/bmag*(1.0 + zi)*sqrt(1.0 - al(il)*bmag) &
               *(1.0 - bmax*al(il))**(0.5*real(ie-2))
       else
          g0(:,1,iglo) = 0.0
       end if
       g0(:,2,iglo) = g0(:,1,iglo)
    end do
    call lintegrate (g0, gint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             do ig = -ntgrid, ntgrid
                do ie = 1, negrid
                   igint = idx(gint_lo,ik,it,ie,is)
                   if (proc0) then
                      if (idx_local(gint_lo,igint)) then
                         dumout(ie) = gint(ig,igint)
                      else
                         call receive (dumout(ie), proc_id(gint_lo,igint))
                      end if
                   else if (idx_local(gint_lo,igint)) then
                      call send (gint(ig,igint), 0)
                   end if
                end do
                write (unit,*) "is,j,fac ", &
                     is, ig, max(0.0,1.0-bmag(ig)/bmag)**0.5
                write (unit,"(20(1x,1pe12.5))") 0.25*real(ie)*dumout
             end do
          end do
       end do
    end do

    write (unit,*)
    write (unit,*)
    write (unit,*) "trapped check:"
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       if (il <= 2*ngauss) then
          g0(:,1,iglo) = 0.0
       else
          g0(:,1,iglo) = (max(0.0,1.0-al(il)*bmag)**0.5)**real(ie-1)
       end if
       g0(:,2,iglo) = g0(:,1,iglo)
    end do
    call lintegrate (g0, gint)
    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             do ig = -ntgrid, ntgrid
                do ie = 1, negrid
                   igint = idx(gint_lo,ik,it,ie,is)
                   if (proc0) then
                      if (idx_local(gint_lo,igint)) then
                         dumout(ie) = gint(ig,igint)
                      else
                         call receive (dumout(ie), proc_id(gint_lo,igint))
                      end if
                   else if (idx_local(gint_lo,igint)) then
                      call send (gint(ig,igint), 0)
                   end if
                end do
                write (unit,*) "is,j,fac ", &
                     is, ig, max(0.0,1.0-bmag(ig)/bmag)**0.5
                write (unit,"(20(1x,1pe12.5))") &
                     0.5*real(ie)*dumout &
                     /(max(0.0,1.0-bmag(ig)/bmag)**0.5)**real(ie)
             end do
          end do
       end do
    end do

    do ig = -ntgrid, ntgrid
       write (unit,*) "j ", ig
       do il = 1, nlambda
          write (unit,*) "il,wl,sq ",il,wl(ig,il), &
               sqrt(max(0.0,1.0-al(il)*bmag(ig)))
       end do
    end do

    call close_output_file (unit)
  end subroutine intcheck

  subroutine fcheck (g, f)
    use species, only: nspec
    use theta_grid, only: ntgrid, gradpar, bmag, delthet
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use mp, only: sum_allreduce
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,:,:), intent (out) :: f
    integer :: iglo, ig, ik, it, il, ie, is

    f = 0.0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       if (il > 1) then
          f(il,ik,it,is) = f(il,ik,it,is) &
               + pi/(al(il)-al(il-1))*w(ie,is) &
                 *sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1) &
                      *(g(:ntgrid-1,1,iglo) + g(:ntgrid-1,2,iglo)) &
                      *sqrt(max(0.0,1.0-al(il)*bmag(:ntgrid-1)))) &
                 /sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1))
       end if
       if (il < nlambda) then
          f(il+1,ik,it,is) = f(il+1,ik,it,is) &
               - pi/(al(il+1)-al(il))*w(ie,is) &
                 *sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1) &
                      *(g(:ntgrid-1,1,iglo) + g(:ntgrid-1,2,iglo)) &
                      *sqrt(max(0.0,1.0-al(il)*bmag(:ntgrid-1)))) &
                 /sum(delthet(:ntgrid-1)/gradpar(:ntgrid-1)/bmag(:ntgrid-1))
       end if
    end do

    do is = 1, nspec
       do it = 1, ntheta0
          do ik = 1, naky
             call sum_allreduce (f(:,ik,it,is))
          end do
       end do
    end do
  end subroutine fcheck

end module le_grids
