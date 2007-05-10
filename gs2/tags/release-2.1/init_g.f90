module init_g
  implicit none

  public :: ginit
  private

  ! knobs
  integer :: ginitopt_switch
  integer, parameter :: ginitopt_default = 1, ginitopt_test1 = 2, &
       ginitopt_xi = 3, ginitopt_xi2 = 4, ginitopt_rh = 5, ginitopt_zero = 6, &
       ginitopt_test3 = 7
  real :: width0, phiinit
  logical :: chop_side

contains

  subroutine ginit
    use mp, only: proc0, broadcast
    implicit none

    if (proc0) call read_parameters

    call broadcast (ginitopt_switch)
    call broadcast (width0)
    call broadcast (phiinit)
    call broadcast (chop_side)

    select case (ginitopt_switch)
    case (ginitopt_default)
       call ginit_default
    case (ginitopt_test1)
       call ginit_test1
    case (ginitopt_xi)
       call ginit_xi
    case (ginitopt_xi2)
       call ginit_xi2
    case (ginitopt_rh)
       call ginit_rh
    case (ginitopt_zero)
       call ginit_zero
    case (ginitopt_test3)
       call ginit_test3
    end select
  end subroutine ginit

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use text_options
    implicit none

    type (text_option), dimension (7), parameter :: ginitopts = &
         (/ text_option('default', ginitopt_default), &
            text_option('test1', ginitopt_test1), &
            text_option('xi', ginitopt_xi), &
            text_option('xi2', ginitopt_xi2), &
            text_option('zero', ginitopt_zero), &
            text_option('test3', ginitopt_test3), &
            text_option('rh', ginitopt_rh) /)
    character(20) :: ginit_option
    namelist /init_g_knobs/ ginit_option, width0, phiinit, chop_side
    integer :: ierr

    ginit_option = "default"
    width0 = 3.5
    phiinit = 1.0
    chop_side = .true.
    read (unit=input_unit("init_g_knobs"), nml=init_g_knobs)

    ierr = error_unit()
    call get_option_value &
         (ginit_option, ginitopts, ginitopt_switch, &
         ierr, "ginit_option in ginit_knobs")
  end subroutine read_parameters

  subroutine ginit_default
    use species, only: nspec, spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0
    use le_grids, only: nlambda, negrid, forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: phi
    integer :: iglo
    integer :: ig, ik, it, il, is

    do it = 1, ntheta0
       do ik = 1, naky
          do ig = -ntgrid, ntgrid
             phi(ig,ik,it) = exp(-((theta(ig)-theta0(ik,it))/width0)**2) &
                  *cmplx(1.0,1.0)
          end do
          if (chop_side) phi(:-1,ik,it) = 0.0
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = -phi(:,ik,it)*spec(is)%z*phiinit
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_default

  subroutine ginit_test1
    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, akr
    use le_grids, only: e, forbid
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    complex, dimension (-ntgrid:ntgrid,naky,ntheta0) :: phi
    integer :: iglo
    integer :: ig, ik, it, il, ie, is

    do it = 1, ntheta0
       do ik = 1, naky
          do ig = -ntgrid, ntgrid
             phi(ig,ik,it) = sin(akr(ig,it)*theta(ig))/(akr(ig,it))*zi
          end do
       end do
    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = -phi(:,ik,it)*spec(is)%z*phiinit*exp(-e(ie,is))
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_test1

  subroutine ginit_xi
    use theta_grid, only: ntgrid, theta, bmag
    use le_grids, only: forbid, al
    use kt_grids, only: theta0
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is
    real, dimension(-ntgrid:ntgrid) :: xi

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       xi = sqrt(max(1.0-bmag*al(il),0.0))
       g(:,1,iglo) = xi*exp(-((theta-theta0(ik,it))/width0)**2)
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = -g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_xi

  subroutine ginit_xi2
    use theta_grid, only: ntgrid, theta, bmag
    use le_grids, only: forbid, al
    use kt_grids, only: theta0
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is
    real, dimension(-ntgrid:ntgrid) :: xi

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       xi = sqrt(max(1.0-bmag*al(il),0.0))
       g(:,1,iglo) = (1.0 - 3.0*xi*xi)*exp(-((theta-theta0(ik,it))/width0)**2)
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_xi2

  subroutine ginit_rh
    use theta_grid, only: ntgrid
    use le_grids, only: forbid, e
    use dist_fn_arrays, only: g, gnew
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use constants
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       g(:,1,iglo) = exp(-e(ie,is))
       where (forbid(:,il)) g(:,1,iglo) = 0.0
       g(:,2,iglo) = g(:,1,iglo)
    end do
    gnew = g
  end subroutine ginit_rh

  subroutine ginit_zero
    use dist_fn_arrays, only: g, gnew
    implicit none
    g = 0.0
    gnew = 0.0
  end subroutine ginit_zero

  subroutine ginit_test3
    use dist_fn_arrays, only: g, gnew, vpa
    use theta_grid, only: ntgrid, delthet, bmag
    use kt_grids, only: akx
    use le_grids, only: forbid
    use theta_grid_params, only: eps, epsl, pk
    use gs2_layouts, only: g_lo, ik_idx, il_idx
    use mp, only: broadcast
    use constants
    implicit none
    integer :: iglo, ik, il
    real :: c1, c2

    call broadcast (epsl)
    call broadcast (eps)
    call broadcast (pk)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       if (any(vpa(-ntgrid:ntgrid-1,:,iglo) == 0.0)) then
          c1 = 0.0
          c2 = 0.0
       else
          c2 = -akx(ik)*epsl/eps/pk &
               *sum(delthet(-ntgrid:ntgrid-1)/bmag(-ntgrid:ntgrid-1))
          c1 = c2/sum(delthet(-ntgrid:ntgrid-1)/vpa(-ntgrid:ntgrid-1,1,iglo))
          c2 = c2/sum(delthet(-ntgrid:ntgrid-1)/vpa(-ntgrid:ntgrid-1,2,iglo))
       end if
       g(:,1,iglo) = -zi*akx(ik)*epsl/eps/pk*vpa(:,1,iglo)/bmag - zi*c1
       g(:,2,iglo) = -zi*akx(ik)*epsl/eps/pk*vpa(:,2,iglo)/bmag - zi*c2
    end do
    gnew = g
  end subroutine ginit_test3

end module init_g
