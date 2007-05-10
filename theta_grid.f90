module theta_grid_params
  implicit none

  public :: init_theta_grid_params

  real, public :: rhoc, rmaj, r_geo, eps, epsl
  real, public :: qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri

  integer, public :: ntheta, nperiod

  private

contains

  subroutine init_theta_grid_params
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    call read_parameters
  end subroutine init_theta_grid_params

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none

    namelist /theta_grid_parameters/ rhoc, rmaj, r_geo, eps, epsl, &
         qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri, &
         ntheta, nperiod

    rhoc = 0.5
    rmaj = 3.0
    r_geo = 3.0
    eps = 0.3
    epsl = 0.3
    qinp = 1.5
    shat = 0.75
    pk = 0.3
    shift = 0.0
    akappa = 1.0
    akappri = 0.0
    tri = 0.0
    tripri = 0.0
    ntheta = 24
    nperiod = 2
    read (unit=input_unit("theta_grid_parameters"), nml=theta_grid_parameters)
  end subroutine read_parameters

end module theta_grid_params

module theta_grid_gridgen
  implicit none

  public :: theta_grid_gridgen_init
  public :: gridgen_get_grids

  private

  ! knobs
  integer :: npadd
  real :: alknob, epsknob, bpknob, extrknob, tension
  real :: thetamax, deltaw, widthw

contains

  subroutine theta_grid_gridgen_init
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    call read_parameters
  end subroutine theta_grid_gridgen_init

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none
    namelist /theta_grid_gridgen_knobs/ &
         npadd, alknob, epsknob, bpknob, extrknob, tension, thetamax, deltaw, widthw

    npadd = 2
    alknob = 0.0
    epsknob = 1e-5
    bpknob = 1.e-8
    extrknob = 0.0
    tension = 1.0
    thetamax = 0.0
    deltaw = 0.0
    widthw = 1.0
    read (unit=input_unit("theta_grid_gridgen_knobs"), &
         nml=theta_grid_gridgen_knobs)
  end subroutine read_parameters

  subroutine gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
       theta, bset, bmag, &
       gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22, grho)
    use gridgen4mod
    use constants
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: theta
    real, dimension (nbset), intent (in out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    integer :: ntheta_old, ntgrid_old, nbset_old
    real, dimension (-ntgrid:ntgrid) :: thetasave
    real, dimension (ntheta+1) :: thetaold, thetanew
    real, dimension (ntheta+1) :: bmagold, bmagnew
    integer :: i

    ntheta_old = ntheta
    ntgrid_old = ntgrid
    nbset_old = nbset

    thetasave = theta
    thetaold = theta(-ntheta/2:ntheta/2)
    bmagold = bmag(-ntheta/2:ntheta/2)

    call gridgen4_2 (1,ntheta_old+1,thetaold,bmagold, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
         ntheta,nbset,thetanew,bmagnew,bset)

    if (ntheta_old /= ntheta) then
       write(*,*) 'Error in theta_grid_gridgen?'
       write(*,*) 'ntheta_old = ',ntheta_old
       write(*,*) 'ntheta_new = ',ntheta
       write(*,*) 'Stopping this run would be wise.'
    end if

    ! interpolate to new grid
    ntgrid = ntheta/2 + (nperiod-1)*ntheta

    theta(-ntheta/2:ntheta/2-1) = thetanew(1:ntheta)
    theta(ntheta/2) = thetanew(1) + real(2)*pi
    bmag(-ntheta/2:ntheta/2-1) = bmagnew(1:ntheta)
    bmag(ntheta/2) = bmagnew(1)
    do i = 1, nperiod-1
       theta(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) &
            = thetanew(1:ntheta) + real(2*i)*pi
       theta(ntheta/2+i*ntheta) = thetanew(1) + real(2*(i+1))*pi
       theta(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) &
            = thetanew(1:ntheta) - real(2*i)*pi
       bmag(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) = bmagnew(1:ntheta)
       bmag(ntheta/2+i*ntheta) = bmagnew(1)
       bmag(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) = bmagnew(1:ntheta)
    end do

    call regrid (ntgrid_old, thetasave, gradpar, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds2, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds21, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds22, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, grho, ntgrid, theta)

  end subroutine gridgen_get_grids

  subroutine regrid (nold, x, y, nnew, xnew)
    use splines
    implicit none
    integer, intent (in) :: nold
    real, dimension (-nold:nold), intent (in) :: x
    real, dimension (-nold:nold), intent (in out) :: y
    integer, intent (in) :: nnew
    real, dimension (-nnew:nnew), intent (in) :: xnew
    type (spline) :: spl
    integer :: i

    call new_spline (2*nold+1, x(-nold:nold), y(-nold:nold), spl)

    do i = -nnew, nnew
       y(i) = splint(xnew(i), spl)
    end do

    call delete_spline (spl)
  end subroutine regrid

end module theta_grid_gridgen

module theta_grid_salpha
  implicit none

  public :: init_theta_grid_salpha
  public :: salpha_get_sizes
  public :: salpha_get_grids

  private

  ! knobs
  real :: alpmhdfac, alpha1

  ! internal variable
  integer :: model_switch
  integer, parameter :: model_salpha = 1, model_alpha1 = 2, &
       model_nocurve = 3, model_ccurv = 4, model_b2 = 5, model_eps = 6

  real :: shift


contains

  subroutine init_theta_grid_salpha
    use theta_grid_params, only: init_theta_grid_params
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .false.

    call init_theta_grid_params
    call read_parameters
  end subroutine init_theta_grid_salpha

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use theta_grid_params, only: shift_in => shift, alpmhd
    use text_options
    implicit none

    character(20) :: model_option
    type (text_option), dimension (7), parameter :: modelopts = &
         (/ text_option('default', model_salpha), &
            text_option('s-alpha', model_salpha), &
            text_option('alpha1', model_alpha1), &
            text_option('rogers', model_eps), &
            text_option('b2', model_b2), &
            text_option('const-curv', model_ccurv), &
            text_option('no-curvature', model_nocurve) /)

    namelist /theta_grid_salpha_knobs/ alpmhdfac, alpha1, model_option
    integer :: ierr

    alpmhdfac = 0.0
    alpha1 = 0.0
    model_option = 'default'
    read (unit=input_unit("theta_grid_salpha_knobs"), &
         nml=theta_grid_salpha_knobs)

    ierr = error_unit()
    call get_option_value &
         (model_option, modelopts, model_switch, &
         ierr, "model_option in theta_grid_salpha_knobs")

    if (alpmhdfac > epsilon(0.0)) then
       shift = - alpmhd*alpmhdfac
    else
       shift = shift_in
    end if

  end subroutine read_parameters

  subroutine salpha_get_sizes (nthetaout, nperiodout, nbsetout)
    use theta_grid_params, only: ntheta, nperiod
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound when alpha1 model is used
  end subroutine salpha_get_sizes

  subroutine salpha_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, shat, drhodpsi, kxfac)
    use constants
    use theta_grid_params, only: eps, epsl, shat_param => shat, pk, qinp, rhoc
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    real, intent (out) :: shat, drhodpsi, kxfac
    integer :: i

    theta = (/ (real(i)*2.0*pi/real(ntheta), i=-ntgrid,ntgrid) /)

    if (model_switch == model_alpha1) then
       bmag = 1.0-eps*cos(theta)-alpha1*cos(3.0*theta)
    else if (model_switch == model_b2) then
       bmag = 1.0 - eps*cos(theta)
    else
       bmag = 1.0/(1.0 + eps*cos(theta))
    end if

    shat = shat_param
    drhodpsi = qinp/rhoc
    kxfac = 1.0
    select case (model_switch)
    case (model_salpha,model_alpha1,model_b2)
       gbdrift = epsl*(cos(theta) + (shat*theta-shift*sin(theta))*sin(theta))
       gbdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (model_switch == model_b2) then
          gbdrift = gbdrift/bmag**2
          gbdrift0 = gbdrift0/bmag**2
       end if
    case (model_eps)
       gbdrift = epsl*(cos(theta) -eps + (shat*theta-shift*sin(theta))*sin(theta))
       gbdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (model_switch == model_b2) then
          gbdrift = gbdrift/bmag**2
          gbdrift0 = gbdrift0/bmag**2
       end if
    case (model_ccurv,model_nocurve)
       gbdrift = epsl
       gbdrift0 = 0.0

       if (model_switch == model_nocurve) then
          gds2 = 1.0 + (shat*theta)**2
          gds21 = -shat*shat*theta
       else
          gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
          gds21 = -shat*shat*theta
       endif

       gds22 = shat*shat
       grho = 1.0
    end select
    cvdrift = gbdrift
    cvdrift0 = gbdrift0
    gradpar = pk/2.0

    if (model_switch /= model_alpha1) then
       bset = bmag(-ntheta/2:0)
    else
       call theta_grid_gridgen_init
       call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
            gds2, gds21, gds22, grho)
    end if
  end subroutine salpha_get_grids

end module theta_grid_salpha

module theta_grid_eik
  implicit none

  public :: init_theta_grid_eik
  public :: eik_get_sizes
  public :: eik_get_grids

  private

contains

  subroutine init_theta_grid_eik
    use geometry, only: init_theta
    use geometry, only: eikcoefs
    use theta_grid_params, only: init_theta_grid_params, ntheta
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_theta_grid_params

    call read_parameters
    call init_theta (ntheta)
    call eikcoefs
  end subroutine init_theta_grid_eik

  subroutine eik_get_sizes (nthetaout, nperiodout, nbsetout)
    use geometry, only: nperiod
    use theta_grid_params, only: ntheta
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound
  end subroutine eik_get_sizes

  subroutine eik_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag,&
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22,&
            grho, shat, drhodpsi, kxfac)
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    use geometry, only: kxfac_out => kxfac
    use geometry, only: theta_out => theta
    use geometry, only: gradpar_out => gradpar
    use geometry, only: bmag_out => bmag
    use geometry, only: cvdrift_out => cvdrift
    use geometry, only: cvdrift0_out => cvdrift0
    use geometry, only: gbdrift_out => gbdrift
    use geometry, only: gbdrift0_out => gbdrift0
    use geometry, only: gds2_out => gds2
    use geometry, only: gds21_out => gds21
    use geometry, only: gds22_out => gds22
    use geometry, only: grho_out => grho
    use geometry, only: s_hat_input, drhodpsin
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    real, intent (out) :: shat, drhodpsi, kxfac

    theta(-ntgrid:ntgrid) = theta_out(-ntgrid:ntgrid)
    gradpar(-ntgrid:ntgrid) = gradpar_out(-ntgrid:ntgrid)
    bmag(-ntgrid:ntgrid) = bmag_out(-ntgrid:ntgrid)
    cvdrift(-ntgrid:ntgrid) = cvdrift_out(-ntgrid:ntgrid)
    cvdrift0(-ntgrid:ntgrid) = cvdrift0_out(-ntgrid:ntgrid)
    gbdrift(-ntgrid:ntgrid) = gbdrift_out(-ntgrid:ntgrid)
    gbdrift0(-ntgrid:ntgrid) = gbdrift0_out(-ntgrid:ntgrid)
    gds2(-ntgrid:ntgrid) = gds2_out(-ntgrid:ntgrid)
    gds21(-ntgrid:ntgrid) = gds21_out(-ntgrid:ntgrid)
    gds22(-ntgrid:ntgrid) = gds22_out(-ntgrid:ntgrid)
    grho(-ntgrid:ntgrid) = grho_out(-ntgrid:ntgrid)

    call theta_grid_gridgen_init
    call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
         theta, bset, bmag, &
         gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22, &
         grho)
    shat = s_hat_input
    drhodpsi = drhodpsin
    kxfac = kxfac_out
  end subroutine eik_get_grids

  subroutine read_parameters
    use file_utils, only: input_unit
    use geometry, only: nperiod
    use geometry, only: rhoc
    use geometry, only: itor, iflux, irho
    use geometry, only: ppl_eq, gen_eq, vmom_eq, efit_eq, eqfile, local_eq, dfit_eq
    use geometry, only: gs2d_eq
    use geometry, only: equal_arc
    use geometry, only: bishop
    use geometry, only: s_hat_input
    use geometry, only: alpha_input, invLp_input, beta_prime_input, dp_mult
    use geometry, only: rmaj, r_geo
    use geometry, only: shift, qinp, akappa, akappri, tri, tripri
    use geometry, only: delrho, rmin, rmax
    use geometry, only: ismooth, ak0, k1, k2
    use geometry, only: isym, in_nt, writelots
    use theta_grid_params, only: nperiod_in => nperiod
    use theta_grid_params, only: rhoc_in => rhoc
    use theta_grid_params, only: rmaj_in => rmaj, r_geo_in => r_geo
    use theta_grid_params, only: qinp_in => qinp, shat_in => shat
    use theta_grid_params, only: alpmhd_in => alpmhd
    use theta_grid_params, only: shift_in => shift
    use theta_grid_params, only: akappa_in => akappa, akappri_in => akappri
    use theta_grid_params, only: tri_in => tri, tripri_in => tripri
    implicit none

    namelist /theta_grid_eik_knobs/ itor, iflux, irho, &
         ppl_eq, gen_eq, vmom_eq, efit_eq, eqfile, dfit_eq, &
         equal_arc, bishop, local_eq, gs2d_eq, &
         s_hat_input, alpha_input, invLp_input, beta_prime_input, dp_mult, &
         delrho, rmin, rmax, ismooth, ak0, k1, k2, isym, writelots

    nperiod = nperiod_in
    rhoc = rhoc_in
    s_hat_input = shat_in
    alpha_input = alpmhd_in
    rmaj = rmaj_in
    r_geo = r_geo_in
    shift = shift_in
    qinp = qinp_in
    akappa = akappa_in
    akappri = akappri_in
    tri = tri_in
    tripri = tripri_in

    itor = 1
    iflux = 0
    irho = 2
    equal_arc = .true.
    bishop = 5
    dp_mult = 1.0
    delrho = 1e-3
    rmin = 1e-3
    rmax = 1.0
    ismooth = 0
    isym = 0
    in_nt = .false.
    writelots = .false.
    local_eq = .true.

    read (unit=input_unit("theta_grid_eik_knobs"), nml=theta_grid_eik_knobs)
  end subroutine read_parameters
end module theta_grid_eik

module theta_grid_file
  implicit none

  public :: init_theta_grid_file
  public :: file_get_sizes
  public :: file_get_grids

  private

  character(200) :: gridout_file
  real :: shat_input, drhodpsi_input, kxfac_input

contains

  subroutine init_theta_grid_file
    use theta_grid_params, only: init_theta_grid_params
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_theta_grid_params
    call read_parameters
  end subroutine init_theta_grid_file

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none
    namelist /theta_grid_file_knobs/ gridout_file

    gridout_file = "grid.out"
    read (unit=input_unit("theta_grid_file_knobs"), nml=theta_grid_file_knobs)
  end subroutine read_parameters

  subroutine file_get_sizes (ntheta, nperiod, nbset)
    use file_utils, only: get_unused_unit
    implicit none
    integer, intent (out) :: ntheta, nperiod, nbset
    integer :: unit
    character(200) :: line
    integer :: i, ntgrid
    real :: rmaj

    call get_unused_unit (unit)
    open (unit=unit, file=gridout_file, status="old")

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) nbset
    read (unit=unit, fmt="(a)") line
    do i = 1, nbset
       read (unit=unit, fmt="(a)") line
    end do

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) ntgrid, nperiod, ntheta, &
         drhodpsi_input, rmaj, shat_input, kxfac_input

    close (unit=unit)
  end subroutine file_get_sizes

  subroutine file_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, shat, drhodpsi, kxfac)
    use file_utils, only: get_unused_unit
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    real, intent (out) :: shat, drhodpsi, kxfac
    integer :: unit
    character(200) :: line
    integer :: i

    shat = shat_input
    drhodpsi = drhodpsi_input
    kxfac = kxfac_input

    call get_unused_unit (unit)
    open (unit=unit, file=gridout_file, status="old")
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line
    do i = 1, nbset
       read (unit=unit, fmt=*) bset(i) ! actually alambda
    end do
    bset = 1.0/bset ! switch alambda to bset

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gbdrift(i), gradpar(i), grho(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift(i), gds2(i), bmag(i), theta(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gds21(i), gds22(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift0(i), gbdrift0(i)
    end do

    close (unit=unit)
  end subroutine file_get_grids

end module theta_grid_file

module theta_grid
  implicit none

  public :: init_theta_grid
  public :: theta, theta2, delthet, delthet2
  public :: bset
  public :: bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0
  public :: gds2, gds21, gds22, kxfac
  public :: grho
  public :: bmin, bmax, eps, shat, drhodpsi, jacob
  public :: ntheta, ntgrid, nperiod, nbset

  private

  real, dimension (:), allocatable :: theta, theta2, delthet, delthet2
  real, dimension (:), allocatable :: bset
  real, dimension (:), allocatable :: bmag, gradpar
  real, dimension (:), allocatable :: gbdrift, gbdrift0, cvdrift, cvdrift0
  real, dimension (:), allocatable :: gds2, gds21, gds22
  real, dimension (:), allocatable :: grho, jacob
  real :: bmin, bmax, eps, shat, drhodpsi, kxfac
  integer :: ntheta, ntgrid, nperiod, nbset

  ! internal variables
  integer :: eqopt_switch
  integer, parameter :: eqopt_eik = 1, eqopt_salpha = 2, eqopt_file = 3

contains

  subroutine init_theta_grid
    use mp, only: proc0
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    if (proc0) then
       call read_parameters
       call get_sizes
       call allocate_arrays
       call get_grids
       call finish_init
    end if
    call broadcast_results
  end subroutine init_theta_grid

  subroutine broadcast_results
    use mp, only: proc0, broadcast
    implicit none

    call broadcast (bmin)
    call broadcast (bmax)
    call broadcast (eps)
    call broadcast (kxfac)
    call broadcast (ntheta)
    call broadcast (ntgrid)
    call broadcast (nperiod)
    call broadcast (nbset)

    if (.not. proc0) then
       call allocate_arrays
       allocate (theta2(-ntgrid:ntgrid))
       allocate (delthet(-ntgrid:ntgrid))
       allocate (delthet2(-ntgrid:ntgrid))
    end if
    call broadcast (theta)
    call broadcast (theta2)
    call broadcast (delthet)
    call broadcast (delthet2)
    call broadcast (bset)
    call broadcast (bmag)
    call broadcast (gradpar)
    call broadcast (gbdrift)
    call broadcast (gbdrift0)
    call broadcast (cvdrift)
    call broadcast (cvdrift0)
    call broadcast (gds2)
    call broadcast (gds21)
    call broadcast (gds22)
    call broadcast (grho)
    call broadcast (shat)
    call broadcast (jacob)
    call broadcast (drhodpsi)
  end subroutine broadcast_results

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use text_options
    implicit none
    type (text_option), dimension (5), parameter :: eqopts = &
         (/ text_option('default', eqopt_eik), &
            text_option('eik', eqopt_eik), &
            text_option('s-alpha', eqopt_salpha), &
            text_option('grid.out', eqopt_file), &
            text_option('file', eqopt_file) /)
    character(20) :: equilibrium_option
    ! 'default' 'eik': call eikcoefs for parameterized equilibrium
    ! 's-alpha': s-alpha
    ! 'grid.out' 'file': read grid from grid.out file generated by rungridgen
    namelist /theta_grid_knobs/ equilibrium_option
    integer :: ierr

    equilibrium_option = 'default'
    read (unit=input_unit("theta_grid_knobs"), nml=theta_grid_knobs)

    ierr = error_unit()
    call get_option_value &
         (equilibrium_option, eqopts, eqopt_switch, &
         ierr, "equilibrium_option in theta_grid_knobs")
  end subroutine read_parameters

  subroutine allocate_arrays
    implicit none
    allocate (theta(-ntgrid:ntgrid))
    allocate (bset(nbset))
    allocate (bmag(-ntgrid:ntgrid))
    allocate (gradpar(-ntgrid:ntgrid))
    allocate (gbdrift(-ntgrid:ntgrid))
    allocate (gbdrift0(-ntgrid:ntgrid))
    allocate (cvdrift(-ntgrid:ntgrid))
    allocate (cvdrift0(-ntgrid:ntgrid))
    allocate (gds2(-ntgrid:ntgrid))
    allocate (gds21(-ntgrid:ntgrid))
    allocate (gds22(-ntgrid:ntgrid))
    allocate (grho(-ntgrid:ntgrid))
    allocate (jacob(-ntgrid:ntgrid))
  end subroutine allocate_arrays

  subroutine finish_init
    implicit none
    real, dimension (nbset) :: bset_save
    real, dimension (-ntgrid:ntgrid) :: eik_save

    ! in case nbset changes after gridgen
    if (nbset /= size(bset)) then
       bset_save = bset(:nbset)
       deallocate (bset)
       allocate (bset(nbset))
       bset = bset_save
    end if

    ! in case ntgrid changes after gridgen
    if (ntgrid*2+1 /= size(theta)) then
       eik_save = theta(-ntgrid:ntgrid); deallocate (theta)
       allocate (theta(-ntgrid:ntgrid)); theta = eik_save

       eik_save = bmag(-ntgrid:ntgrid); deallocate (bmag)
       allocate (bmag(-ntgrid:ntgrid)); bmag = eik_save

       eik_save = gradpar(-ntgrid:ntgrid); deallocate (gradpar)
       allocate (gradpar(-ntgrid:ntgrid)); gradpar = eik_save

       eik_save = gbdrift(-ntgrid:ntgrid); deallocate (gbdrift)
       allocate (gbdrift(-ntgrid:ntgrid)); gbdrift = eik_save

       eik_save = gbdrift0(-ntgrid:ntgrid); deallocate (gbdrift0)
       allocate (gbdrift0(-ntgrid:ntgrid)); gbdrift0 = eik_save

       eik_save = cvdrift(-ntgrid:ntgrid); deallocate (cvdrift)
       allocate (cvdrift(-ntgrid:ntgrid)); cvdrift = eik_save

       eik_save = cvdrift0(-ntgrid:ntgrid); deallocate (cvdrift0)
       allocate (cvdrift0(-ntgrid:ntgrid)); cvdrift0 = eik_save

       eik_save = gds2(-ntgrid:ntgrid); deallocate (gds2)
       allocate (gds2(-ntgrid:ntgrid)); gds2 = eik_save

       eik_save = gds21(-ntgrid:ntgrid); deallocate (gds21)
       allocate (gds21(-ntgrid:ntgrid)); gds21 = eik_save

       eik_save = gds22(-ntgrid:ntgrid); deallocate (gds22)
       allocate (gds22(-ntgrid:ntgrid)); gds22 = eik_save

       eik_save = grho(-ntgrid:ntgrid); deallocate (grho)
       allocate (grho(-ntgrid:ntgrid)); grho = eik_save
    end if

    bmax = maxval(bmag)
    bmin = minval(bmag)
! ?? check collision operator coding
    eps = 1.0 - sqrt(bmin/bmax)
    eps = 1.0 - 1.0/bmax

    allocate (theta2(-ntgrid:ntgrid))
    allocate (delthet(-ntgrid:ntgrid))
    allocate (delthet2(-ntgrid:ntgrid))

    theta2 = theta*theta
    delthet(:ntgrid-1) = theta(-ntgrid+1:) - theta(:ntgrid-1)
    delthet(ntgrid) = delthet(-ntgrid)
    delthet2 = delthet*delthet

    jacob = 1.0/(drhodpsi*gradpar*bmag)
  end subroutine finish_init

  subroutine get_sizes
    use theta_grid_eik, only: eik_get_sizes, init_theta_grid_eik
    use theta_grid_salpha, only: salpha_get_sizes, init_theta_grid_salpha
    use theta_grid_file, only: file_get_sizes, init_theta_grid_file
    implicit none
    select case (eqopt_switch)
    case (eqopt_eik)
       call init_theta_grid_eik
       call eik_get_sizes (ntheta, nperiod, nbset)
    case (eqopt_salpha)
       call init_theta_grid_salpha
       call salpha_get_sizes (ntheta, nperiod, nbset)
    case (eqopt_file)
       call init_theta_grid_file
       call file_get_sizes (ntheta, nperiod, nbset)
    end select
    ntgrid = ntheta/2 + (nperiod-1)*ntheta
  end subroutine get_sizes

  subroutine get_grids
    use theta_grid_eik, only: eik_get_grids
    use theta_grid_salpha, only: salpha_get_grids
    use theta_grid_file, only: file_get_grids
    implicit none
    select case (eqopt_switch)
    case (eqopt_eik)
       call eik_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
            gds2, gds21, gds22, grho, &
            shat, drhodpsi, kxfac)
    case (eqopt_salpha)
       call salpha_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
            gds2, gds21, gds22, grho, &
            shat, drhodpsi, kxfac)
    case (eqopt_file)
       call file_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
            gds2, gds21, gds22, grho, &
            shat, drhodpsi, kxfac)
    end select
    kxfac = abs(kxfac)
  end subroutine get_grids

end module theta_grid
