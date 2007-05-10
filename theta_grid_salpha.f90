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
       model_nocurve = 3, model_b2 = 4

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
    type (text_option), dimension (5), parameter :: modelopts = &
         (/ text_option('default', model_salpha), &
            text_option('s-alpha', model_salpha), &
            text_option('alpha1', model_alpha1), &
            text_option('b2', model_b2), &
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
       gds2, gds21, gds22, shat, drhodpsi)
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
         gds2, gds21, gds22
    real, intent (out) :: shat, drhodpsi
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
    if (model_switch /= model_nocurve) then
       gbdrift = epsl*(cos(theta) + (shat*theta-shift*sin(theta))*sin(theta))
       gbdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       if (model_switch == model_b2) then
          gbdrift = gbdrift/bmag**2
          gbdrift0 = gbdrift0/bmag**2
       end if
    else
       gbdrift = epsl
       gbdrift0 = 0.0
       gds2 = 1.0 + (shat*theta)**2
       gds21 = -shat*shat*theta
       gds22 = shat*shat
    end if
    cvdrift = gbdrift
    cvdrift0 = gbdrift0
    gradpar = pk/2.0

    if (model_switch /= model_alpha1) then
       bset = bmag(-ntheta/2:0)
    else
       call theta_grid_gridgen_init
       call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22)
    end if
  end subroutine salpha_get_grids

end module theta_grid_salpha
