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
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    real :: kp = -1.
    integer :: in_file
    logical :: exist

    namelist /theta_grid_parameters/ rhoc, rmaj, r_geo, eps, epsl, &
         qinp, shat, alpmhd, pk, shift, akappa, akappri, tri, tripri, &
         ntheta, nperiod, kp

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
    in_file = input_unit_exist("theta_grid_parameters", exist)
!    if (exist) read (unit=input_unit("theta_grid_parameters"), nml=theta_grid_parameters)
    if (exist) read (unit=in_file, nml=theta_grid_parameters)

    if (kp > 0.) pk = 2.*kp

  end subroutine read_parameters

end module theta_grid_params
