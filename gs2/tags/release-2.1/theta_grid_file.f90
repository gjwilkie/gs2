module theta_grid_file
  implicit none

  public :: init_theta_grid_file
  public :: file_get_sizes
  public :: file_get_grids

  private

  character(200) :: gridout_file
  real :: shat_input, drhodpsi_input

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
         drhodpsi_input, rmaj, shat_input

    close (unit=unit)
  end subroutine file_get_sizes

  subroutine file_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
       gds2, gds21, gds22, grho, shat, drhodpsi)
    use file_utils, only: get_unused_unit
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    real, intent (out) :: shat, drhodpsi
    integer :: unit
    character(200) :: line
    integer :: i

    shat = shat_input
    drhodpsi = drhodpsi_input

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
