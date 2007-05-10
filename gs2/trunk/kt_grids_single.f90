module kt_grids_single
  implicit none

  public :: init_kt_grids_single, single_get_sizes, single_get_grids

  private

  real :: aky, theta0, akx

contains

  subroutine init_kt_grids_single
    use file_utils, only: input_unit
    implicit none
    namelist /kt_grids_single_parameters/ aky, theta0, akx

    aky = 0.4
    theta0 = 0.0
    akx = 0.0
    read (unit=input_unit("kt_grids_single_parameters"), &
         nml=kt_grids_single_parameters)
  end subroutine init_kt_grids_single

  subroutine single_get_sizes (naky, ntheta0)
    implicit none
    integer, intent (out) :: naky, ntheta0
    naky = 1
    ntheta0 = 1
  end subroutine single_get_sizes

  subroutine single_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    call get_grids (aky, theta0, akx)
  end subroutine single_get_grids

  subroutine get_grids (aky_out, theta0_out, akx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out
    real, dimension (:,:), intent (out) :: theta0_out
    real, dimension (:), intent (out) :: akx_out
    aky_out(1) = aky
    theta0_out(1,1) = theta0
    akx_out(1) = akx
  end subroutine get_grids

end module kt_grids_single
