module kt_grids_xbox
  implicit none

  public :: init_kt_grids_xbox, xbox_get_sizes, xbox_get_grids

  private

  integer :: ntheta0_private, nx_private
  real :: lx, aky_private

contains

  subroutine init_kt_grids_xbox
    use file_utils, only: input_unit
    implicit none
    integer :: ntheta0, nx
    real :: aky
    namelist /kt_grids_xbox_parameters/ ntheta0, lx, aky, nx

    ntheta0 = 1
    lx = 1.0
    aky = 0.2
    nx = 0
    read (unit=input_unit("kt_grids_xbox_parameters"), &
         nml=kt_grids_xbox_parameters)
    ntheta0_private = ntheta0
    aky_private = aky
    nx_private = nx
  end subroutine init_kt_grids_xbox

  subroutine xbox_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny
    naky = 1
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = 0
  end subroutine xbox_get_sizes

  subroutine xbox_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer :: i, naky, ntheta0

    aky(1) = aky_private

    ntheta0 = size(akx)
    akx(:(ntheta0+1)/2) = (/ (real(2*(i-1))*pi/lx, i=1,(ntheta0+1)/2) /)
    akx((ntheta0+1)/2+1:) &
         = (/ (real(2*(i-ntheta0-1))*pi/lx, i=(ntheta0+1)/2+1,ntheta0) /)
    theta0(1,:) = akx(:)/(aky(1)*shat)
  end subroutine xbox_get_grids

end module kt_grids_xbox
