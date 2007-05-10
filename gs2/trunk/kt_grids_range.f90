module kt_grids_range
  implicit none

  public :: init_kt_grids_range, range_get_sizes, range_get_grids

  private

  integer :: naky_private, ntheta0_private
  real :: aky_min, aky_max, theta0_min, theta0_max

contains

  subroutine init_kt_grids_range
    use file_utils, only: input_unit
    implicit none
    integer :: naky, ntheta0
    namelist /kt_grids_range_parameters/ naky, ntheta0, &
         aky_min, aky_max, theta0_min, theta0_max

    naky = 1
    ntheta0 = 1
    aky_min = 0.0
    aky_max = 0.0
    theta0_min = 0.0
    theta0_max = 0.0
    read (unit=input_unit("kt_grids_range_parameters"), &
         nml=kt_grids_range_parameters)
    naky_private = naky
    ntheta0_private = ntheta0
  end subroutine init_kt_grids_range

  subroutine range_get_sizes (naky, ntheta0)
    implicit none
    integer, intent (out) :: naky, ntheta0
    naky = naky_private
    ntheta0 = ntheta0_private
  end subroutine range_get_sizes

  subroutine range_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    real :: daky, dtheta0
    integer :: i, j, naky, ntheta0

    naky = size(aky)
    ntheta0 = size(akx)

    daky = 0.0
    if (naky > 1) daky = (aky_max - aky_min)/real(naky - 1)
    dtheta0 = 0.0
    if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

    aky = (/ (aky_min + daky*real(i), i=0,naky_private-1) /)
    do j = 1, naky
       theta0(j,:) &
            = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0_private-1) /)
    end do
    akx = 0.0
  end subroutine range_get_grids

end module kt_grids_range
