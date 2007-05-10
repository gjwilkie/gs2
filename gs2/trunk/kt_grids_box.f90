module kt_grids_box
  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids

  private

  integer :: naky_private, ntheta0_private
  real :: ly, lx

contains

  subroutine init_kt_grids_box
    use file_utils, only: input_unit
    implicit none
    integer :: naky, ntheta0
    namelist /kt_grids_box_parameters/ naky, ntheta0, ly, lx

    read (unit=input_unit("kt_grids_box_parameters"), &
         nml=kt_grids_box_parameters)
    naky_private = naky
    ntheta0_private = ntheta0
  end subroutine init_kt_grids_box

  subroutine box_get_sizes (naky, ntheta0)
    implicit none
    integer, intent (out) :: naky, ntheta0
    naky = naky_private
    ntheta0 = ntheta0_private
  end subroutine box_get_sizes

  subroutine box_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer :: i, naky, ntheta0

    naky = size(aky)
    ntheta0 = size(akx)

    aky = (/ (real(2*(i-1))*pi/ly, i=1,naky) /)
    akx(:ntheta0/2) = (/ (real(2*(i-1))*pi/lx, i=1,ntheta0/2) /)
    akx(ntheta0/2+1:) &
         = (/ (real(2*(i-ntheta0-1))*pi/lx, i=ntheta0/2+1,ntheta0) /)

    do i = 1, ntheta0
       theta0(2:,i) = akx(i)/(aky(2:)*shat)
       theta0(1,i) = 0.0
    end do
  end subroutine box_get_grids

end module kt_grids_box
