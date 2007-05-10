module kt_grids_box
  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids

  private

  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: jtwist
  real :: ly

contains

  subroutine init_kt_grids_box
    use theta_grid, only: init_theta_grid
    use file_utils, only: input_unit
    implicit none
    integer :: naky, ntheta0, nx, ny
    namelist /kt_grids_box_parameters/ naky, ntheta0, ly, nx, ny, jtwist

    call init_theta_grid

    naky = 0
    ntheta0 = 0
    ly = 1.0
    nx = 0
    ny = 0
    jtwist = 1
    read (unit=input_unit("kt_grids_box_parameters"), &
         nml=kt_grids_box_parameters)
    if (naky == 0) naky = (ny-1)/3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_box

  subroutine box_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
  end subroutine box_get_sizes

  subroutine box_get_grids (aky, theta0, akx)
    use theta_grid, only: shat, nperiod
    use constants
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer :: i, naky, ntheta0
    real :: dkx

    naky = size(aky)
    ntheta0 = size(akx)

    dkx = 2.0*pi/real(jtwist)* 2.0*pi/ly*shat

    do i = 1, naky
       aky(i) = real(i-1)*2.*pi/ly
    end do

    do i = 1, (ntheta0+1)/2
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       akx(i) = real(i-ntheta0-1)*dkx
    end do

    do i = 1, ntheta0
       theta0(2:,i) = akx(i)/(aky(2:)*shat)
       theta0(1,i) = 0.0
    end do
  end subroutine box_get_grids

end module kt_grids_box
