module mu_grid

  implicit none

  public :: init_mu_grid, finish_mu_grid
  public :: mu

  real, dimension (:), allocatable :: mu, wgts_mu

contains

  subroutine init_mu_grid

    use gs3_input, only: read_mu_grid_parameters
    use gs3_input, only: nmu, mu_max, ntgrid
    use theta_grid, only: ntheta, bmag, bmax

    implicit none

    integer :: imu

    call read_mu_grid_parameters

    ! allocate arrays and initialize to zero
    if (.not. allocated(mu)) then
       allocate (mu(nmu)) ; mu = 0.0
       allocate (wgts_mu(nmu)) ; wgts_mu = 0.0
    end if

    ! construct mu grid
    do imu = 1, nmu
       mu(imu) = real(imu-1)*mu_max/(nmu-1)
    end do

    ! wgts_mu are the integration weights associated with mu grid points
    wgts_mu = mu_max/nmu

  end subroutine init_mu_grid

  subroutine finish_mu_grid

    implicit none

    if (allocated(mu)) deallocate (mu)
    if (allocated(wgts_mu)) deallocate (wgts_mu)

  end subroutine finish_mu_grid

end module mu_grid
