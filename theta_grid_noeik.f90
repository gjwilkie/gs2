module theta_grid_eik
  implicit none

  public :: init_theta_grid_eik
  public :: eik_get_sizes
  public :: eik_get_grids

  private

contains

  subroutine init_theta_grid_eik
    implicit none
    print *, "theta_grid_eik not available"
    stop
  end subroutine init_theta_grid_eik

  subroutine eik_get_sizes (nthetaout, nperiodout, nbsetout)
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout
    stop
  end subroutine eik_get_sizes

  subroutine eik_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag,&
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22,&
            grho, shat, drhodpsi)
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22, grho
    real, intent (out) :: shat, drhodpsi
    stop
  end subroutine eik_get_grids

end module theta_grid_eik
