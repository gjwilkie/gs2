module theta_grid_gridgen
  implicit none

  public :: theta_grid_gridgen_init
  public :: gridgen_get_grids

  private

  ! knobs
  integer :: npadd
  real :: alknob, epsknob, extrknob, tension
  real :: thetamax, deltaw, widthw

contains

  subroutine theta_grid_gridgen_init
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    call read_parameters
  end subroutine theta_grid_gridgen_init

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none
    namelist /theta_grid_gridgen_knobs/ &
         npadd, alknob, epsknob, extrknob, tension, thetamax, deltaw, widthw

    npadd = 2
    alknob = 0.0
    epsknob = 1e-5
    extrknob = 0.0
    tension = 1.0
    thetamax = 0.0
    deltaw = 0.0
    widthw = 1.0
    read (unit=input_unit("theta_grid_gridgen_knobs"), &
         nml=theta_grid_gridgen_knobs)
  end subroutine read_parameters

  subroutine gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
       theta, bset, bmag, &
       gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22)
    use constants
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: theta
    real, dimension (nbset), intent (in out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, &
         gds2, gds21, gds22
    integer :: ntheta_old, ntgrid_old, nbset_old
    real, dimension (-ntgrid:ntgrid) :: thetasave
    real, dimension (ntheta+1) :: thetaold, thetanew
    real, dimension (ntheta+1) :: bmagold, bmagnew
    integer :: nbmag
    integer :: i

    ntheta_old = ntheta
    ntgrid_old = ntgrid
    nbset_old = nbset

    thetasave = theta
    thetaold = theta(-ntheta/2:ntheta/2)
    bmagold = bmag(-ntheta/2:ntheta/2)

    call gridgen4_2 (ntheta_old+1,thetaold,bmagold, npadd, &
         alknob,epsknob,extrknob,thetamax,deltaw,widthw,tension, &
         ntheta,nbset,thetanew,bmagnew,bset)

    ! interpolate to new grid
    ntgrid = ntheta/2 + (nperiod-1)*ntheta

    theta(-ntheta/2:ntheta/2-1) = thetanew(1:ntheta)
    theta(ntheta/2) = theta(-ntheta/2)+2.*pi
    bmag(-ntheta/2:ntheta/2-1) = bmagnew(1:ntheta)
    bmag(ntheta/2) = bmag(-ntheta/2)

    do i = 1, nperiod-1
       theta(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) &
            = thetanew(1:ntheta) + real(2*i)*pi
       theta(ntheta/2+i*ntheta) = thetanew(1) + real(2*(i+1))*pi
       theta(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) &
            = thetanew(1:ntheta) - real(2*i)*pi
       bmag(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) = bmagnew(1:ntheta)
       bmag(ntheta/2+i*ntheta) = bmagnew(1)
       bmag(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) = bmagnew(1:ntheta)
    end do

    call regrid (ntgrid_old, thetasave, gradpar, ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, gbdrift, ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, gbdrift0, ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, cvdrift, ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, cvdrift0, ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, gds2,  ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, gds21, ntgrid, theta(-ntgrid:ntgrid))
    call regrid (ntgrid_old, thetasave, gds22, ntgrid, theta(-ntgrid:ntgrid))

  end subroutine gridgen_get_grids

  subroutine regrid (nold, x, y, nnew, xnew)
    use splines
    implicit none
    integer, intent (in) :: nold
    real, dimension (-nold:nold), intent (in) :: x
    real, dimension (-nold:nold), intent (in out) :: y
    integer, intent (in) :: nnew
    real, dimension (-nnew:nnew), intent (in) :: xnew
    type (spline) :: spl
    integer :: i

    call new_spline (2*nold+1, x(-nold:nold), y(-nold:nold), spl)

    do i = -nnew, nnew
       y(i) = splint(xnew(i), spl)
    end do

    call delete_spline (spl)
  end subroutine regrid

end module theta_grid_gridgen
