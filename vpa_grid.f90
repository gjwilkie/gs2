module vpa_grid

  implicit none
  
  public :: init_vpa_grid, finish_vpa_grid
  public :: vpa, wgts_vpa, vpac, dvpa

  real, dimension (:), allocatable :: vpa, wgts_vpa, vpac, dvpa

contains

  subroutine init_vpa_grid

    use gs3_input, only: read_vpa_grid_parameters, nvgrid, vpa_max, vpa_impfac
    use constants, only: pi => dpi
    use interp, only: get_cell_value

    implicit none

    integer :: iv, i

    call read_vpa_grid_parameters

    if (.not. allocated(vpa)) then
       ! vpa is the parallel velocity at grid points
       allocate (vpa(-nvgrid:nvgrid)) ; vpa = 0.0
       ! wgts_vpa are the integration weights assigned
       ! to the parallel velocity grid points
       allocate (wgts_vpa(-nvgrid:nvgrid)) ; wgts_vpa = 0.0
       ! vpac is the parallel velocity at cells
       allocate (vpac(-nvgrid:nvgrid)) ; vpac = 0.0
       ! dvpa is the grid spacing in vpa
       allocate (dvpa(-nvgrid:nvgrid)) ; dvpa = 0.0
    end if

    ! velocity grid goes from -vpa_max to vpa_max
    ! with a point at vpa = 0

    ! obtain vpa grid for vpa >= 0
    do iv = 0, nvgrid
       vpa(iv) = real(iv)*vpa_max/nvgrid
    end do
    ! fill in vpa grid for vpa < 0
    vpa(-nvgrid:-1) = -vpa(nvgrid:1:-1)

    dvpa(-nvgrid:nvgrid-1) = (/ (vpa(i+1)-vpa(i), i=-nvgrid,nvgrid-1) /)
    ! dvpa(nvgrid) should never be needed, but give it a nonzero
    ! value so division by zero never occurs when doing matrix operations
    dvpa(nvgrid) = dvpa(nvgrid-1)

    ! get vpa at cell centers
    call get_cell_value (vpa_impfac, vpa, vpac)
    ! vpac(nvgrid) should not be needed but set it nonzero to avoid possible
    ! divide by zero
    vpac(nvgrid) = vpa(nvgrid)

    ! get integration weights corresponding to vpa grid points
    wgts_vpa = 2.*vpa_max/(2.*nvgrid+1)

  end subroutine init_vpa_grid

  subroutine finish_vpa_grid

    implicit none

    if (allocated(vpa)) deallocate (vpa)
    if (allocated(wgts_vpa)) deallocate (wgts_vpa)
    if (allocated(vpac)) deallocate (vpac)
    if (allocated(dvpa)) deallocate (dvpa)

  end subroutine finish_vpa_grid

end module vpa_grid
