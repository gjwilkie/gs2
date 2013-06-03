module vpamu_grids

  implicit none

  public :: init_vpa_grid, finish_vpa_grid
  public :: init_mu_grid, finish_mu_grid
  public :: vpa, nvgrid, wgts_vpa, dvpa
  public :: mu, nmu, wgts_mu

  integer :: nvgrid
  integer :: nmu

  real, dimension (:), allocatable :: mu, wgts_mu
  real, dimension (:), allocatable :: vpa, wgts_vpa, dvpa

contains

  subroutine init_vpamu_grids

    use mp, only: proc0

    implicit none

    if (proc0) then
       call read_parameters

       call init_vpa_grid
       call init_mu_grid
    end if

    call broadcast_results
    
  end subroutine init_vpamu_grids

  subroutine read_parameters

    use file_utils, only: input_unit_exist

    implicit none

    namelist /vpamu_grids_parameters/ nvgrid, nmu, vpa_max, mu_max

    integer :: in_file
    logical :: exist

    nvgrid = 8
    vpa_max = 3.0
    nmu = 16
    mu_max = 9.0

    in_file = input_unit_exist("vpamu_grids_parameters", exist)
    if (exist) read (unit=in_file, nml=vpamu_grids_parameters)

  end subroutine read_parameters

  subroutine init_vpa_grid

    implicit none

    integer :: iv, i

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

    ! get integration weights corresponding to vpa grid points
    wgts_vpa = 2.*vpa_max/(2.*nvgrid+1)

  end subroutine init_vpa_grid

  subroutine finish_vpa_grid

    implicit none

    if (allocated(vpa)) deallocate (vpa)
    if (allocated(wgts_vpa)) deallocate (wgts_vpa)
    if (allocated(dvpa)) deallocate (dvpa)

  end subroutine finish_vpa_grid

  subroutine init_mu_grid

    implicit none

    integer :: imu

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

  subroutine broadcast_results

    use mp, only: broadcast

    implicit none

    call broadcast (nvgrid)
    call broadcast (vpa)
    call broadcast (wgts_vpa)
    call broadcast (dvpa)
    call broadcast (nmu)
    call broadcast (mu)
    call broadcast (wgts_mu)

  end subroutine broadcast_results

end module vpamu_grids
