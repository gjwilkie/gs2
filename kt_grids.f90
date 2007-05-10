module kt_grids
  implicit none

  public :: init_kt_grids
  public :: aky, theta0, akx, akr
  public :: naky, ntheta0

  private

  integer :: naky, ntheta0
  real, dimension (:), allocatable :: aky
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: akx
  real, dimension (:,:), allocatable :: akr

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4

contains

  subroutine init_kt_grids
    use theta_grid, only: init_theta_grid, shat, gds22
    use mp, only: proc0, broadcast
    implicit none
    logical, save :: initialized = .false.
    integer :: i

    if (initialized) return
    initialized = .true.

    call init_theta_grid

    if (proc0) then
       call read_parameters
       call get_sizes
    end if
    call broadcast (naky)
    call broadcast (ntheta0)
    call allocate_arrays
    if (proc0) call get_grids
    call broadcast (aky)
    call broadcast (akx)
    do i = 1, ntheta0
       call broadcast (theta0(:,i))
    end do

    do i = 1, ntheta0
       akr(:,i) = akx(i)*sqrt(abs(gds22))/abs(shat)
    end do
  end subroutine init_kt_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use text_options
    implicit none
    type (text_option), dimension (6), parameter :: gridopts = &
         (/ text_option('default', gridopt_single), &
            text_option('single', gridopt_single), &
            text_option('range', gridopt_range), &
            text_option('specified', gridopt_specified), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box) /)
    character(20) :: grid_option
    ! 'default' 'specified': specify grid in namelists

    namelist /kt_grids_knobs/ grid_option
    integer :: ierr

    grid_option = 'default'
    read (unit=input_unit("kt_grids_knobs"), nml=kt_grids_knobs)

    ierr = error_unit()
    call get_option_value &
         (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")
  end subroutine read_parameters

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    implicit none
    allocate (aky(naky))
    allocate (theta0(naky,ntheta0))
    allocate (akx(ntheta0))
    allocate (akr(-ntgrid:ntgrid,ntheta0))
  end subroutine allocate_arrays

  subroutine get_sizes
    use kt_grids_single, only: init_kt_grids_single, single_get_sizes
    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_specified, only: init_kt_grids_specified, specified_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call init_kt_grids_single
       call single_get_sizes (naky, ntheta0)
    case (gridopt_range)
       call init_kt_grids_range
       call range_get_sizes (naky, ntheta0)
    case (gridopt_specified)
       call init_kt_grids_specified
       call specified_get_sizes (naky, ntheta0)
    case (gridopt_box)
       call init_kt_grids_box
       call box_get_sizes (naky, ntheta0)
    end select
  end subroutine get_sizes

  subroutine get_grids
    use kt_grids_single, only: single_get_grids
    use kt_grids_range, only: range_get_grids
    use kt_grids_specified, only: specified_get_grids
    use kt_grids_box, only: box_get_grids
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call single_get_grids (aky, theta0, akx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx)
    case (gridopt_box)
       call box_get_grids (aky, theta0, akx)
    end select
  end subroutine get_grids

end module kt_grids
