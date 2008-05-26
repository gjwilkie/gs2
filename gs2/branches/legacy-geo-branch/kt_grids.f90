module kt_grids_single
  implicit none

  public :: init_kt_grids_single, single_get_sizes, single_get_grids

  private

  real :: aky, theta0, akx

contains

  subroutine init_kt_grids_single
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist
    namelist /kt_grids_single_parameters/ aky, theta0, akx

    aky = 0.4
    theta0 = 0.0
    akx = 0.0
    in_file = input_unit_exist ("kt_grids_single_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_single_parameters"), &
         nml=kt_grids_single_parameters)
  end subroutine init_kt_grids_single

  subroutine single_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny
    naky = 1
    ntheta0 = 1
    nx = 0
    ny = 0
  end subroutine single_get_sizes

  subroutine single_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    call get_grids (aky, theta0, akx)
  end subroutine single_get_grids

  subroutine get_grids (aky_out, theta0_out, akx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out
    real, dimension (:,:), intent (out) :: theta0_out
    real, dimension (:), intent (out) :: akx_out
    aky_out(1) = aky
    theta0_out(1,1) = theta0
    akx_out(1) = akx
  end subroutine get_grids

end module kt_grids_single

module kt_grids_range
  implicit none

  public :: init_kt_grids_range, range_get_sizes, range_get_grids

  private

  integer :: naky_private, ntheta0_private
  real :: aky_min, aky_max, theta0_min, theta0_max

contains

  subroutine init_kt_grids_range
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: naky, ntheta0
    integer :: in_file
    logical :: exist
    namelist /kt_grids_range_parameters/ naky, ntheta0, &
         aky_min, aky_max, theta0_min, theta0_max

    naky = 1
    ntheta0 = 1
    aky_min = 0.0
    aky_max = 0.0
    theta0_min = 0.0
    theta0_max = 0.0
    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_range_parameters"), &
         nml=kt_grids_range_parameters)
    naky_private = naky
    ntheta0_private = ntheta0
  end subroutine init_kt_grids_range

  subroutine range_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = 0
    ny = 0
  end subroutine range_get_sizes

  subroutine range_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
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
       theta0(:,j) &
            = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0_private-1) /)
    end do

    akx = theta0(:,1) * shat * aky(1)

  end subroutine range_get_grids

end module kt_grids_range

module kt_grids_specified
  implicit none

  public :: init_kt_grids_specified, specified_get_sizes, specified_get_grids

  private

  integer :: naky_private, ntheta0_private, nx_private, ny_private

contains

  subroutine init_kt_grids_specified
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: naky, ntheta0, nx, ny
    integer :: in_file
    logical :: exist
    namelist /kt_grids_specified_parameters/ naky, ntheta0, nx, ny

    naky = 1
    ntheta0 = 1
    nx = 0
    ny = 0
    in_file = input_unit_exist("kt_grids_specified_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_specified_parameters"), &
         nml=kt_grids_specified_parameters)
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_specified

  subroutine specified_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny

    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
  end subroutine specified_get_sizes

  subroutine specified_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer :: i, naky, ntheta0
    real :: aky_dummy, theta0_dummy, akx_dummy

    naky = size(aky)
    ntheta0 = size(akx)
    do i = 1, max(naky,ntheta0)
       call read_element (i, aky_dummy, theta0_dummy, akx_dummy)
       if (i <= naky) aky(i) = aky_dummy
       if (i <= ntheta0) theta0(i,:) = theta0_dummy
       if (i <= ntheta0) akx(i) = akx_dummy
    end do
  end subroutine specified_get_grids

  subroutine read_element (i, aky_dummy, theta0_dummy, akx_dummy)
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer, intent (in) :: i
    real, intent (out) :: aky_dummy, theta0_dummy, akx_dummy
    integer :: unit
    real :: aky, theta0, akx
    namelist /kt_grids_specified_element/ aky, theta0, akx

    aky = 0.4
    theta0 = 0.0
    akx = 0.0
    call get_indexed_namelist_unit (unit, "kt_grids_specified_element", i)
    read (unit=unit, nml=kt_grids_specified_element)
    close (unit)
    aky_dummy = aky
    theta0_dummy = theta0
    akx_dummy = akx
  end subroutine read_element

end module kt_grids_specified

module kt_grids_box
  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids

  private

  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: nkpolar_private
  integer :: jtwist
  real :: ly, y0, x0, rtwist

contains

  subroutine init_kt_grids_box
    use theta_grid, only: init_theta_grid
    use file_utils, only: input_unit, input_unit_exist
    use constants
    implicit none
    integer :: naky, ntheta0, nx, ny, nkpolar
    integer :: in_file
    logical :: exist
    namelist /kt_grids_box_parameters/ naky, ntheta0, ly, nx, ny, jtwist, &
	y0, rtwist, x0, nkpolar

    call init_theta_grid

    nkpolar = 0
    naky = 0
    ntheta0 = 0
    ly = 0.0
    y0 = 2.0
    nx = 0
    ny = 0
    jtwist = 1
    rtwist = 0.0
    x0 = 0.
    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_box_parameters"), nml=kt_grids_box_parameters)
    if (y0 < 0) then
       y0 = -1./y0
    end if
    if (ly == 0.) ly = 2.0*pi*y0
    if (naky == 0) naky = (ny-1)/3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
    if (rtwist == 0.) rtwist = real(jtwist)
!    if (nkpolar == 0) nkpolar = naky  ! should be generalized later.  
                                       ! For now, assuming square domain with nx=ny
    if (nkpolar == 0) nkpolar = int(real(naky-1.)*sqrt(2.))
    
    nkpolar_private = nkpolar
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_box

  subroutine box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny, nkpolar
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
    nkpolar = nkpolar_private
  end subroutine box_get_sizes

!  subroutine box_get_grids (aky, theta0, akx, akpolar)
  subroutine box_get_grids (aky, theta0, akx, ikx, iky)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: aky
    real, dimension (:,:), intent (out) :: theta0
    real, dimension (:), intent (out) :: akx
    integer, dimension (:), intent (out) :: ikx, iky
!    real, dimension (:), intent (out) :: akpolar
    integer :: i, naky, ntheta0
    real :: dkx

    naky = size(aky)
    ntheta0 = size(akx)

    if(abs(shat) <=  1.e-5) then
       if (x0 == 0.) then
          if (rtwist > 0) then
             dkx = 2.0*pi/ly/rtwist
          else
             dkx = 2.0*pi/ly*rtwist
          endif
       else
          if (x0 > 0.) then
             dkx = 1./x0
          else
             dkx = -x0
          end if
       end if
    else
       dkx = 2.0*pi/real(jtwist)* 2.0*pi/ly*shat
    endif

    do i = 1, naky
       iky(i) = i-1
       aky(i) = real(i-1)*2.*pi/ly
    end do

    do i = 1, (ntheta0+1)/2
       ikx(i) = i-1
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       ikx(i) = i-ntheta0-1
       akx(i) = real(i-ntheta0-1)*dkx
    end do

!!! Is there a negative sign missing here?  For up-down symmetric
!!! equilibria, it would not matter, but otherwise it could??
!!! This is mainly used to define wdrift and kperp2
    do i = 1, ntheta0
       theta0(i,2:) = akx(i)/(aky(2:)*shat)
       theta0(i,1) = 0.0
    end do
    
! currently assuming square domain with nx=ny for MHD energy diagnostics:
!    akpolar = aky
    
  end subroutine box_get_grids

end module kt_grids_box

module kt_grids_xbox
  implicit none

  public :: init_kt_grids_xbox, xbox_get_sizes, xbox_get_grids

  private

  integer :: ntheta0_private, nx_private
  real :: lx, aky_private

contains

  subroutine init_kt_grids_xbox
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: ntheta0, nx
    integer :: in_file
    logical :: exist
    real :: aky
    namelist /kt_grids_xbox_parameters/ ntheta0, lx, aky, nx

    ntheta0 = 1
    lx = 1.0
    aky = 0.2
    nx = 0
    in_file = input_unit_exist ("kt_grids_xbox_parameters", exist)
    if (exist) read (unit=input_unit("kt_grids_xbox_parameters"), &
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
    integer :: i, ntheta0

    aky(1) = aky_private

    ntheta0 = size(akx)
    akx(:(ntheta0+1)/2) = (/ (real(2*(i-1))*pi/lx, i=1,(ntheta0+1)/2) /)
    akx((ntheta0+1)/2+1:) &
         = (/ (real(2*(i-ntheta0-1))*pi/lx, i=(ntheta0+1)/2+1,ntheta0) /)
    theta0(:,1) = akx(:)/(aky(1)*shat)
  end subroutine xbox_get_grids

end module kt_grids_xbox

module kt_grids
  implicit none

  public :: init_kt_grids, box
  public :: aky, theta0, akx, akr
  public :: aky_out, akx_out, akr_out
  public :: naky, ntheta0, nx, ny, reality
  public :: nkpolar 
  public :: ikx, iky
  private

  integer :: naky, ntheta0, nx, ny, nkpolar
  integer, dimension(:), allocatable :: ikx, iky
  real, dimension (:), allocatable :: aky, aky_out
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: akx, akx_out
  real, dimension (:,:), allocatable :: akr, akr_out
!  real, dimension (:), allocatable :: akpolar, akpolar_out

  ! internal variables
  integer :: gridopt_switch, normopt_switch
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4, gridopt_xbox = 5
  integer, parameter :: normopt_mtk = 1, normopt_bd = 2
  logical :: reality = .false.
  logical :: box = .false.

contains

  subroutine init_kt_grids (tnorm)
    use theta_grid, only: init_theta_grid, shat, gds22
    use mp, only: proc0, broadcast
    implicit none

    real, optional, intent (out) :: tnorm
    integer :: ik, it
    real :: tfac = 1.0
    logical, save :: initialized = .false.

    if (present(tnorm)) tnorm = tfac

    if (initialized) return
    initialized = .true.

    call init_theta_grid

    if (proc0) then
       nkpolar = 0   ! will be set to non-zero value only in box case; only used for an MHD diagnostic
       call read_parameters
       call get_sizes

       select case (normopt_switch)
       case (normopt_mtk)
          tfac = 1.
       case (normopt_bd)
          tfac = sqrt(2.)
       end select
    end if

    call broadcast (tfac)
    if (present(tnorm)) tnorm = tfac

    call broadcast (reality)
    call broadcast (box)
    call broadcast (naky)
    call broadcast (nkpolar)
    call broadcast (ntheta0)
    call broadcast (ny)
    call broadcast (nx)
    call allocate_arrays
    if (proc0) call get_grids
    call broadcast (aky)
    call broadcast (akx)
    call broadcast (ikx)
    do ik = 1, naky
       call broadcast (theta0(:,ik))
    end do
!    if (nkpolar > 0) call broadcast (akpolar)

    if (abs(shat) > epsilon(0.)) then
       do it = 1, ntheta0
          akr(:,it) = akx(it)*sqrt(abs(gds22))/abs(shat)
       end do
    else
       akr = 1.
    end if

    select case (normopt_switch)
    case (normopt_mtk)
       akr_out = akr
       akx_out = akx
       aky_out = aky
!       if (nkpolar > 0) akpolar_out = akpolar
    case (normopt_bd)
       akr_out = akr / sqrt(2.)
       akx_out = akx / sqrt(2.)
       aky_out = aky / sqrt(2.)
!       if (nkpolar > 0) akpolar_out = akpolar / sqrt(2.)
    end select

  end subroutine init_kt_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (7), parameter :: gridopts = &
         (/ text_option('default', gridopt_single), &
            text_option('single', gridopt_single), &
            text_option('range', gridopt_range), &
            text_option('specified', gridopt_specified), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box), &
            text_option('xbox', gridopt_xbox) /)
    character(20) :: grid_option
    ! 'default' 'specified': specify grid in namelists

    type (text_option), dimension(6), parameter :: normopts = &
         (/ text_option('default', normopt_mtk), &
            text_option('with_root_2', normopt_mtk), &
            text_option('mtk', normopt_mtk), &
            text_option('no_root_2', normopt_bd), &
            text_option('bd', normopt_bd), &
            text_option('t_over_m', normopt_bd) /)
    character(20) :: norm_option

    namelist /kt_grids_knobs/ grid_option, norm_option
    integer :: ierr, in_file
    logical :: exist

    norm_option = 'default'
    grid_option = 'default'
    in_file = input_unit_exist ("kt_grids_knobs", exist)
    if (exist) read (unit=input_unit("kt_grids_knobs"), nml=kt_grids_knobs)

    ierr = error_unit()
    call get_option_value &
         (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")

    ierr = error_unit()
    call get_option_value &
         (norm_option, normopts, normopt_switch, &
         ierr, "norm_option in kt_grids_knobs")

  end subroutine read_parameters

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    implicit none
    allocate (aky(naky))
    allocate (aky_out(naky))
    allocate (theta0(ntheta0,naky))
    allocate (akx(ntheta0))
    allocate (akx_out(ntheta0))
    allocate (akr(-ntgrid:ntgrid,ntheta0))
    allocate (akr_out(-ntgrid:ntgrid,ntheta0))

    allocate (ikx(ntheta0))
    allocate (iky(naky))

!    if (nkpolar > 0) &
!         allocate (akpolar(nkpolar), akpolar_out(nkpolar))
  end subroutine allocate_arrays

  subroutine get_sizes
    use kt_grids_single, only: init_kt_grids_single, single_get_sizes
    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_specified, only: init_kt_grids_specified, specified_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes
    use kt_grids_xbox, only: init_kt_grids_xbox, xbox_get_sizes
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call init_kt_grids_single
       call single_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_range)
       call init_kt_grids_range
       call range_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_specified)
       call init_kt_grids_specified
       call specified_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_box)
       call init_kt_grids_box
       call box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
       reality = .true.
       box = .true.
    case (gridopt_xbox)
       call init_kt_grids_xbox
       call xbox_get_sizes (naky, ntheta0, nx, ny)
    end select
  end subroutine get_sizes

  subroutine get_grids
    use kt_grids_single, only: single_get_grids
    use kt_grids_range, only: range_get_grids
    use kt_grids_specified, only: specified_get_grids
    use kt_grids_box, only: box_get_grids
    use kt_grids_xbox, only: xbox_get_grids
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call single_get_grids (aky, theta0, akx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx)
    case (gridopt_box)
!       call box_get_grids (aky, theta0, akx, akpolar)
       call box_get_grids (aky, theta0, akx, ikx, iky)
    case (gridopt_xbox)
       call xbox_get_grids (aky, theta0, akx)
    end select

    select case (normopt_switch)
    case (normopt_mtk)
       ! nothing -- this is how the original code is designed
    case (normopt_bd)
       aky = aky * sqrt(2.)
       akx = akx * sqrt(2.)
!       if (nkpolar > 0) akpolar = akpolar * sqrt(2.)
    end select

  end subroutine get_grids

end module kt_grids

