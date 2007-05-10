module kt_grids_specified
  implicit none

  public :: init_kt_grids_specified, specified_get_sizes, specified_get_grids

  private

  integer :: naky_private, ntheta0_private

contains

  subroutine init_kt_grids_specified
    use file_utils, only: input_unit
    implicit none
    integer :: naky, ntheta0
    namelist /kt_grids_specified_parameters/ naky, ntheta0

    naky = 1
    ntheta0 = 1
    read (unit=input_unit("kt_grids_specified_parameters"), &
         nml=kt_grids_specified_parameters)
    naky_private = naky
    ntheta0_private = ntheta0
  end subroutine init_kt_grids_specified

  subroutine specified_get_sizes (naky, ntheta0)
    implicit none
    integer, intent (out) :: naky, ntheta0

    naky = naky_private
    ntheta0 = ntheta0_private
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
       if (i <= ntheta0) theta0(:,i) = theta0_dummy
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
