module prof
  implicit none

  public :: init_prof
  public :: prof_entering
  public :: prof_leaving

  private

  integer, parameter :: maxregions = 20

  character (len=80), dimension (maxregions,2) :: regions
  integer :: nregions

contains

  subroutine init_prof
    implicit none
    nregions = 0
  end subroutine init_prof

  subroutine find_region (name, group, iregion)
    implicit none
    character(*), intent (in) :: name, group
    integer, intent (out) :: iregion
    integer :: i, ierr

    do i = 1, nregions
       if (name == trim(regions(i,1)) .and. group == trim(regions(i,2))) then
          iregion = i
          return
       end if
    end do
    if (nregions >= maxregions) then
       iregion = 0
       return
    end if
    nregions = nregions + 1
    regions(nregions,1) = name
    regions(nregions,2) = group
    iregion = nregions
    call vtsymdef (iregion, name, group, ierr)
    return
  end subroutine find_region

  subroutine prof_entering (name, group)
    implicit none
    character(*), intent (in) :: name
    character(*), intent (in), optional :: group
    integer :: i, ierr

    if (present(group)) then
       call find_region (name, group, i)
    else
       call find_region (name, "user", i)
    end if
    if (i == 0) return
    call vtbegin (i, ierr)
  end subroutine prof_entering

  subroutine prof_leaving (name, group)
    implicit none
    character(*), intent (in) :: name
    character(*), intent (in), optional :: group
    integer :: i, ierr

    if (present(group)) then
       call find_region (name, group, i)
    else
       call find_region (name, "user", i)
    end if
    if (i == 0) return
    call vtend (i, ierr)
  end subroutine prof_leaving

end module prof
