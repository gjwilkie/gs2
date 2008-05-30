# include "define.inc"

! RN 2008/05/26:
!  give fpp -DVAMPIR to use Vampire analyzer
!  Vampire analyzer is now replaced by Intel Trace Collector and Analyzer,
!  but this may work with them.

module prof
  implicit none

  public :: init_prof
  public :: prof_entering
  public :: prof_leaving

  private

# ifdef VAMPIR
  integer, parameter :: maxregions = 20

  character (len=80), dimension (maxregions,2) :: regions
  integer :: nregions
# endif

contains

  subroutine init_prof
# ifdef VAMPIR
    implicit none
    nregions = 0
# endif
  end subroutine init_prof

# ifdef VAMPIR
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
#endif

  subroutine prof_entering (name, group)
    implicit none
    character(*), intent (in) :: name
    character(*), intent (in), optional :: group
# ifdef VAMPIR
    integer :: i, ierr

    if (present(group)) then
       call find_region (name, group, i)
    else
       call find_region (name, "user", i)
    end if
    if (i == 0) return
    call vtbegin (i, ierr)
#endif
  end subroutine prof_entering

  subroutine prof_leaving (name, group)
    implicit none
    character(*), intent (in) :: name
    character(*), intent (in), optional :: group
# ifdef VAMPIR
    integer :: i, ierr

    if (present(group)) then
       call find_region (name, group, i)
    else
       call find_region (name, "user", i)
    end if
    if (i == 0) return
    call vtend (i, ierr)
#endif
  end subroutine prof_leaving

end module prof
