module antenna_data

! quick fix for now; put amplitudes here and thus be 
! certain that there are no weird module dependencies
  implicit none

  private

  public :: init_antenna_data, finish_antenna_data
  public :: a_ant, b_ant, nk_stir, ant_on

  complex, dimension(:), allocatable :: a_ant, b_ant
  ! (nk_stir)
  integer :: nk_stir
  logical :: ant_on = .false., initialized = .false.

contains

  subroutine init_antenna_data (nk_stir_in)

    integer, intent (in) :: nk_stir_in

! do not reallocate this array on this processor...

    if (initialized) return
    initialized = .true.

    if (nk_stir_in <= 0) then
       ant_on = .false.
       return
    else
       ant_on = .true.
    end if

    nk_stir = nk_stir_in

    allocate (a_ant(nk_stir))
    allocate (b_ant(nk_stir))

  end subroutine init_antenna_data

  subroutine finish_antenna_data

    implicit none

    ant_on = .false.
!CMR,17/7/2009: add allocated test to avoid severe runtime library error
    if (allocated(a_ant) .and. allocated(b_ant)) deallocate (a_ant, b_ant)
    initialized = .false.

  end subroutine finish_antenna_data

end module antenna_data
