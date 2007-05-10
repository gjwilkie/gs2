module prof
  implicit none

  public :: init_prof
  public :: prof_entering
  public :: prof_leaving

  private

contains

  subroutine init_prof
  end subroutine init_prof

  subroutine prof_entering (name, group)
    implicit none
    character(*), intent (in) :: name
    character(*), intent (in), optional :: group
  end subroutine prof_entering

  subroutine prof_leaving (name, group)
    implicit none
    character(*), intent (in) :: name
    character(*), intent (in), optional :: group
  end subroutine prof_leaving

end module prof
