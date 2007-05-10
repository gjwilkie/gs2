module additional_terms
  implicit none

  public :: init_additional_terms
  public :: reset_init

  private

  logical :: phi0_term
  logical :: initialized = .false.

contains

  subroutine init_additional_terms
    use mp, only: proc0
    implicit none

    if (initialized) return
    initialized = .true.

    if (proc0) call read_parameters
    call broadcast_parameters

    if (phi0_term) call init_phi0_term
  end subroutine init_additional_terms

  subroutine read_parameters
    use file_utils, only: input_unit
    implicit none
    namelist /additional_terms_knobs/ phi0_term

    phi0_term = .false.
    read (unit=input_unit("additional_terms_knobs"), &
         nml=additional_terms_knobs)
  end subroutine read_parameters

  subroutine broadcast_parameters
    use mp, only: broadcast
    implicit none

    call broadcast (phi0_term)
  end subroutine broadcast_parameters

  subroutine reset_init

    initialized = .false.

  end subroutine reset_init

end module additional_terms
