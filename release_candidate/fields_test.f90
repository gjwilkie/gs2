module fields_test
  implicit none

  private

  public :: init_fields_test
  public :: advance_test
  public :: init_phi_test
  public :: reset_init

  logical :: initialized = .false.

contains

  subroutine init_fields_test
    implicit none

    if (initialized) return
    initialized = .true.

  end subroutine init_fields_test

  subroutine init_phi_test
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    implicit none

    phi = 0.0
    apar = 0.0
    bpar = 0.0
    phinew = 0.0
    aparnew = 0.0
    bparnew = 0.0
  end subroutine init_phi_test

  subroutine advance_test (istep)
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn, only: timeadv
    use dist_fn_arrays, only: g, gnew
    implicit none
    integer, intent (in) :: istep

    g = gnew
    call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep)
  end subroutine advance_test

  subroutine reset_init
    initialized = .false.
  end subroutine reset_init

end module fields_test
