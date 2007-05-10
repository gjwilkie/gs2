module fields_test
  implicit none

  public :: init_fields_test
  public :: advance_test
  public :: init_phi_test

  private

contains

  subroutine init_fields_test
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

  end subroutine init_fields_test

  subroutine init_phi_test
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    implicit none

    phi = 0.0
    apar = 0.0
    aperp = 0.0
    phinew = 0.0
    aparnew = 0.0
    aperpnew = 0.0
  end subroutine init_phi_test

  subroutine advance_test (istep)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use dist_fn, only: timeadv
    use dist_fn_arrays, only: g, gnew
    implicit none
    integer, intent (in) :: istep

    g = gnew
    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep)
  end subroutine advance_test

end module fields_test
