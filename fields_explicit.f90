module fields_explicit

! NOT UP TO DATE ... DO NOT USE

  implicit none

  public :: init_fields_explicit
  public :: advance_explicit
  public :: init_phi_explicit
  public :: reset_init

  private

  logical :: initialized = .false.

contains

  subroutine init_fields_explicit
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    call init_kt_grids
    
  end subroutine init_fields_explicit

  subroutine init_phi_explicit
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use dist_fn, only: getfieldexp

    call getfieldexp (phinew, aparnew, aperpnew)
    phi = phinew; apar = aparnew; aperp = aperpnew
  end subroutine init_phi_explicit

  subroutine advance_explicit (istep, dt_cfl)
    use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew
    use dist_fn, only: timeadv, getfieldexp
    use dist_fn_arrays, only: g, gnew
    implicit none
    integer, intent (in) :: istep
    real :: dt_cfl

    phi = phinew; apar = aparnew; aperp = aperpnew
    g = gnew
    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep, dt_cfl, -1)
    call getfieldexp (phinew, aparnew, aperpnew)
    phi = phinew
!    call timeadv (phi, apar, aperp, phinew, aparnew, aperpnew, istep, dt_cfl)
!    call getfieldexp (phinew, aparnew, aperpnew)
  end subroutine advance_explicit

  subroutine reset_init

    initialized = .false.
    
  end subroutine reset_init

end module fields_explicit
