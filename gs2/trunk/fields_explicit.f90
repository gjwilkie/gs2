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
    use fields_arrays, only: phitmp, apartmp, aperptmp
    use fields_arrays, only: phitmp1, apartmp1, aperptmp1
    use dist_fn, only: timeadv, getfieldexp
    use dist_fn_arrays, only: g, gnew, gold
    implicit none
    integer, intent (in) :: istep
    real :: dt_cfl

!    gold = g ; phitmp1 = phi ; apartmp1 = apar ; aperptmp1 = aperp

    g = gnew ; phi = phinew ; apar = aparnew ; aperp = aperpnew

    call getfieldexp (phinew, aparnew, aperpnew)

!    phitmp = (8.*phinew - 7.*phi + 2.*phitmp1)/3.
!    apartmp = (1.5*aparnew - 0.5*apar)
!    aperptmp = (1.5*aperpnew - 0.5*aperp)

    phitmp = (1.5*phinew - 0.5*phi)
!    apartmp = (1.5*aparnew - 0.5*apar)
!    aperptmp = (1.5*aperpnew - 0.5*aperp)

    call timeadv (phitmp, apartmp, aperptmp, phitmp, apartmp, aperptmp, istep, dt_cfl, -1)

  end subroutine advance_explicit

  subroutine reset_init

    initialized = .false.
    
  end subroutine reset_init

end module fields_explicit
