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
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn, only: getfieldexp

    call getfieldexp (phinew, aparnew, bparnew)
    phi = phinew; apar = aparnew; bpar = bparnew
  end subroutine init_phi_explicit

  subroutine advance_explicit (istep)
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use fields_arrays, only: phitmp1, apartmp1, bpartmp1
    use dist_fn, only: timeadv, getfieldexp
    use dist_fn_arrays, only: g, gnew, gold
    implicit none
    integer, intent (in) :: istep

!    gold = g ; phitmp1 = phi ; apartmp1 = apar ; bpartmp1 = bpar

    g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew

    call getfieldexp (phinew, aparnew, bparnew)

!    phitmp = (8.*phinew - 7.*phi + 2.*phitmp1)/3.
!    apartmp = (1.5*aparnew - 0.5*apar)
!    bpartmp = (1.5*bparnew - 0.5*bpar)

    phitmp = (1.5*phinew - 0.5*phi)
!    apartmp = (1.5*aparnew - 0.5*apar)
!    bpartmp = (1.5*bparnew - 0.5*bpar)

    call timeadv (phitmp, apartmp, bpartmp, phitmp, apartmp, bpartmp, istep, -1)

  end subroutine advance_explicit

  subroutine reset_init

    initialized = .false.
    
  end subroutine reset_init

end module fields_explicit
