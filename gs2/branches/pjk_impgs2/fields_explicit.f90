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
!+PJK
    use fields_arrays, only: phiold, aparold, bparold
    use dist_fn_arrays, only: gold, g
!-PJK
    use dist_fn, only: getfieldexp

!+PJK
!    call getfieldexp (phinew, aparnew, bparnew)
!    phi = phinew; apar = aparnew; bpar = bparnew
!-PJK
    gold = g
    call getfieldexp (phiold, aparold, bparold)
    phi = phiold; apar = aparold; bpar = bparold

  end subroutine init_phi_explicit

  subroutine advance_explicit (istep)
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phitmp, apartmp, bpartmp
    use fields_arrays, only: phitmp1, apartmp1, bpartmp1
!+PJK
    use fields_arrays, only: phiold, aparold, bparold
!-PJK
    use dist_fn, only: timeadv, getfieldexp
    use dist_fn_arrays, only: g, gnew, gold
!+PJK
    use dist_fn_arrays, only: gwork
    use gs2_time, only: user2code, code_dt, user_dt
!-PJK
    implicit none
    integer, intent (in) :: istep
!+PJK
    real :: dtnew
!-PJK

!    gold = g ; phitmp1 = phi ; apartmp1 = apar ; bpartmp1 = bpar

!+PJK Following two lines are present in advance_implicit...
!    call antenna_amplitudes (apar_ext)
!    if (allocated(kx_shift)) call exb_shear (gnew, phinew, aparnew, bparnew) 
!-PJK

    !  Original code...
    g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew
    call getfieldexp (phinew, aparnew, bparnew)
    call timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, -1)

    !  New code...
goto 10
    adaptive_dt: do

       !  It is not clear where these four variables are used within the code,
       !  but they have to be set here, otherwise there is no convergence.
       g = gold ; phi = phiold ; apar = aparold ; bpar = bparold

       !  Given gnew as input, calculate fields
       call getfieldexp (phinew, aparnew, bparnew)
       !  Given old and new fields, calculate updated g, =gnew
       call timeadv (phiold, aparold, bparold, phinew, aparnew, bparnew, istep, -1)

       gwork = gnew ! to be compared with gnew after second timeadv to determine new dt

       !  Calculate fields using this 1st order guess for gnew
       call getfieldexp (phitmp, apartmp, bpartmp)

       !  Refine the guess for the fields and recalculate gnew
       phinew  = 0.5*(phinew + phitmp)
       aparnew = 0.5*(aparnew + apartmp)
       bparnew = 0.5*(bparnew + bpartmp)
       call timeadv (phiold, aparold, bparold, phinew, aparnew, bparnew, istep, -1)

       !  Modify the timestep according to how well the 1st order and 2nd order
       !  estimates for gnew agree

       dtnew = 0.9 * code_dt * sqrt( 1.0E-2 * maxval(abs(gnew)) / maxval(abs(gnew-gwork)) )

       if ( abs(dtnew/code_dt - 1.0) < 0.3 ) then
          code_dt = dtnew ! don't forget user_dt
          write(*,*) 'dt changed to ',code_dt
          exit adaptive_dt
       else
          code_dt = dtnew
          write(*,*) 'dt changed to ',code_dt,' - restarting timestep'
       end if

    end do adaptive_dt
10 continue
    !  For now, set gold = g, although this should only be done if we
    !  are happy with the convergence (so this is temporary coding)

    g = gnew
    gold = gnew ; phiold = phinew ; aparold = aparnew ; bparold = bparnew

  end subroutine advance_explicit

  subroutine reset_init

    initialized = .false.
    
  end subroutine reset_init

end module fields_explicit
