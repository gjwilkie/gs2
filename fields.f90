module fields
  use fields_arrays, only: phi, apar, aperp, phinew, aparnew, aperpnew

  implicit none

  public :: init_fields
  public :: advance
  public :: phinorm, kperp, fieldlineavgphi
  public :: phi, apar, aperp, phinew, aparnew, aperpnew
  public :: reset_init

  private

  ! knobs
  integer :: fieldopt_switch
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2

  logical :: initialized = .false.

contains

  subroutine init_fields
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters, quick
    use dist_fn, only: init_dist_fn, ginit
    use fields_implicit, only: init_fields_implicit, init_phi_implicit
    use fields_test, only: init_fields_test, init_phi_test
    implicit none
    logical :: restarted

    if (initialized) return
    initialized = .true.

    call init_theta_grid
! next two lines switched...2.8.98
    call init_run_parameters
    call init_kt_grids
    call init_dist_fn

    call read_parameters
    call allocate_arrays

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call init_fields_implicit (quick)
    case (fieldopt_test)
       call init_fields_test
    end select

    call ginit (restarted)
    if (restarted) return

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call init_phi_implicit
    case (fieldopt_test)
       call init_phi_test
    end select
  end subroutine init_fields

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use text_options
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (3), parameter :: fieldopts = &
         (/ text_option('default', fieldopt_implicit), &
            text_option('implicit', fieldopt_implicit), &
            text_option('test', fieldopt_test) /)
    character(20) :: field_option
    namelist /fields_knobs/ field_option
    integer :: ierr

    if (proc0) then
       field_option = 'default'
       
       read (unit=input_unit("fields_knobs"), nml=fields_knobs)

       ierr = error_unit()
       call get_option_value &
            (field_option, fieldopts, fieldopt_switch, &
            ierr, "field_option in fields_knobs")

    end if

    call broadcast (fieldopt_switch)
  end subroutine read_parameters

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
    logical :: alloc = .true.

    if (alloc) then
       allocate (phi(-ntgrid:ntgrid,naky,ntheta0))
       allocate (apar(-ntgrid:ntgrid,naky,ntheta0))
       allocate (aperp(-ntgrid:ntgrid,naky,ntheta0))
       allocate (phinew(-ntgrid:ntgrid,naky,ntheta0))
       allocate (aparnew(-ntgrid:ntgrid,naky,ntheta0))
       allocate (aperpnew(-ntgrid:ntgrid,naky,ntheta0))
    endif
    phi = 0.; phinew = 0.
    apar = 0.; aparnew = 0.
    aperp = 0.; aperpnew = 0.

    alloc = .false.
  end subroutine allocate_arrays

  subroutine advance (istep, dt_cfl)
    use fields_implicit, only: advance_implicit
    use fields_test, only: advance_test
    implicit none
    integer, intent (in) :: istep
    real :: dt_cfl

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call advance_implicit (istep, dt_cfl)
    case (fieldopt_test)
       call advance_test (istep)
    end select
  end subroutine advance

  subroutine phinorm (phitot)
    use theta_grid, only: ntgrid, ntheta, delthet
    use kt_grids, only: naky, ntheta0
    use constants
    implicit none
    real, dimension (:,:), intent (out) :: phitot
    integer :: ik, it

    do it = 1, ntheta0
       do ik = 1, naky
          phitot(ik,it) = 0.5/pi &
           *(sum((abs(phinew(:,ik,it))**2 + abs(aparnew(:,ik,it))**2 &
                  + abs(aperpnew(:,ik,it))**2) &
                 *delthet))
       end do
    end do
  end subroutine phinorm

  subroutine kperp (ntgrid_output, akperp)
    use theta_grid, only: delthet
    use kt_grids, only: naky, aky, ntheta0
    use dist_fn_arrays, only: kperp2
    implicit none
    integer, intent (in) :: ntgrid_output
    real, dimension (:,:), intent (out) :: akperp
    real :: anorm
    integer :: ik, it

    do it = 1, ntheta0
       do ik = 1, naky
          anorm = sum(abs(phinew(-ntgrid_output:ntgrid_output,ik,it) &
                         + aparnew(-ntgrid_output:ntgrid_output,ik,it) &
                         + aperpnew(-ntgrid_output:ntgrid_output,ik,it))**2 &
                      *delthet(-ntgrid_output:ntgrid_output))
          if (anorm < 2.0*epsilon(0.0) .or. aky(ik) == 0.0) then
             akperp(ik,it) = 0.0
          else
             akperp(ik,it) &
                  = sqrt(sum(kperp2(-ntgrid_output:ntgrid_output,ik,it) &
                     *abs(phinew(-ntgrid_output:ntgrid_output,ik,it) &
                          + aparnew(-ntgrid_output:ntgrid_output,ik,it) &
                          + aperpnew(-ntgrid_output:ntgrid_output,ik,it))**2 &
                     *delthet(-ntgrid_output:ntgrid_output))/anorm)
          end if
       end do
    end do
  end subroutine kperp

  subroutine fieldlineavgphi (ntgrid_output, ik, it, phiavg)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    implicit none
    integer, intent (in) :: ntgrid_output, ik, it
    complex, intent (out) :: phiavg
    real, dimension (-ntgrid:ntgrid) :: jac

    jac = 1.0/abs(drhodpsi*gradpar*bmag)
    phiavg = sum(phi(-ntgrid_output:ntgrid_output,ik,it) &
                 *jac(-ntgrid_output:ntgrid_output) &
                 *delthet(-ntgrid_output:ntgrid_output)) &
            /sum(delthet(-ntgrid_output:ntgrid_output) &
                 *jac(-ntgrid_output:ntgrid_output))
  end subroutine fieldlineavgphi

  subroutine reset_init

    initialized  = .false.

  end subroutine reset_init
end module fields
