module fields
  use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
  use fields_arrays, only: phitmp, apartmp, bpartmp
  use fields_arrays, only: phitmp1, apartmp1, bpartmp1
  use fields_arrays, only: phi_ext, apar_ext

  implicit none

  public :: init_fields
  public :: advance
  public :: phinorm, kperp, fieldlineavgphi
  public :: phi, apar, bpar, phinew, aparnew, bparnew
  public :: reset_init

  private

  interface fieldlineavgphi
     module procedure fieldlineavgphi_loc
     module procedure fieldlineavgphi_tot
  end interface

  ! knobs
  integer :: fieldopt_switch
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2, fieldopt_explicit = 3

  logical :: initialized = .false.

contains

  subroutine init_fields
    use mp, only: proc0
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn, ginit
    use init_g, only: init_init_g
    use fields_implicit, only: init_fields_implicit, init_phi_implicit
    use fields_explicit, only: init_fields_explicit, init_phi_explicit
    use fields_test, only: init_fields_test, init_phi_test
    use nonlinear_terms, only: nl_finish_init => finish_init
    use antenna, only: init_antenna
    implicit none
    logical :: restarted

    if (initialized) return
    initialized = .true.

    call init_theta_grid
    
    call init_init_g
    call init_run_parameters
    call init_dist_fn
    call read_parameters
    call allocate_arrays

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call init_fields_implicit
    case (fieldopt_explicit)
       call init_fields_explicit
    case (fieldopt_test)
       call init_fields_test
    end select

! Turn on nonlinear terms.
    call nl_finish_init

    call ginit (restarted)
    call init_antenna
    if (restarted) return

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call init_phi_implicit
    case (fieldopt_explicit)
       call init_phi_explicit
    case (fieldopt_test)
       call init_phi_test
    end select
    
  end subroutine init_fields

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (4), parameter :: fieldopts = &
         (/ text_option('default', fieldopt_implicit), &
            text_option('implicit', fieldopt_implicit), &
            text_option('explicit', fieldopt_explicit), &
            text_option('test', fieldopt_test) /)
    character(20) :: field_option
    namelist /fields_knobs/ field_option
    integer :: ierr, in_file
    logical :: exist

    if (proc0) then
       field_option = 'default'

       in_file = input_unit_exist ("fields_knobs", exist)
       if (exist) read (unit=input_unit("fields_knobs"), nml=fields_knobs)

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
       allocate (     phi (-ntgrid:ntgrid,ntheta0,naky))
       allocate (    apar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (   bpar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (  phinew (-ntgrid:ntgrid,ntheta0,naky))
       allocate ( aparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (  phitmp (-ntgrid:ntgrid,ntheta0,naky))
       allocate ( apartmp (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bpartmp (-ntgrid:ntgrid,ntheta0,naky))
!       allocate (  phitmp1(-ntgrid:ntgrid,ntheta0,naky))
!       allocate ( apartmp1(-ntgrid:ntgrid,ntheta0,naky))
!       allocate (bpartmp1(-ntgrid:ntgrid,ntheta0,naky))
!       allocate ( phi_ext (-ntgrid:ntgrid,ntheta0,naky))
       allocate (apar_ext (-ntgrid:ntgrid,ntheta0,naky))
    endif
    phi = 0.; phinew = 0.; phitmp = 0. 
    apar = 0.; aparnew = 0.; apartmp = 0. 
    bpar = 0.; bparnew = 0.; bpartmp = 0.
!    phitmp1 = 0. ; apartmp1 = 0. ; bpartmp1 = 0.
!    phi_ext = 0.
    apar_ext = 0.

    alloc = .false.
  end subroutine allocate_arrays

  subroutine advance (istep)
    use fields_implicit, only: advance_implicit
    use fields_explicit, only: advance_explicit
    use fields_test, only: advance_test
    implicit none
    integer, intent (in) :: istep

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call advance_implicit (istep)
    case (fieldopt_explicit)
       call advance_explicit (istep)
    case (fieldopt_test)
       call advance_test (istep)
    end select
  end subroutine advance

  subroutine phinorm (phitot)
    use theta_grid, only: delthet
    use kt_grids, only: naky, ntheta0
    use constants
    implicit none
    real, dimension (:,:), intent (out) :: phitot
    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          phitot(it,ik) = 0.5/pi &
           *(sum((abs(phinew(:,it,ik))**2 + abs(aparnew(:,it,ik))**2 &
                  + abs(bparnew(:,it,ik))**2) &
                 *delthet))
       end do
    end do
  end subroutine phinorm

  subroutine kperp (ntgrid_output, akperp)
    use theta_grid, only: delthet
    use kt_grids, only: naky, aky, ntheta0
    use run_parameters, only: fphi, fapar, fbpar
    use dist_fn_arrays, only: kperp2
    implicit none
    integer, intent (in) :: ntgrid_output
    real, dimension (:,:), intent (out) :: akperp
    real :: anorm
    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          anorm = sum(abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
                         + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
                         + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
                      *delthet(-ntgrid_output:ntgrid_output))
          if (anorm < 2.0*epsilon(0.0) .or. aky(ik) == 0.0) then
             akperp(it,ik) = 0.0
          else
             akperp(it,ik) &
                  = sqrt(sum(kperp2(-ntgrid_output:ntgrid_output,it,ik) &
                     *abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
                          + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
                          + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
                     *delthet(-ntgrid_output:ntgrid_output))/anorm)
          end if
       end do
    end do
  end subroutine kperp

  subroutine fieldlineavgphi_loc (ntgrid_output, it, ik, phiavg)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    implicit none
    integer, intent (in) :: ntgrid_output, ik, it
    complex, intent (out) :: phiavg
    real, dimension (-ntgrid:ntgrid) :: jac

    jac = 1.0/abs(drhodpsi*gradpar*bmag)
    phiavg = sum(phi(-ntgrid_output:ntgrid_output,it,ik) &
                 *jac(-ntgrid_output:ntgrid_output) &
                 *delthet(-ntgrid_output:ntgrid_output)) &
            /sum(delthet(-ntgrid_output:ntgrid_output) &
                 *jac(-ntgrid_output:ntgrid_output))
  end subroutine fieldlineavgphi_loc

  subroutine fieldlineavgphi_tot (phiavg)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    implicit none
    complex, dimension (:,:), intent (out) :: phiavg
    real, dimension (-ntgrid:ntgrid) :: jac
    integer :: ntg

    ntg = ntgrid
    jac = 1.0/abs(drhodpsi*gradpar*bmag)
!    phiavg = sum(phi(-ntg:ntg,:,:)*jac(-ntg:ntg)*delthet(-ntg:ntg)) &
!         /sum(delthet(-ntg:ntg)*jac(-ntg:ntg))
    phiavg = 0.
    write(*,*) 'error in fields'
  end subroutine fieldlineavgphi_tot

  subroutine reset_init

    initialized  = .false.
    phi = 0.
    phinew = 0.

  end subroutine reset_init

  subroutine timer
    
    character (len=10) :: zdate, ztime, zzone
    integer, dimension(8) :: ival
    real, save :: told=0., tnew=0.
    
    call date_and_time (zdate, ztime, zzone, ival)
    tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
    if (told > 0.) then
       print *, 'Fields: Time since last called: ',tnew-told,' seconds'
    end if
    told = tnew
  end subroutine timer

end module fields
