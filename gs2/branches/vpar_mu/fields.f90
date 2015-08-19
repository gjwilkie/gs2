
module fields
  use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
  use fields_arrays, only: phip, aparp, bparp, phipnew, aparpnew, bparpnew
  use fields_arrays, only: phih, aparh, bparh, phihnew, aparhnew, bparhnew
  use fields_arrays, only: phitmp, apartmp, bpartmp
  use fields_arrays, only: phitmp1, apartmp1, bpartmp1
  use fields_arrays, only: phi_ext, apar_ext

  implicit none

  public :: init_fields, finish_fields
  public :: read_parameters, wnml_fields, check_fields
  public :: advance
  public :: phinorm, kperp
  public :: phi, apar, bpar, phinew, aparnew, bparnew
  public :: phip, aparp, bparp, phipnew, aparpnew, bparpnew
  public :: phih, aparh, bparh, phihnew, aparhnew, bparhnew
  public :: reset_init

  private

  ! knobs
  integer :: fieldopt_switch
  logical :: remove_zonal_flows_switch
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2

  logical :: initialized = .false.
  logical :: exist

contains

  subroutine check_fields(report_unit)
  implicit none
  integer :: report_unit
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly.')")
    case (fieldopt_test)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The field equations will only be tested.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end select
  end subroutine check_fields


  subroutine wnml_fields(unit)
  implicit none
  integer :: unit
     if (.not. exist) return 
       write (unit, *)
       write (unit, fmt="(' &',a)") "fields_knobs"
       select case (fieldopt_switch)
       case (fieldopt_implicit)
          write (unit, fmt="(' field_option = ',a)") '"implicit"'
       case (fieldopt_test)
          write (unit, fmt="(' field_option = ',a)") '"test"'
       end select
       write (unit, fmt="(' /')")
  end subroutine wnml_fields

  subroutine init_fields
!CMR,18/2/2011:
! add optional logical arg "noalloc" to avoid array allocations in ingen
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn, write_mpdist
    use init_g, only: ginit, init_init_g
    use fields_implicit, only: init_fields_implicit, init_phi_implicit
    use fields_test, only: init_fields_test, init_phi_test
    use nonlinear_terms, only: nl_finish_init => finish_init
    use antenna, only: init_antenna
    implicit none
    logical :: restarted
    logical:: debug=.false.

    if (initialized) return
    initialized = .true.
    
    if (debug) write(6,*) "init_fields: init_theta_grid"
    call init_theta_grid
    
!CMR,30/3/2009:
! call init_init_g before init_run_parameters to read delt from restart file

    if (debug) write(6,*) "init_fields: init_init_g"
    call init_init_g
    if (debug) write(6,*) "init_fields: init_run_parameters"
    call init_run_parameters
    if (debug) write(6,*) "init_fields: init_dist_fn"
    call init_dist_fn
    !if (debug) write(6,*) "init_fields: init_parameter_scan"
    !call init_parameter_scan
    if (debug) write(6,*) "init_fields: read_parameters"
    call read_parameters
    if (debug) write(6,*) "init_fields: allocate_arrays"
    call allocate_arrays

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       if (debug) write(6,*) "init_fields: init_fields_implicit"
       call init_fields_implicit
    case (fieldopt_test)
       if (debug) write(6,*) "init_fields: init_fields_test"
       call init_fields_test
    end select

! Turn on nonlinear terms.
    if (debug) write(6,*) "init_fields: nl_finish_init"
    call nl_finish_init

    if (debug) write(6,*) "init_fields: ginit"
    call ginit (restarted)

    if (debug) write(6,*) "init_fields: init_antenna"
    call init_antenna
    if (restarted) return

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       if (debug) write(6,*) "init_fields: init_phi_implicit"
       call init_phi_implicit
    case (fieldopt_test)
       if (debug) write(6,*) "init_fields: init_phi_test"
       call init_phi_test
    end select
    
  end subroutine init_fields

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (3), parameter :: fieldopts = &
         (/ text_option('default', fieldopt_implicit), &
            text_option('implicit', fieldopt_implicit), &
            text_option('test', fieldopt_test) /)
    character(20) :: field_option
    namelist /fields_knobs/ field_option, remove_zonal_flows_switch
    integer :: ierr, in_file

    if (proc0) then
       field_option = 'default'
       remove_zonal_flows_switch = .false.

       in_file = input_unit_exist ("fields_knobs", exist)
!       if (exist) read (unit=input_unit("fields_knobs"), nml=fields_knobs)
       if (exist) read (unit=in_file, nml=fields_knobs)

       ierr = error_unit()
       call get_option_value &
            (field_option, fieldopts, fieldopt_switch, &
            ierr, "field_option in fields_knobs")

    end if

    call broadcast (fieldopt_switch)
    call broadcast (remove_zonal_flows_switch)
  end subroutine read_parameters

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    implicit none
!    logical :: alloc = .true.

!    if (alloc) then
    if (.not. allocated(phi)) then
       allocate (     phi (-ntgrid:ntgrid,ntheta0,naky))
       allocate (    apar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (   bpar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (  phinew (-ntgrid:ntgrid,ntheta0,naky))
       allocate ( aparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (phip (-ntgrid:ntgrid,ntheta0,naky))
       allocate (aparp (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparp (-ntgrid:ntgrid,ntheta0,naky))
       allocate (phipnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (aparpnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparpnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (phih (-ntgrid:ntgrid,ntheta0,naky))
       allocate (aparh (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparh (-ntgrid:ntgrid,ntheta0,naky))
       allocate (phihnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (aparhnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparhnew (-ntgrid:ntgrid,ntheta0,naky))
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
    phip = 0. ; phipnew = 0.
    aparp = 0. ; aparpnew = 0.
    bparp = 0. ; bparpnew = 0.
    phih = 0. ; phihnew = 0.
    aparh = 0. ; aparhnew = 0.
    bparh = 0. ; bparhnew = 0.
!    phitmp1 = 0. ; apartmp1 = 0. ; bpartmp1 = 0.
!    phi_ext = 0.
    apar_ext = 0.

!    alloc = .false.
  end subroutine allocate_arrays

  subroutine advance (istep)
    use fields_implicit, only: advance_implicit
    use fields_test, only: advance_test
    implicit none
    integer, intent (in) :: istep

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call advance_implicit (istep, remove_zonal_flows_switch)
    case (fieldopt_test)
       call advance_test (istep)
    end select
  end subroutine advance

  subroutine phinorm (phifnc, aparfnc, bparfnc, phitot)
    use theta_grid, only: delthet, ntgrid
    use kt_grids, only: naky, ntheta0
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc, aparfnc, bparfnc
    real, dimension (:,:), intent (out) :: phitot
    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          phitot(it,ik) = 0.5/pi &
           *(sum((abs(phifnc(:,it,ik))**2 + abs(aparfnc(:,it,ik))**2 &
                  + abs(bparfnc(:,it,ik))**2) &
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
          if (anorm < 2.0*epsilon(0.0) .or. aky(ik) < epsilon(0.0)) then
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

  !!> This generates a flux surface average of phi. 

  !subroutine flux_surface_average_phi (phi_in, phi_average)
    !use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    !use kt_grids, only: ntheta0, naky

    !implicit none
    !complex, intent (in) :: phi_in
    !complex, intent (out) :: phi_average
    !complex, dimension(-ntgrid:ntgrid,1:ntheta0,1:naky) :: phi_fieldline_avg
    !integer it, ig

    !call fieldline_average_phi(phi_in, phi_fieldline_avg)
    !do it = 1,ntheta0
      !do ig = -ntgrid,ntgrid
        !phi_average(ig, it, :) = sum(phi_fieldline_avg(ig, it, :))/real(naky)
      !end do
    !end do

  !end subroutine fieldline_average_phi


  subroutine reset_init

    initialized  = .false.
    phi = 0.
    phinew = 0.
    phip = 0.
    phipnew = 0.
    phih = 0.
    phihnew = 0.

  end subroutine reset_init

  ! subroutine timer
    
  !   character (len=10) :: zdate, ztime, zzone
  !   integer, dimension(8) :: ival
  !   real, save :: told=0., tnew=0.
    
  !   call date_and_time (zdate, ztime, zzone, ival)
  !   tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
  !   if (told > 0.) then
  !      print *, 'Fields: Time since last called: ',tnew-told,' seconds'
  !   end if
  !   told = tnew
  ! end subroutine timer

  subroutine finish_fields

    use fields_implicit, only: implicit_reset => reset_init
    use fields_test, only: test_reset => reset_init
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phitmp, apartmp, bpartmp, apar_ext
    use fields_arrays, only: phip, aparp, bparp, phipnew, aparpnew, bparpnew
    use fields_arrays, only: phih, aparh, bparh, phihnew, aparhnew, bparhnew

    implicit none

    call reset_init
    call implicit_reset
    call test_reset

    if (allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew, &
         phitmp, apartmp, bpartmp, apar_ext)
    if (allocated(phip)) deallocate (phip, aparp, bparp, phipnew, aparpnew, bparpnew)
    if (allocated(phih)) deallocate (phih, aparh, bparh, phihnew, aparhnew, bparhnew)

  end subroutine finish_fields

end module fields
