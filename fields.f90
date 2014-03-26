module fields
  use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
  use fields_arrays, only: phitmp, apartmp, bpartmp
  use fields_arrays, only: phitmp1, apartmp1, bpartmp1
  use fields_arrays, only: phi_ext, apar_ext

  implicit none

  public :: init_fields, finish_fields
  public :: read_parameters, wnml_fields, check_fields
  public :: advance
  public :: phinorm, kperp, fieldlineavgphi
  public :: phi, apar, bpar, phinew, aparnew, bparnew
  public :: reset_init, set_init_fields

  !> Made public for unit tests
  public :: fields_pre_init
  public :: remove_zonal_flows_switch

  private

  interface fieldlineavgphi
     module procedure fieldlineavgphi_loc
     module procedure fieldlineavgphi_tot
  end interface

  ! knobs
  integer :: fieldopt_switch
  logical :: remove_zonal_flows_switch
  logical :: force_maxwell_reinit
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2, fieldopt_local = 3
  logical :: dump_response, read_response
  logical :: initialized = .false.
  logical :: exist

contains

  subroutine check_fields(report_unit)
    use fields_local, only: minNrow, do_smart_update
    implicit none
    integer :: report_unit
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly.')")
       if(dump_response) write (report_unit, fmt="('The response matrix will be dumped to file.')")
       if(read_response) write (report_unit, fmt="('The response matrix will be read from file.')")
    case (fieldopt_test)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The field equations will only be tested.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    case (fieldopt_local)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly with decomposition respecting g_lo layout.')")
       if(dump_response) write (report_unit, fmt="('The response matrix will be dumped to file.')")
       if(read_response) write (report_unit, fmt="('The response matrix will be read from file.')")
       write(report_unit, fmt="('Using a min block size of ',I0)") minNrow
       if(do_smart_update) write(report_unit, fmt="('Using optimised field update.')")
    end select
  end subroutine check_fields


  subroutine wnml_fields(unit)
    use fields_local, only: minNrow, do_smart_update
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
       case (fieldopt_local)
          write (unit, fmt="(' field_option = ',a)") '"local"'
          write (unit, fmt="(' minNrow = ',I0)") minNrow
          write (unit, fmt="(' do_smart_update = ',L1)") do_smart_update
       end select
       if(dump_response) write (unit, fmt="(' dump_response = ',L1)") dump_response
       if(read_response) write (unit, fmt="(' read_response = ',L1)") read_response
       write (unit, fmt="(' /')")
  end subroutine wnml_fields

  !> Calls all initialisations required for init_fields_implicit/local, 
  !! reads parameters and allocates field arrays
  subroutine fields_pre_init
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn
    use init_g, only: ginit, init_init_g
    implicit none
    logical, parameter :: debug=.false.

    
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
  end subroutine fields_pre_init

  subroutine init_fields
!CMR,18/2/2011:
! add optional logical arg "noalloc" to avoid array allocations in ingen
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn
    use init_g, only: ginit, init_init_g
    use fields_implicit, only: init_fields_implicit
    use fields_test, only: init_fields_test
    use fields_local, only: init_fields_local
    use nonlinear_terms, only: nl_finish_init => finish_init
    use antenna, only: init_antenna
    use kt_grids, only: gridopt_switch, gridopt_box, kwork_filter
    use mp, only: iproc
    implicit none
    logical :: restarted
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    call fields_pre_init

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       if (debug) write(6,*) "init_fields: init_fields_implicit"
       call init_fields_implicit
    case (fieldopt_test)
       if (debug) write(6,*) "init_fields: init_fields_test"
       call init_fields_test
    case (fieldopt_local)
       if (debug) write(6,*) "init_fields: init_fields_local"
       call init_fields_local
    end select

! Turn on nonlinear terms.
    if (debug) write(6,*) "init_fields: nl_finish_init"
    call nl_finish_init

    if (debug) write(6,*) "init_fields: ginit"
    call ginit (restarted)
    if (debug) write(6,*) "init_fields: init_antenna"
    call init_antenna
    if (restarted .and. .not. force_maxwell_reinit) return

    !Set the initial fields
    call set_init_fields

    !If running in flux tube disable evolution of ky=kx=0 mode
    if(gridopt_switch.eq.gridopt_box) kwork_filter(1,1)=.true.
  end subroutine init_fields

  subroutine set_init_fields
    use fields_implicit, only: init_allfields_implicit
    use fields_test, only: init_phi_test
    use mp, only: proc0
    use fields_local, only: init_allfields_local
    implicit none
    logical, parameter :: debug=.false.
    if(proc0.and.debug) write(6,*) "Syncing fields with g."
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       if (debug) write(6,*) "init_fields: init_allfields_implicit"
       call init_allfields_implicit
    case (fieldopt_test)
       if (debug) write(6,*) "init_fields: init_phi_test"
       call init_phi_test
    case (fieldopt_local)
       if (debug) write(6,*) "init_fields: init_allfields_local"
       call init_allfields_local
    end select
  end subroutine set_init_fields

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    use fields_implicit, only: field_subgath, dump_response_imp=>dump_response, read_response_imp=>read_response
    use fields_local, only: dump_response_local=>dump_response, read_response_local=>read_response, minNrow
    use fields_local, only: do_smart_update
    implicit none
    type (text_option), dimension (5), parameter :: fieldopts = &
         (/ text_option('default', fieldopt_implicit), &
            text_option('implicit', fieldopt_implicit), &
            text_option('test', fieldopt_test),&
            text_option('local', fieldopt_local),&
            text_option('implicit_local', fieldopt_local)/)
    character(20) :: field_option
    namelist /fields_knobs/ field_option, remove_zonal_flows_switch, field_subgath, force_maxwell_reinit,&
         dump_response, read_response, minNrow, do_smart_update
    integer :: ierr, in_file

    if (proc0) then
       field_option = 'default'
       remove_zonal_flows_switch = .false.
       field_subgath=.false.
       force_maxwell_reinit=.true.
       dump_response=.false.
       read_response=.false.
       minnrow=64 !Tuning this can influence both init and advance times
       do_smart_update=.false.
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
    call broadcast (field_subgath)
    call broadcast (force_maxwell_reinit)
    call broadcast (dump_response)
    call broadcast (read_response)

    !Set the solve type specific flags
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       dump_response_imp=dump_response
       read_response_imp=read_response
    case (fieldopt_test)
    case (fieldopt_local)
       call broadcast (minNrow)
       call broadcast (do_smart_update)
       dump_response_local=dump_response
       read_response_local=read_response
    end select

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

!    alloc = .false.
  end subroutine allocate_arrays

  subroutine advance (istep)
    use fields_implicit, only: advance_implicit
    use fields_test, only: advance_test
    use fields_local, only: advance_local
    implicit none
    integer, intent (in) :: istep

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call advance_implicit (istep, remove_zonal_flows_switch)
    case (fieldopt_test)
       call advance_test (istep)
    case (fieldopt_local)
       call advance_local (istep, remove_zonal_flows_switch)
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
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag
!    use theta_grid, only: delthet
    implicit none
    complex, dimension (:,:), intent (out) :: phiavg
    real, dimension (-ntgrid:ntgrid) :: jac
!    integer :: ntg

!    ntg = ntgrid
    jac = 1.0/abs(drhodpsi*gradpar*bmag)
!    phiavg = sum(phi(-ntg:ntg,:,:)*jac(-ntg:ntg)*delthet(-ntg:ntg)) &
!         /sum(delthet(-ntg:ntg)*jac(-ntg:ntg))
    phiavg = 0.
    write(*,*) 'error in fields'
  end subroutine fieldlineavgphi_tot


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
    use fields_implicit, only: fi_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use fields_local, only: fl_reset => reset_fields_local
    implicit none
    initialized  = .false.
    phi = 0.
    phinew = 0.
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call fi_reset
    case (fieldopt_test)
       call ft_reset
    case (fieldopt_local)
       call fl_reset
    end select
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

  subroutine finish_fields

    use fields_implicit, only: implicit_reset => reset_init
    use fields_test, only: test_reset => reset_init
    use fields_local, only: finish_fields_local
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: phitmp, apartmp, bpartmp, apar_ext

    implicit none

    call reset_init
    
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call implicit_reset
    case (fieldopt_test)
       call test_reset
    case (fieldopt_local)
       call finish_fields_local
    end select

    if (allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew, &
         phitmp, apartmp, bpartmp, apar_ext)

  end subroutine finish_fields

end module fields
