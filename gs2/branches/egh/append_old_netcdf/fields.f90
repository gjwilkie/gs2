module fields

  implicit none

  private 

  public :: init_fields, finish_fields
  public :: read_parameters, wnml_fields, check_fields
  public :: advance, force_maxwell_reinit
  public :: reset_init, set_init_fields
  public :: fields_init_response, set_dump_and_read_response
  public :: dump_response_to_file

  !> Made public for unit tests
  public :: fields_pre_init
  public :: remove_zonal_flows_switch
  !> Made public for replay
  public :: allocate_arrays

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
    integer, intent(in) :: report_unit
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
    integer, intent(in) :: unit
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
    use antenna, only: init_antenna
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    
    call debug_message(verb, "init_fields: init_theta_grid")
    call init_theta_grid
    
!CMR,30/3/2009:
! call init_init_g before init_run_parameters to read delt from restart file

    call debug_message(verb, "init_fields: init_init_g")
    call init_init_g
    call debug_message(verb, "init_fields: init_run_parameters")
    call init_run_parameters
    call debug_message(verb, "init_fields: init_dist_fn")
    call init_dist_fn
    !call debug_message(verb, "init_fields: init_parameter_scan")
    !call init_parameter_scan
    call debug_message(verb, "init_fields: init_antenna")
    call init_antenna !Must come before allocate_arrays so we know if we need apar_ext
    call debug_message(verb, "init_fields: read_parameters")
    call read_parameters
    call debug_message(verb, "init_fields: allocate_arrays")
    call allocate_arrays
  end subroutine fields_pre_init

  subroutine fields_init_response
    use fields_implicit, only: init_fields_implicit
    use fields_test, only: init_fields_test
    use fields_local, only: init_fields_local
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    logical, parameter :: debug = .false.
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_implicit")
       call init_fields_implicit
    case (fieldopt_test)
       call debug_message(verb, "fields::fields_init_response init_fields_test")
       call init_fields_test
    case (fieldopt_local)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_local")
       call init_fields_local
    case default
       !Silently ignore unsupported field options
    end select
  end subroutine fields_init_response

  subroutine init_fields
!CMR,18/2/2011:
! add optional logical arg "noalloc" to avoid array allocations in ingen
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn
    use init_g, only: ginit, init_init_g
    use nonlinear_terms, only: nl_finish_init => finish_init
    use antenna, only: init_antenna
    use kt_grids, only: gridopt_switch, gridopt_box, kwork_filter
    implicit none
    logical :: restarted
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    call fields_pre_init

    call fields_init_response

! Turn on nonlinear terms.
    if (debug) write(6,*) "init_fields: nl_finish_init"
    call nl_finish_init

    ! EGH Commented out the following lines as they are now
    ! handled by gs2_init

    !if (debug) write(6,*) "init_fields: ginit"
    !call ginit (restarted)
    !if (restarted .and. .not. force_maxwell_reinit) return
    !if (debug) write(6,*) "init_fields: init_antenna"
    !call init_antenna

    !!Set the initial fields
    !call set_init_fields

    !If running in flux tube disable evolution of ky=kx=0 mode
    if(gridopt_switch.eq.gridopt_box) kwork_filter(1,1)=.true.
  end subroutine init_fields

  !>Force the current 
  subroutine dump_response_to_file(suffix)
    use fields_implicit, only: dump_response_to_file_imp
    use fields_local, only: dump_response_to_file_local
    implicit none
    character(len=*), intent(in), optional :: suffix 
    !Note can pass optional straight through as long as also optional
    !in called routine (and not different routines combined in interface)
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call dump_response_to_file_imp(suffix)
    case (fieldopt_local)
       call dump_response_to_file_local(suffix)
    case default
       !Silently ignore unsupported field options
    end select
  end subroutine dump_response_to_file

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
    use fields_implicit, only: field_subgath
    use fields_local, only: minNrow
    use fields_local, only: do_smart_update, field_local_allreduce, field_local_allreduce_sub
    use fields_arrays, only: response_file
    use file_utils, only: run_name
    implicit none
    type (text_option), dimension (5), parameter :: fieldopts = &
         (/ text_option('default', fieldopt_implicit), &
            text_option('implicit', fieldopt_implicit), &
            text_option('test', fieldopt_test),&
            text_option('local', fieldopt_local),&
            text_option('implicit_local', fieldopt_local)/)
    character(20) :: field_option
    character(len=256) :: response_dir
    namelist /fields_knobs/ field_option, remove_zonal_flows_switch, field_subgath, force_maxwell_reinit,&
         dump_response, read_response, minNrow, do_smart_update, field_local_allreduce, field_local_allreduce_sub, response_dir
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
       field_local_allreduce=.false.
       field_local_allreduce_sub=.false.
       response_dir=''
       in_file = input_unit_exist ("fields_knobs", exist)
       if (exist) read (unit=in_file, nml=fields_knobs)

       ierr = error_unit()
       call get_option_value &
            (field_option, fieldopts, fieldopt_switch, &
            ierr, "field_option in fields_knobs",.true.)

       if(trim(response_dir).eq.'')then
          write(response_file,'(A)') trim(run_name)
       else
          write(response_file,'(A,"/",A)') trim(response_dir),trim(run_name)
       endif

    end if

    call broadcast (fieldopt_switch)
    call broadcast (remove_zonal_flows_switch)
    call broadcast (field_subgath)
    call broadcast (force_maxwell_reinit)
    call broadcast (dump_response)
    call broadcast (read_response)

    !Setup response file location
    call broadcast(response_dir)
    call broadcast(response_file)

    !Set the solve type specific flags
    call set_dump_and_read_response(dump_response, read_response)
    select case (fieldopt_switch)
    case (fieldopt_implicit)
    case (fieldopt_test)
    case (fieldopt_local)
       call broadcast (minNrow)
       call broadcast (do_smart_update)
       call broadcast (field_local_allreduce)
       call broadcast (field_local_allreduce_sub)
    end select
  end subroutine read_parameters

  subroutine set_dump_and_read_response(dump_flag, read_flag)
    use fields_implicit, only: dump_response_imp => dump_response, read_response_imp=>read_response
    use fields_local, only: dump_response_loc => dump_response, read_response_loc=>read_response
    implicit none
    logical, intent(in) :: dump_flag, read_flag
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       dump_response_imp=dump_flag
       read_response_imp=read_flag
    case (fieldopt_local)
       dump_response_loc=dump_flag
       read_response_loc=read_flag
    case default
       !Silently ignore unsupported field types
    end select
  end subroutine set_dump_and_read_response

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use antenna, only: no_driver
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew, apar_ext
    implicit none

    if (.not. allocated(phi)) then
       allocate (     phi (-ntgrid:ntgrid,ntheta0,naky))
       allocate (    apar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (   bpar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (  phinew (-ntgrid:ntgrid,ntheta0,naky))
       allocate ( aparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparnew (-ntgrid:ntgrid,ntheta0,naky))
    endif
    phi = 0.; phinew = 0.
    apar = 0.; aparnew = 0.
    bpar = 0.; bparnew = 0.
    if(.not.allocated(apar_ext).and.(.not.no_driver))then
       allocate (apar_ext (-ntgrid:ntgrid,ntheta0,naky))
       apar_ext = 0.
    endif
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

  !This routine has a potentially misleading name and isn't used anywhere
  subroutine kperp (ntgrid_output, akperp)
    use theta_grid, only: delthet
    use kt_grids, only: naky, aky, ntheta0, kperp2
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phinew, aparnew, bparnew
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

  !This routine isn't used anywhere
  subroutine fieldlineavgphi_loc (ntgrid_output, it, ik, phiavg)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    use fields_arrays, only: phi
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

  !This doesn't look like a useful routine and isn't used anywhere
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
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    implicit none
    initialized  = .false.
    phi = 0.
    phinew = 0.
    apar = 0.
    aparnew = 0.
    bpar = .0
    bparnew = 0.
    !What about apar_ext?
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call fi_reset
    case (fieldopt_test)
       call ft_reset
    case (fieldopt_local)
       call fl_reset
    end select
  end subroutine reset_init

  subroutine finish_fields

    use fields_implicit, only: implicit_reset => reset_init
    use fields_test, only: test_reset => reset_init
    use fields_local, only: finish_fields_local
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: apar_ext

    implicit none

    initialized  = .false.
    phi = 0.
    phinew = 0.
    apar = 0.
    aparnew = 0.
    bpar = .0
    bparnew = 0.
    
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call implicit_reset
    case (fieldopt_test)
       call test_reset
    case (fieldopt_local)
       call finish_fields_local
    end select

    if (allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew)
    if (allocated(apar_ext)) deallocate (apar_ext)

  end subroutine finish_fields

end module fields
