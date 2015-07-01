!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters. In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.
 

module run_parameters
  implicit none

  public :: init_run_parameters, finish_run_parameters
  public :: check_run_parameters, wnml_run_parameters
  public :: write_trinity_parameters


  public :: beta, zeff, tite, reset, immediate_reset
  public :: fphi, fapar, fbpar
!  public :: delt, delt_max, wunits, woutunits, tunits
  public :: code_delt_max, wunits, woutunits, tunits
  public :: nstep, wstar_units, eqzip, margin
  public :: secondary, tertiary, harris
  public :: ieqzip
  public :: k0
  public :: vnm_init
  public :: avail_cpu_time, margin_cpu_time
!  public :: include_lowflow, rhostar, neo_test
  public :: rhostar, neo_test
  public :: do_eigsolve
  !> If true and nonlinear_mode is "off", return
  !! simple diffusive estimates of fluxes to trinity
  public :: trinity_linear_fluxes
  public :: user_comments

  private

  real :: beta, zeff, tite
  real :: fphi, fapar, fbpar, faperp
  real :: delt, code_delt_max, user_delt_max, margin
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  real, dimension (2) :: vnm_init
  real :: avail_cpu_time, margin_cpu_time
  integer :: nstep
  logical :: reset=.false.
  logical :: immediate_reset
  logical :: wstar_units, eqzip
  logical :: secondary, tertiary, harris
  real :: k0
  integer, public :: delt_option_switch
  integer, public, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: initialized = .false.
  logical :: rpexist, knexist
  real :: rhostar
!  logical :: include_lowflow, neo_test
  logical :: neo_test

  logical :: trinity_linear_fluxes, do_eigsolve
  character(len=100000) :: user_comments

  integer, allocatable :: ieqzip(:,:)
  integer :: eqzip_option_switch
  integer, parameter :: &
       eqzip_option_none = 1, &
       eqzip_option_secondary = 2, &
       eqzip_option_tertiary = 3, &
       eqzip_option_equilibrium = 4

contains

  subroutine check_run_parameters(report_unit)
  implicit none
  integer :: report_unit
    if (fphi /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fphi in the knobs namelist = ',e11.4)") fphi
       write (report_unit, fmt="('fphi is a scale factor of all instances of Phi (the electrostatic potential).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (fapar == 0.) then
       write (report_unit, fmt="('A_parallel will not be included in the calculation.')")
    end if
    if (fapar == 1.) then
       write (report_unit, fmt="('A_parallel will be included in the calculation.')")
    end if
    if (fapar /= 0. .and. fapar /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fapar in the knobs namelist = ',e11.4)") fapar
       write (report_unit, fmt="('fapar is a scale factor of all instances of A_parallel (the parallel vector potential).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (fbpar == 0.) then
       write (report_unit, fmt="('B_parallel will not be included in the calculation.')")
    end if
    if (fbpar == 1.) then
       write (report_unit, fmt="('B_parallel will be included in the calculation.')")
    end if
    if (fbpar /= 0. .and. fbpar /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fbpar in the knobs namelist = ',e11.4)") fbpar
       write (report_unit, fmt="('fbpar is a scale factor of all instances of B_parallel &
           & (the perturbed parallel magnetic field).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (eqzip) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('eqzip = T in the knobs namelist.')")
       write (report_unit, fmt="('This freezes some modes in time for a secondary stability analysis.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
       if (secondary) write (report_unit, fmt="('Mode with kx = 0, ky = ky_min fixed in time')")
       if (tertiary)  write (report_unit, fmt="('Mode with ky = 0, kx = kx_min fixed in time')")
    end if

    write (report_unit, *) 
    if(immediate_reset)then
       write (report_unit, fmt="('The time step will be reset immediately after cfl violation detected.')") 
    else
       write (report_unit, fmt="('The time step will be reset just before the next time step after cfl violation detected.')") 
    endif
    write (report_unit, *) 
  end subroutine check_run_parameters

  subroutine wnml_run_parameters(unit,electrons,collisions)
  implicit none
  integer :: unit
  logical :: electrons, collisions
     if (rpexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "parameters"
       write (unit, fmt="(' beta = ',e17.10)") beta       ! if zero, fapar, fbpar should be zero
       if (collisions) write (unit, fmt="(' zeff = ',e17.10)") zeff
       if (.not. electrons)  write (unit, fmt="(' tite = ',e17.10)") tite
!CMR, 10/2/2011: zip not in this namelist, so removing it!
!       if (zip) write (unit, fmt="(' zip = ',L1)") zip
       write (unit, fmt="(' /')")
     endif
     if (knexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "knobs"
       write (unit, fmt="(' fphi   = ',f6.3)") fphi
       write (unit, fmt="(' fapar  = ',f6.3)") fapar
       write (unit, fmt="(' fbpar = ',f6.3)") fbpar
       write (unit, fmt="(' delt = ',e17.10)") delt
       write (unit, fmt="(' nstep = ',i8)") nstep
       write (unit, fmt="(' wstar_units = ',L1)") wstar_units
       if (eqzip) then
          write (unit, fmt="(' eqzip = ',L1)") eqzip
          write (unit, fmt="(' secondary = ',L1)") secondary
          write (unit, fmt="(' tertiary = ',L1)") tertiary
       end if
       write (unit, fmt="(' margin = ',e17.10)") margin
       select case (delt_option_switch)
       case (delt_option_auto)
          write (unit, fmt="(' delt_option = ',a)") '"check_restart"'
       case (delt_option_hand)
          ! nothing
       end select
       write (unit, fmt="(' immediate_reset = ',L1)") immediate_reset
       write (unit, fmt="(' /')")
     endif
  end subroutine wnml_run_parameters

  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky, nakx => ntheta0
    use gs2_time, only: init_delt, user2code
    
    implicit none
!    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_kt_grids
    call init_delt (delt)
    call user2code (user_delt_max, code_delt_max)

    if(.not. allocated(wunits)) allocate (wunits(naky))
    if(.not. allocated(woutunits)) allocate (woutunits(naky))
    if(.not. allocated(tunits)) allocate (tunits(naky))

! omega_* normalization of time: 
    call adjust_time_norm

    if(.not.allocated(ieqzip)) allocate(ieqzip(nakx,naky))
    ieqzip(1:nakx,1:naky)=1
    select case (eqzip_option_switch)
    case (eqzip_option_secondary)
       ! suppress evolution of secondary mode
       ieqzip(1,2) = 0
    case (eqzip_option_tertiary)
       ! suppress evolution of tertiary mode
       ieqzip(2,1) = 0
       ieqzip(nakx,1) = 0
    case (eqzip_option_equilibrium)
       ! suppress evolution of 1D equilibrium (x dependent)
       ieqzip(1:nakx,1) = 0
    end select
  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use gs2_save, only: init_dt, init_vnm
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (4), parameter :: eqzipopts = &
         (/ text_option('none', eqzip_option_none), &
            text_option('secondary', eqzip_option_secondary), &
            text_option('tertiary', eqzip_option_tertiary), &
            text_option('equilibrium', eqzip_option_equilibrium) /)
    character (len=20) :: eqzip_option
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option
    integer :: ierr, istatus, in_file
    real :: delt_saved
    real, dimension (2) :: vnm_saved

    real :: teti  ! for back-compatibility
    namelist /parameters/ beta, zeff, tite, teti, k0, rhostar, user_comments
    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, wstar_units, eqzip, &
         delt_option, margin, secondary, tertiary, faperp, harris, &
!         avail_cpu_time, eqzip_option, include_lowflow, neo_test
         avail_cpu_time, margin_cpu_time, eqzip_option, neo_test, &
         trinity_linear_fluxes, do_eigsolve, immediate_reset

    if (proc0) then
       fbpar = -1.0
       faperp = 0.0
       beta = 0.0
       zeff = 1.0
       tite = 1.0
       teti = -100.0
       rhostar = 3.e-3
!       include_lowflow = .false.
       neo_test = .false.
       wstar_units = .false.
       eqzip_option = 'none'
       eqzip = .false.
       secondary = .true.
       tertiary = .false.
       harris = .false.
       k0 = 1.
       delt_option = 'default'
       margin = 0.05
       avail_cpu_time = 1.e10
       margin_cpu_time = 300.
       trinity_linear_fluxes = .false.
       do_eigsolve = .false.
       immediate_reset = .true.
       user_comments = ''
       
       in_file = input_unit_exist("parameters", rpexist)
!       if (rpexist) read (unit=input_unit("parameters"), nml=parameters)
       if (rpexist) read (unit=in_file,nml=parameters)

       in_file = input_unit_exist("knobs", knexist)
!       if (knexist) read (unit=input_unit("knobs"), nml=knobs)
       if (knexist) read (unit=in_file, nml=knobs)

       if (teti /= -100.0) tite = teti

! Allow faperp-style initialization for backwards compatibility.
! Only fbpar is used outside of this subroutine.
       if (fbpar == -1.) then
          fbpar = faperp
       end if

       if (eqzip) then
          if (secondary .and. tertiary) then
             ierr = error_unit()
             write (ierr, *) 'Forcing secondary = FALSE'
             write (ierr, *) 'because you have chosen tertiary = TRUE'
             secondary = .false.
          end if
          if (secondary .and. harris) then
             ierr = error_unit()
             write (ierr, *) 'Forcing secondary = FALSE'
             write (ierr, *) 'because you have chosen harris = TRUE'
             secondary = .false.
          end if
          if (tertiary .and. harris) then
             ierr = error_unit()
             write (ierr, *) 'Forcing tertiary = FALSE'
             write (ierr, *) 'because you have chosen harris = TRUE'
             tertiary = .false.
          end if
       endif

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs",.true.)

       call get_option_value ( &
            eqzip_option, eqzipopts, eqzip_option_switch, error_unit(), &
            "eqzip_option in knobs",.true.)

!!$       ! eqzip_option replaces eqzip, secondary, tertiary, harris
!!$       if (eqzip .and. eqzip_option_switch == eqzip_option_none) then
!!$          if (harris) then
!!$             eqzip_option_switch = eqzip_option_equilibrium
!!$             write(error_unit(),*) 'eqzip_option is set to equilibrium'
!!$          else
!!$             if (tertiary) then
!!$                eqzip_option_switch = eqzip_option_tertiary
!!$                write(error_unit(),*) 'eqzip_option is set to tertiary'
!!$             else
!!$                eqzip_option_switch = eqzip_option_secondary
!!$                write(error_unit(),*) 'eqzip_option is set to secondary'
!!$             end if
!!$          end if
!!$       else if (eqzip .and. eqzip_option_switch /= eqzip_option_none) then
!!$          write(error_unit(),*) 'eqzip, secondary, tertiary, harris are ignored'
!!$          write(error_unit(),*) 'because eqzip_option exists'
!!$       end if

    end if

    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (beta)
    call broadcast (zeff)
    call broadcast (tite)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (fbpar)
    call broadcast (nstep)
    call broadcast (wstar_units)
    call broadcast (eqzip)
    call broadcast (secondary)
    call broadcast (tertiary)
    call broadcast (harris)
    call broadcast (margin)
    call broadcast (k0)
    call broadcast (avail_cpu_time)
    call broadcast (margin_cpu_time)
    call broadcast (eqzip_option_switch)
!    call broadcast (include_lowflow)
    call broadcast (rhostar)
    call broadcast (neo_test)
    call broadcast (trinity_linear_fluxes)
    call broadcast (do_eigsolve)
    call broadcast (immediate_reset)

    user_delt_max = delt

    delt_saved = delt
    if (delt_option_switch == delt_option_auto) then
       vnm_init = 1.0
       call init_vnm (vnm_saved, istatus)
       if (istatus == 0) vnm_init = vnm_saved
       call init_dt (delt_saved, istatus)
       if (istatus == 0) delt  = delt_saved
    endif

  end subroutine read_parameters
  
  subroutine write_trinity_parameters(trinpars_unit)
    integer, intent(in) :: trinpars_unit
    write(trinpars_unit, "(A15)") '&run_parameters'
    write (trinpars_unit, *) ' beta = ', beta
    write (trinpars_unit, "(A1)") '/'

  end subroutine write_trinity_parameters

  subroutine adjust_time_norm
    use file_utils, only: error_unit
    use mp, only: proc0
    use kt_grids, only: aky
    implicit none
!CMR: Sep 2010
! Attempt to understand time normalisation variables, which are arrays(naky)
!    TUNITS: DT(KY)=TUNITS(KY).CODE_DT
!            This is a generally very useful variable to store ky dependent 
!            timestep in the code time normalisation.
!            Used to multiply ky independent source terms on RHS of GKE.
!    WUNITS: WUNITS(KY)=AKY(KY)*TUNITS(KY)/2
!            Auxiliary variable.  Used to save compute operations when 
!            evaluating source terms on RHS of GKE that are proportional to ky.
!            !! The Mysterious factor 1/2 Explained !!
!            The factor of 1/2 arises because those source terms were first
!            specified using the normalisation Tref=mref vtref^2 
! [R Numata et al, "AstroGK: Astrophysical gyrokinetics code", JCP, 2010].
!CMRend
    if (wstar_units) then
       wunits = 1.0
       where (aky /= 0.0)
          tunits = 2.0/aky
       elsewhere
          tunits = 0.0
       end where
       if (any(tunits == 0.0) .and. proc0) then
          write (error_unit(), *) &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
          print *, &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
       end if
!CMR: Sep 2010
!  Changes to allow wstar_units to be used
!CMRend
    else
       tunits = 1.0
       wunits = aky/2.0
    end if
    woutunits = 1.0/tunits
  end subroutine adjust_time_norm

  subroutine finish_run_parameters

    implicit none

    if (allocated(wunits)) deallocate (wunits, woutunits, tunits)
    if (allocated(ieqzip)) deallocate (ieqzip)

    initialized = .false.

  end subroutine finish_run_parameters

end module run_parameters
