module run_parameters
  implicit none

  public :: init_run_parameters

  public :: beta, zeff, tite, rhostar
  public :: fphi, fapar, fbpar
!  public :: delt, delt_max, wunits, woutunits, tunits, funits, tnorm
  public :: code_delt_max, wunits, woutunits, tunits, funits, tnorm
  public :: nstep, wstar_units, eqzip, margin
  public :: secondary, tertiary, harris
  public :: k0
  public :: vnm_init

  private

  real :: beta, zeff, tite, rhostar
  real :: fphi, fapar, fbpar, faperp
  real :: delt, code_delt_max, user_delt_max, funits, tnorm, margin
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  real, dimension (2) :: vnm_init
  integer :: nstep
  logical :: wstar_units, eqzip
  logical :: secondary, tertiary, harris
  real :: k0
  integer :: delt_option_switch
  integer, parameter :: delt_option_hand = 1, delt_option_auto = 2

contains

  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky
    use gs2_time, only: init_delt, user2code
    
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_kt_grids (tnorm)
    call init_delt (delt, tnorm)
    call user2code (user_delt_max, code_delt_max)

!    delt = delt * tnorm

    allocate (wunits(naky))
    allocate (woutunits(naky))
    allocate (tunits(naky))

! omega_* normalization of time: 
    call adjust_time_norm

  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use gs2_save, only: init_dt, init_vnm
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option
    integer :: ierr, istatus, in_file
    real :: delt_saved
    real, dimension (2) :: vnm_saved

    real :: teti  ! for back-compatibility
    logical :: exist
    namelist /parameters/ beta, zeff, tite, rhostar, teti, k0
    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, wstar_units, eqzip, &
         delt_option, margin, secondary, tertiary, faperp, harris

    if (proc0) then
       fbpar = -1.0
       faperp = 0.0
       beta = 0.0
       zeff = 1.0
       tite = 1.0
       teti = -100.0
       rhostar = 0.1
       wstar_units = .false.
       eqzip = .false.
       secondary = .true.
       tertiary = .false.
       harris = .false.
       k0 = 1.
       delt_option = 'default'
       margin = 0.05
       in_file = input_unit_exist("parameters", exist)
!       if (exist) read (unit=input_unit("parameters"), nml=parameters)
       if (exist) read (unit=in_file,nml=parameters)

       in_file = input_unit_exist("knobs", exist)
!       if (exist) read (unit=input_unit("knobs"), nml=knobs)
       if (exist) read (unit=in_file, nml=knobs)

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
       end if

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs")
    end if

    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (beta)
    call broadcast (zeff)
    call broadcast (tite)
    call broadcast (rhostar)
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

  subroutine adjust_time_norm
    use file_utils, only: error_unit
    use mp, only: proc0
    use kt_grids, only: aky
    implicit none

    if (wstar_units) then
       funits = 1.0
       wunits = 1.0
       woutunits = aky/sqrt(2.0)
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
    else
       tunits = 1.0
       wunits = aky/2.0
       funits = tnorm
       woutunits = tnorm
    end if

  end subroutine adjust_time_norm

end module run_parameters
