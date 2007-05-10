module run_parameters
  implicit none

  public :: init_run_parameters

  public :: beta, zeff, tite, rhostar
  public :: fphi, fapar, faperp
  public :: delt, delt_max, wunits, woutunits, tunits, funits, tnorm
  public :: nstep, wstar_units, eqzip, margin

  private

  real :: beta, zeff, tite, rhostar
  real :: fphi, fapar, faperp
  real :: delt, delt_max, funits, tnorm, margin
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  integer :: nstep
  logical :: wstar_units, eqzip

  integer :: delt_option_switch
  integer, parameter :: delt_option_hand = 1, delt_option_auto = 2

contains

  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters
    call init_kt_grids (tnorm)

! possible internal sqrt(2.) normalization of time:
    delt = delt * tnorm

    allocate (wunits(naky))
    allocate (woutunits(naky))
    allocate (tunits(naky))

! omega_* normalization of time: 
    call adjust_time_norm

  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit
    use mp, only: proc0, broadcast
    use gs2_save, only: init_dt
    use text_options
    implicit none
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option
    integer :: ierr, istatus
    real :: delt_saved

    real :: teti  ! for back-compatibility
    namelist /parameters/ beta, zeff, tite, rhostar, teti
    namelist /knobs/ fphi, fapar, faperp, delt, nstep, wstar_units, eqzip, &
         delt_option, margin

    if (proc0) then
       beta = 0.0
       zeff = 1.0
       tite = 1.0
       teti = -100.0
       rhostar = 0.1
       wstar_units = .false.
       eqzip = .false.
       delt_option = 'default'
       margin = 0.05
       read (unit=input_unit("parameters"), nml=parameters)
       read (unit=input_unit("knobs"), nml=knobs)
       if (teti /= -100.0) tite = teti

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
    call broadcast (faperp)
    call broadcast (nstep)
    call broadcast (wstar_units)
    call broadcast (eqzip)
    call broadcast (margin)
    delt_max = delt

    delt_saved = delt
    if (delt_option_switch == delt_option_auto) then
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
