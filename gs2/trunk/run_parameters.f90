module run_parameters
  implicit none

  public :: init_run_parameters

  public :: beta, zeff, tite, rhostar
  public :: fphi, fapar, faperp
  public :: delt, wunits, woutunits, tunits
  public :: nstep, wstar_units, quick

  private

  real :: beta, zeff, tite, rhostar
  real :: fphi, fapar, faperp
  real :: delt
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  integer :: nstep
  logical :: wstar_units, quick

contains

  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky
    use mp, only: proc0
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters
    call init_kt_grids

    allocate (wunits(naky))
    allocate (woutunits(naky))
    allocate (tunits(naky))

    call adjust_time_norm

  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit
    use mp, only: proc0, broadcast
    implicit none

    real :: teti  ! for back-compatibility
    namelist /parameters/ beta, zeff, tite, rhostar, teti
    namelist /knobs/ fphi, fapar, faperp, delt, nstep, wstar_units, quick

    if (proc0) then
       beta = 0.0
       zeff = 1.0
       tite = 1.0
       teti = -100.0
       rhostar = 0.1
       wstar_units = .false.
       quick = .false.
       read (unit=input_unit("parameters"), nml=parameters)
       read (unit=input_unit("knobs"), nml=knobs)
       if (teti /= -100.0) tite = teti
    endif
    call broadcast (beta)
    call broadcast (zeff)
    call broadcast (tite)
    call broadcast (rhostar)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (faperp)
    call broadcast (delt)
    call broadcast (nstep)
    call broadcast (wstar_units)
    call broadcast (quick)

  end subroutine read_parameters

  subroutine adjust_time_norm
    use file_utils, only: error_unit
    use mp, only: proc0
    use kt_grids, only: aky
    implicit none

    if (wstar_units) then
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
       wunits = aky/2.0
       woutunits = sqrt(2.0)
       tunits = 1.0
    end if

  end subroutine adjust_time_norm

end module run_parameters
