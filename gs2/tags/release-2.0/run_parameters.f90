module run_parameters
  implicit none

  public :: init_run_parameters

  public :: beta, zeff, teti
  public :: fphi, fapar, faperp
  public :: delt, wunits, woutunits, tunits
  public :: nstep

  private

  real :: beta, zeff, teti
  real :: fphi, fapar, faperp
  real :: delt
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  integer :: nstep

contains

  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky
    use mp, only: proc0
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_kt_grids

    allocate (wunits(naky))
    allocate (woutunits(naky))
    allocate (tunits(naky))
    if (proc0) call read_parameters
    call broadcast_parameters
  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit
    use kt_grids, only: aky
    implicit none
    namelist /parameters/ beta, zeff, teti
    logical :: wstar_units
    namelist /knobs/ fphi, fapar, faperp, delt, nstep, wstar_units

    beta = 0.0
    zeff = 1.0
    teti = 1.0
    wstar_units = .false.
    read (unit=input_unit("parameters"), nml=parameters)
    read (unit=input_unit("knobs"), nml=knobs)

    if (wstar_units) then
       wunits = 1.0
       woutunits = aky/sqrt(2.0)
       where (aky /= 0.0)
          tunits = 2.0/aky
       elsewhere
          tunits = 0.0
       end where
    else
       wunits = aky/2.0
       woutunits = sqrt(2.0)
       tunits = 1.0
    end if
  end subroutine read_parameters

  subroutine broadcast_parameters
    use mp, only: broadcast
    implicit none
    call broadcast (beta)
    call broadcast (zeff)
    call broadcast (teti)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (faperp)
    call broadcast (delt)
    call broadcast (nstep)
    call broadcast (wunits)
    call broadcast (woutunits)
    call broadcast (tunits)
  end subroutine broadcast_parameters

end module run_parameters
