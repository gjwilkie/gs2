module species
  implicit none

  public :: init_species
  public :: nspec, specie, spec
  public :: ion_species, electron_species, slowing_down_species
  public :: stm, zstm, tz, smz

  type :: specie
     real :: z
     real :: mass
     real :: dens
     real :: temp
     real :: tprim
     real :: fprim
     real :: uprim, uprim2
     real :: vnewk
     real :: e_flat
     real :: temp_flat
     integer :: type
  end type specie

  private

  integer, parameter :: ion_species = 1
  integer, parameter :: electron_species = 2 ! for collision operator
  integer, parameter :: slowing_down_species = 3 ! slowing-down distn

  integer :: nspec
  type (specie), dimension (:), allocatable :: spec

contains

  subroutine init_species
    implicit none
    logical, save :: initialized = .false.
    if (initialized) return
    initialized = .true.

    call read_parameters
  end subroutine init_species

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, get_indexed_namelist_unit
    use text_options
    use mp, only: proc0, broadcast
    implicit none
    real :: z, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, &
         e_flat, temp_flat
    character(20) :: type
    integer :: unit
    integer :: is
    namelist /species_knobs/ nspec
    namelist /species_parameters/ &
         z, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, type, &
         e_flat, temp_flat
    integer :: ierr

    type (text_option), dimension (8), parameter :: typeopts = &
         (/ text_option('default', ion_species), &
            text_option('ion', ion_species), &
            text_option('electron', electron_species), &
            text_option('e', electron_species), &
            text_option('beam', slowing_down_species), &
            text_option('fast', slowing_down_species), &
            text_option('alpha', slowing_down_species), &
            text_option('slowing-down', slowing_down_species) /)

    if (proc0) then
       nspec = 2
       read (unit=input_unit("species_knobs"), nml=species_knobs)
       if (nspec < 1) then
          ierr = error_unit()
          write (unit=ierr, &
               fmt="('Invalid nspec in species_knobs: ', i5)") nspec
          stop
       end if
    end if

    call broadcast (nspec)
    allocate (spec(nspec))

    if (proc0) then
       do is = 1, nspec
          call get_indexed_namelist_unit (unit, "species_parameters", is)
          e_flat = -1.
          temp_flat = 1.
          uprim = 0.0
          uprim2 = 0.0
          vnewk = 0.0
          type = "default"
          read (unit=unit, nml=species_parameters)
          close (unit=unit)

          spec(is)%z = z
          spec(is)%mass = mass
          spec(is)%dens = dens
          spec(is)%temp = temp
          spec(is)%tprim = tprim
          spec(is)%fprim = fprim
          spec(is)%uprim = uprim
          spec(is)%uprim2 = uprim2
          spec(is)%vnewk = vnewk
          spec(is)%e_flat = e_flat
          spec(is)%temp_flat = temp_flat
          ierr = error_unit()
          call get_option_value &
               (type, typeopts, spec(is)%type, &
               ierr, "type in species_parameters_x")
       end do
    end if

    do is = 1, nspec
       call broadcast (spec(is)%z)
       call broadcast (spec(is)%mass)
       call broadcast (spec(is)%dens)
       call broadcast (spec(is)%temp)
       call broadcast (spec(is)%tprim)
       call broadcast (spec(is)%fprim)
       call broadcast (spec(is)%uprim)
       call broadcast (spec(is)%uprim2)
       call broadcast (spec(is)%vnewk)
       call broadcast (spec(is)%type)
       call broadcast (spec(is)%e_flat)
       call broadcast (spec(is)%temp_flat)
    end do
  end subroutine read_parameters

  pure function stm (spec, ispec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    integer, intent (in) :: ispec
    real :: stm
    stm = sqrt(spec(ispec)%temp/spec(ispec)%mass)
  end function stm

  pure function zstm (spec, ispec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    integer, intent (in) :: ispec
    real :: zstm
    zstm = spec(ispec)%z/sqrt(spec(ispec)%temp*spec(ispec)%mass)
  end function zstm

  pure function tz (spec, ispec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    integer, intent (in) :: ispec
    real :: tz
    tz = spec(ispec)%temp/spec(ispec)%z
  end function tz

  pure function smz (spec, ispec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    integer, intent (in) :: ispec
    real :: smz
    smz = abs(sqrt(spec(ispec)%temp*spec(ispec)%mass)/spec(ispec)%z)
  end function smz

end module species
