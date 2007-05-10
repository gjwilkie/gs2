module species
  implicit none

  public :: init_species
  public :: nspec, specie, spec
  public :: ion_species, electron_species, slowing_down_species
  public :: has_electron_species, has_slowing_down_species

  type :: specie
     real :: z
     real :: mass
     real :: dens
     real :: temp
     real :: tprim
     real :: fprim
     real :: uprim, uprim2
     real :: vnewk, vnewk4
     real :: stm, zstm, tz, smz, zt
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
    use file_utils, only: input_unit, error_unit, get_indexed_namelist_unit, input_unit_exist
    use text_options
    use mp, only: proc0, broadcast
    implicit none
    real :: z, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, vnewk4
    character(20) :: type
    integer :: unit
    integer :: is
    namelist /species_knobs/ nspec
    namelist /species_parameters/ &
         z, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, vnewk4, type
    integer :: ierr, in_file
    logical :: exist

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
       in_file = input_unit_exist("species_knobs", exist)
       if (exist) read (unit=input_unit("species_knobs"), nml=species_knobs)
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
          uprim = 0.0
          uprim2 = 0.0
          vnewk = 0.0
          vnewk4 = 0.0
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
          spec(is)%vnewk4 = vnewk4

          spec(is)%stm = sqrt(temp/mass)
          spec(is)%zstm = z/sqrt(temp*mass)
          spec(is)%tz = temp/z
          spec(is)%zt = z/temp
          spec(is)%smz = abs(sqrt(temp*mass)/z)

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
       call broadcast (spec(is)%vnewk4)
       call broadcast (spec(is)%stm)
       call broadcast (spec(is)%zstm)
       call broadcast (spec(is)%tz)
       call broadcast (spec(is)%zt)
       call broadcast (spec(is)%smz)
       call broadcast (spec(is)%type)
    end do
  end subroutine read_parameters

  pure function has_electron_species (spec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    logical :: has_electron_species
    has_electron_species = any(spec%type == electron_species)
  end function has_electron_species

  pure function has_slowing_down_species (spec)
    implicit none
    type (specie), dimension (:), intent (in) :: spec
    logical :: has_slowing_down_species
    has_slowing_down_species = any(spec%type == slowing_down_species)
  end function has_slowing_down_species

end module species
