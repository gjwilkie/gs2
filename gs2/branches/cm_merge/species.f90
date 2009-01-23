module species
  implicit none

  public :: init_species
  public :: nspec, specie, spec
  public :: ion_species, electron_species, slowing_down_species, tracer_species
  public :: has_electron_species, has_slowing_down_species

  type :: specie
     real :: z
     real :: mass
     real :: dens, dens0, u0
     real :: temp
     real :: tprim
     real :: fprim
     real :: uprim, uprim2
     real :: vnewk, nustar
     real :: nu, nu_h  ! nu will be the preferred collisionality parameter moving forward
                       ! as it will have a consistent v_t scaling
                       ! nu_h controls a hyperviscous term embedded in the collision operator
     real :: stm, zstm, tz, smz, zt
     integer :: type
  end type specie

  private

  integer, parameter :: ion_species = 1
  integer, parameter :: electron_species = 2 ! for collision operator
  integer, parameter :: slowing_down_species = 3 ! slowing-down distn
  integer, parameter :: tracer_species = 4 ! for test particle diffusion studies

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
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    real :: z, mass, dens, dens0, u0, temp, tprim, fprim, uprim, uprim2, vnewk, nustar, nu, nu_h
    character(20) :: type
    integer :: unit
    integer :: is
    namelist /species_knobs/ nspec
    namelist /species_parameters/ z, mass, dens, dens0, u0, temp, &
         tprim, fprim, uprim, uprim2, vnewk, nustar, type, nu, nu_h
    integer :: ierr, in_file
    logical :: exist

    type (text_option), dimension (9), parameter :: typeopts = &
         (/ text_option('default', ion_species), &
            text_option('ion', ion_species), &
            text_option('electron', electron_species), &
            text_option('e', electron_species), &
            text_option('beam', slowing_down_species), &
            text_option('fast', slowing_down_species), &
            text_option('alpha', slowing_down_species), &
            text_option('slowing-down', slowing_down_species), &
            text_option('trace', tracer_species) /)

    if (proc0) then
       nspec = 2
       in_file = input_unit_exist("species_knobs", exist)
!       if (exist) read (unit=input_unit("species_knobs"), nml=species_knobs)
       if (exist) read (unit=in_file, nml=species_knobs)
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
          z = 1
          mass = 1.0
          dens = 1.0
          dens0 = 1.0
          u0 = 1.0
          temp = 1.0
          tprim = 6.9
          fprim = 2.2
          uprim = 0.0
          uprim2 = 0.0
          nustar = -1.0
          vnewk = 0.0
          nu = -1.0
          nu_h = 0.0
          type = "default"
          read (unit=unit, nml=species_parameters)
          close (unit=unit)

          spec(is)%z = z
          spec(is)%mass = mass
          spec(is)%dens = dens
          spec(is)%dens0 = dens0
          spec(is)%u0 = u0
          spec(is)%temp = temp
          spec(is)%tprim = tprim
          spec(is)%fprim = fprim
          spec(is)%uprim = uprim
          spec(is)%uprim2 = uprim2
          spec(is)%vnewk = vnewk
          spec(is)%nustar = nustar
          spec(is)%nu = nu
          spec(is)%nu_h = nu_h

          spec(is)%stm = sqrt(temp/mass)
          spec(is)%zstm = z/sqrt(temp*mass)
          spec(is)%tz = temp/z
          spec(is)%zt = z/temp
          spec(is)%smz = abs(sqrt(temp*mass)/z)

          ierr = error_unit()
          call get_option_value (type, typeopts, spec(is)%type, ierr, "type in species_parameters_x")
       end do
    end if

    do is = 1, nspec
       call broadcast (spec(is)%z)
       call broadcast (spec(is)%mass)
       call broadcast (spec(is)%dens)
       call broadcast (spec(is)%dens0)
       call broadcast (spec(is)%u0)
       call broadcast (spec(is)%temp)
       call broadcast (spec(is)%tprim)
       call broadcast (spec(is)%fprim)
       call broadcast (spec(is)%uprim)
       call broadcast (spec(is)%uprim2)
       call broadcast (spec(is)%vnewk)
       call broadcast (spec(is)%nu)
       call broadcast (spec(is)%nu_h)
       call broadcast (spec(is)%nustar)
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
