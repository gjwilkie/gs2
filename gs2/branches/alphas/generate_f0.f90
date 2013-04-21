
!> A module which allows the simulation of species with a
!! arbitrary background distribution function (e.g. alphas, 
!! beam ions).
!! 
!! Depending on the species type and the input parameters,
!! it will return either a Maxwellian or some other analytic
!! function or a solution for F_0 returned from an external
!! solver. 
!!
!! This is free software released under the GPLv3
!! Written by
!!    Edmund Highcock (edmundhighcock@sourceforge.net)
!!    George Wilkie  (gwilkie@umd.edu)
!!  

module general_f0


  !> Takes in the gridpoints in energy, as calculated 
  !! for the quadrature, and calculates values of f0
  !! on that grid
  public :: calculate_f0_grids


  !> Initialises the module, chiefly reading the parameters.
  !! NB does not allocate arrays, as negrid must be provided
  !! by the egrid module.
  public :: init_general_f0

  !> Deallocate arrays and close output files. 
  public :: finish_general_f0

  !> Grid of f0 as function of energy and species
  public :: f0_grid
  real, dimension (:,:), allocatable :: f0_grid

  !> Grid of generalised temperature 
  !! = -1/T^*_sN d(F_s/d energy_s)_N
  !! = - T_r/T^*_s d(F_s/d energy_s) T*_s/F_s  
  !! (where T*_s is just temperature for Maxwellian species)
  !!  as function of energy and species.
  !! For Maxwellian species this quantity is just equal to T_r/T_s
  !! For alphas, T^*_s = E_alpha, the injection energy
  !! and this grid is calculated in this module

  public :: generalised_temperature
  real, dimension (:,:), allocatable :: generalised_temperature

  ! These arrays below replace the species quantites of the same name

  !> Equal to sqrt(generalised_temperature/mass) as a 
  !! function of energy and species
  public :: stm
  real, dimension (:,:), allocatable :: stm
  
  !> Generalized dF0/drho  for Maxwellian speices
  public :: f0prim
  real,dimension(:,:), allocatable:: f0prim

  ! get rid of this!
  !> Equal to Z/sqrt(generalised_temperature/mass) as a 
  !! function of energy and species
  public :: zstm
  real, dimension (:,:), allocatable :: zstm

  !> Equal to generalised_temperature/Z as a 
  !! function of energy and species
  public :: gtempoz
  real, dimension (:,:), allocatable :: gtempoz

  !> Equal to Z/generalised_temperature as a 
  !! function of energy and species
  public :: zogtemp
  real, dimension (:,:), allocatable :: zogtemp

  ! probably get rid of this!
  !> Equal to abs(sqrt(generalised_temperature/mass)/Z) as a 
  !! function of energy and species
  public :: smz
  real, dimension (:,:), allocatable :: smz

  !> Arrays are initially only calculated from proc0 as 
  !! calculate_f0_grids is called from setvgrid. Later
  !! this function is called to put them on all procs
  public :: broadcast_arrays

  private

  integer :: alpha_f0_switch, &
             beam_f0_switch


  integer, parameter :: alpha_f0_maxwellian = 1, &
                        alpha_f0_analytic = 2, &
                        alpha_f0_external =3

  integer, parameter :: beam_f0_maxwellian = 1, &
                        beam_f0_analytic = 2, &
                        beam_f0_external =3

  integer :: gen_f0_output_file
  
  integer :: negrid

  real, dimension(:,:), allocatable :: egrid

  real, dimension(:), allocatable :: egrid_maxwell

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module admin
!!
!! This subsection contains routines for initializing/finalizing 
!! the module and reading parameters. 

  subroutine init_general_f0
    use file_utils, only: open_output_file
    use mp, only: proc0
    logical, save :: initialized = .false.
    !write (*,*) "initializing parameter_scan"

    if (initialized) return
    initialized = .true.


    call read_parameters


    if (proc0) call open_output_file(gen_f0_output_file, ".general_f0")
      

  end subroutine init_general_f0



  subroutine finish_general_f0
    use file_utils, only: close_output_file

    call close_output_file(gen_f0_output_file)
  end subroutine finish_general_f0
  


  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (3), parameter :: alpha_f0_opts = &
         (/ text_option('maxwellian', alpha_f0_maxwellian), &
            text_option('analytic', alpha_f0_analytic), &
            text_option('external', alpha_f0_external) /)
    character(20) :: alpha_f0
    type (text_option), dimension (3), parameter :: beam_f0_opts = &
         (/ text_option('maxwellian', beam_f0_maxwellian), &
            text_option('analytic', beam_f0_analytic), &
            text_option('external', beam_f0_external) /)
    character(20) :: beam_f0
    namelist /general_f0_parameters/ &
            alpha_f0, &
            beam_f0

    integer :: ierr, in_file
    logical :: exist

    if (proc0) then

       alpha_f0 = 'maxwellian'
       beam_f0 = 'maxwellian'

       in_file = input_unit_exist ("general_f0_parameters", exist)
       if (exist) read (unit=in_file, nml=general_f0_parameters)

       ierr = error_unit()
       call get_option_value &
            (alpha_f0, alpha_f0_opts, alpha_f0_switch, &
            ierr, "alpha_f0 in general_f0_parameters")
       call get_option_value &
            (beam_f0, beam_f0_opts, beam_f0_switch, &
            ierr, "beam_f0 in general_f0_parameters")

    end if

    call broadcast (alpha_f0_switch)
    call broadcast (beam_f0_switch)

  end subroutine read_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General functions
!  Functions which apply to all species

  !> NB this only gets called by proc0
  !! Arrays are broadcast later
  subroutine calculate_f0_grids(epoints)
    use species, only: nspec, spec
    use species, only: ion_species, electron_species, alpha_species
    use species, only: beam_species
    real, dimension(:,:), intent(inout) :: epoints
    integer :: is
    negrid = size(epoints)
    call allocate_arrays
    egrid = epoints
    egrid_maxwell(:) = epoints(:,1)
    do is = 1,nspec
      select case (spec(is)%type)
      case (ion_species)
        call calculate_f0_grids_maxwellian(is)
      case (electron_species)
        call calculate_f0_grids_maxwellian(is)
      case (alpha_species)
        select case (alpha_f0_switch)
        case (alpha_f0_maxwellian)
          call calculate_f0_grids_maxwellian(is)
        end select
      case (beam_species)
        select case (beam_f0_switch)
        case (beam_f0_maxwellian)
          call calculate_f0_grids_maxwellian(is)
        end select
      end select
    end do
    epoints = egrid
    
  end subroutine calculate_f0_grids

  subroutine allocate_arrays
    use species, only: nspec

    allocate(egrid(negrid,nspec))
    allocate(egrid_maxwell(negrid))
    allocate(f0_grid(negrid,nspec))
    allocate(generalised_temperature(negrid,nspec))
    allocate(stm(negrid,nspec))
    allocate(zstm(negrid,nspec))
    allocate(zogtemp(negrid,nspec))
    allocate(gtempoz(negrid,nspec))
    allocate(smz(negrid,nspec))
    allocate(F0prim(negrid,nspec))

  end subroutine allocate_arrays

  subroutine broadcast_arrays
    use mp, only: proc0, broadcast
    
    call broadcast(negrid)
    if (.not. proc0) call allocate_arrays
    call broadcast(egrid)
    call broadcast(egrid_maxwell)
    call broadcast(f0_grid)
    call broadcast(generalised_temperature)
    call broadcast(stm)
    call broadcast(zstm)
    call broadcast(zogtemp)
    call broadcast(gtempoz)
    call broadcast(smz)
    call broadcast(F0prim)
  end subroutine broadcast_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maxwellian
!! 
!! This subsection contains routines for handling 
!! Maxwellian distributions. 


  subroutine calculate_f0_grids_maxwellian(is)
    use species, only: spec
    integer, intent(in) :: is
    !integer :: ie
    egrid(:,is) = egrid_maxwell(:)
    f0_grid(:, is) = exp(-egrid(:,is))
    !do ie = 1,negrid
      generalised_temperature(:,is) = spec(is)%temp
    !end do
    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
    
    f0prim(:,is) = -( spec(is)%fprim + (egrid(:,is) - 1.5)*spec(is)%tprim)
  end subroutine calculate_f0_grids_maxwellian


end module general_f0
