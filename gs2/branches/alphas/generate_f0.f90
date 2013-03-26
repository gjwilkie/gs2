
!> A module which allows the simulation of species with a
!! arbitrary background distribution function (e.g. alphas, 
!! beam ions).
!! 
!! Depending on the species type and the input parameters,
!! it will return either a Maxwellian or some other analytic
!! function or a solution for F_0 returned from an external
!! solver. 

module generate_f0


  !> Takes in the gridpoints in energy, as calculated 
  !! for the quadrature, and calculates values of f0
  !! on that grid
  public :: calculate_f0_grids


  !> Initialises the module, chiefly reading the parameters.
  !! NB does not allocate arrays, as negrid must be provided
  !! by the egrid module.
  public :: init_generate_f0

  !> Deallocate arrays and close output files. 
  public :: finish_generate_f0

  public :: f0_grid
  real, dimension (:,:), allocatable :: f0_grid

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

  real, dimension(:), allocatable :: egrid

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module admin
!!
!! This subsection contains routines for initializing/finalizing 
!! the module and reading parameters. 

  subroutine init_generate_f0
    use file_utils, only: open_output_file
    use mp, only: proc0
    logical, save :: initialized = .false.
    !write (*,*) "initializing parameter_scan"

    if (initialized) return
    initialized = .true.


    call read_parameters


    if (proc0) call open_output_file(gen_f0_output_file, ".generate_f0")
      

  end subroutine init_generate_f0



  subroutine finish_generate_f0
    use file_utils, only: close_output_file

    call close_output_file(gen_f0_output_file)
  end subroutine finish_generate_f0
  


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
    namelist /generate_f0_parameters/ &
            alpha_f0, &
            beam_f0

    integer :: ierr, in_file
    logical :: exist

    if (proc0) then

       alpha_f0 = 'maxwellian'
       beam_f0 = 'maxwellian'

       in_file = input_unit_exist ("generate_f0_parameters", exist)
       if (exist) read (unit=in_file, nml=generate_f0_parameters)

       ierr = error_unit()
       call get_option_value &
            (alpha_f0, alpha_f0_opts, alpha_f0_switch, &
            ierr, "alpha_f0 in generate_f0_parameters")
       call get_option_value &
            (beam_f0, beam_f0_opts, beam_f0_switch, &
            ierr, "beam_f0 in generate_f0_parameters")

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
    real, dimension(:), intent(in) :: epoints
    integer :: is
    negrid = size(epoints)
    allocate(egrid(negrid))
    egrid = epoints
    allocate (f0_grid(negrid,nspec))
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
    
  end subroutine calculate_f0_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maxwellian
!! 
!! This subsection contains routines for handling 
!! Maxwellian distributions. 


  subroutine calculate_f0_grids_maxwellian(is)
    integer, intent(in) :: is
    f0_grid(is, :) = exp(-egrid)
  end subroutine calculate_f0_grids_maxwellian



end module generate_f0
