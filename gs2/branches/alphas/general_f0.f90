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
  !! for the Maxwellian species (which will be in the first
  !! column of the input array, epoints(:,1)), 
  !! calculates any other grids required
  !! and calculates values of f0
  !! and all other arrays on those grids
  public :: calculate_f0_arrays

  !> A single local function that calculates f0. Takes one argument
  public :: eval_f0

  !> To facilitate the need to only call one function with one argument
  !! in eval_f0, user must first call this function to set the species
  public :: set_current_f0_species


  !> Initialises the module, chiefly reading the parameters.
  !! NB does not allocate arrays, as negrid must be provided
  !! by the egrid module.
  public :: init_general_f0

  !> Deallocate arrays and close output files. 
  public :: finish_general_f0

  !> Grid of f0 as function of energy and species
  public :: f0_values
  real, dimension (:,:), allocatable :: f0_values

  !> Grid of generalised temperature 
  !! = -1/T^*_sN d(F_s/d energy_s)_N
  !! = - T_r/T^*_s d(F_s/d energy_s) T*_s/F_s  
  !! (where T*_s is just temperature for Maxwellian species)
  !!  as function of energy and species.
  !! For Maxwellian species this grid is just equal to T_r/T_s
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
  !! calculate_f0_arrays is called from setvgrid. Later
  !! this function is called to put them on all procs
  public :: broadcast_arrays

  !> Unit tests
  public :: general_f0_unit_test_init_general_f0
  public :: general_f0_unit_test_calculate_f0_arrays

  private

  logical :: genquad_flag = .false.

  integer :: is_local = -1

  integer :: alpha_f0_switch, &
             beam_f0_switch
  
  logical :: print_egrid

  !> Flag that controls whether or not an externally-supplied f0 
  !! is rescaled to fit the given species parameters.
  !! f0_rescale = T -- F0 is rescaled to fit spec(is)%dens
  !!            = F -- F0 is taken literally, spec(is)%dens is changed                     
  !! This option does not rescale according to spec(is)%temp!
  logical :: rescale_f0

  integer, parameter :: alpha_f0_maxwellian = 1, &
                        alpha_f0_analytic = 2, &
                        alpha_f0_semianalytic = 3, &
                        alpha_f0_split = 4, &
                        alpha_f0_equivmaxw = 5, &
                        alpha_f0_external =6, &
                        alpha_f0_kappa = 7

  integer, parameter :: beam_f0_maxwellian = 1, &
                        beam_f0_analytic = 2, &
                        beam_f0_external = 3

  integer :: gen_f0_output_file
  
  integer :: negrid

  real :: vcut
  real :: energy_min

  real, dimension(:,:), allocatable :: egrid
  real, dimension(:,:), allocatable :: weights

  real, dimension(:), allocatable :: egrid_maxwell
  real, dimension(:), allocatable :: weights_maxwell

  !> Which species to use as the main ions in falpha
  integer :: main_ion_species, vcrit_opt

  !> "Critical speed" that parameterizes the analytic slowing-down distribution 
  !! Can, in principle, be calculated from ion and electron properties, but for now
  !! just an input parameter
  real:: vcrit 

  !> Input parameter that gives a simple gradient in F0alpha (not a function of velocity - simply 1/Lnalpha)
  !! Not physical, but could diagnose problems.
  logical :: simple_gradient

  real :: kappa

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


  function general_f0_unit_test_init_general_f0()
    use unit_tests
    use species
    logical :: general_f0_unit_test_init_general_f0

    call init_general_f0

    call announce_check('alpha_f0_switch')
    general_f0_unit_test_init_general_f0 = &
      agrees_with(alpha_f0_switch, alpha_f0_split)
    call process_check(general_f0_unit_test_init_general_f0, 'alpha_f0_switch')

    call announce_check('spec 1 type')
    general_f0_unit_test_init_general_f0 = &
      general_f0_unit_test_init_general_f0 .and. &
      agrees_with(spec(1)%type, ion_species)
    call process_check(general_f0_unit_test_init_general_f0, 'spec 1 type')
    call announce_check('spec 2 type')
    general_f0_unit_test_init_general_f0 = &
      general_f0_unit_test_init_general_f0 .and. &
      agrees_with(spec(2)%type, electron_species)
    call process_check(general_f0_unit_test_init_general_f0, 'spec 2 type')
    call announce_check('spec 3 type')
    general_f0_unit_test_init_general_f0 = &
      general_f0_unit_test_init_general_f0 .and. &
      agrees_with(spec(3)%type, alpha_species)
    call process_check(general_f0_unit_test_init_general_f0, 'spec 3 type')
    
  end function general_f0_unit_test_init_general_f0

  subroutine finish_general_f0
    use file_utils, only: close_output_file

    call close_output_file(gen_f0_output_file)
  end subroutine finish_general_f0

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type (text_option), dimension (7), parameter :: alpha_f0_opts = &
         (/ text_option('maxwellian', alpha_f0_maxwellian), &
            text_option('analytic', alpha_f0_analytic), &
            text_option('semianalytic', alpha_f0_semianalytic), &
            text_option('split', alpha_f0_split), &
            text_option('equivmaxw', alpha_f0_equivmaxw), &
            text_option('external', alpha_f0_external) ,&
            text_option('kappa', alpha_f0_kappa) /)
    character(20) :: alpha_f0
    type (text_option), dimension (3), parameter :: beam_f0_opts = &
         (/ text_option('maxwellian', beam_f0_maxwellian), &
            text_option('analytic', beam_f0_analytic), &
            text_option('external', beam_f0_external) /)
    character(20) :: beam_f0
    namelist /general_f0_parameters/ &
            alpha_f0, &
            beam_f0, &
            rescale_f0,&
            main_ion_species,&
            energy_min,&
            print_egrid,&
            vcrit, simple_gradient, kappa, vcrit_opt

    integer :: ierr, in_file
    logical :: exist

    if (proc0) then

       alpha_f0 = 'maxwellian'
       beam_f0 = 'maxwellian'
       main_ion_species = -1 
       energy_min = 0.1
       vcrit = -1.0
       simple_gradient = .false.
       kappa = -1.0
       vcrit_opt = 1

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
    call broadcast (main_ion_species)
    call broadcast (energy_min)
    call broadcast (vcrit)
    call broadcast (kappa)
    call broadcast (vcrit_opt)

  end subroutine read_parameters

  subroutine set_current_f0_species(is)
    implicit none
    integer,intent(in):: is
    is_local = is
  end subroutine

  real function eval_f0(v)
    use mp, only: mp_abort
    use species, only: ion_species, electron_species, alpha_species,beam_species,spec
    use constants, only: pi
    implicit none
    real,intent(in):: v
    integer:: is    
    real:: f0, dummy
  
    if ( is_local .LT. 0 ) then 
       write(*,*) "ERROR: Species needs to be set with set_current_f0_species to call eval_f0" 
       call mp_abort('')
    end if

    is = is_local

    select case (spec(is)%type)
      case (ion_species)
        f0 = exp(-v**2)/(pi**1.5)
      case (electron_species)
        f0 = exp(-v**2)/(pi**1.5)
      case (alpha_species)
        select case (alpha_f0_switch)
        case (alpha_f0_maxwellian)
          f0 = exp(-v**2)/(pi**1.5)
        case (alpha_f0_analytic)
          call eval_f0_analytic(is,v,f0,dummy,dummy) 
        case (alpha_f0_kappa)
          call eval_f0_kappa(is,v,f0,dummy,dummy) 
        case (alpha_f0_semianalytic)
          write(*,*) "ERROR: eval_f0 cannot yet handle semianalytic option."
          call mp_abort('')
        case (alpha_f0_split)
          write(*,*) "ERROR: eval_f0 cannot yet handle split option."
          call mp_abort('')
        case (alpha_f0_equivmaxw)
          f0 = exp(-v**2)/(pi**1.5)
        case (alpha_f0_external)
          write(*,*) "ERROR: eval_f0 cannot yet handle external option."
          call mp_abort('')
        end select
      case (beam_species)
        select case (beam_f0_switch)
        case (beam_f0_maxwellian)
          call calculate_f0_arrays_maxwellian(is)
        end select
    end select

    eval_f0 = f0
    return
  end function eval_f0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General functions
!  Functions which apply to all species

  !> NB this only gets called by proc0
  !! Arrays are broadcast later
  subroutine calculate_f0_arrays(epoints, wgts, vcut_in, genquad_flag_in)
    use species, only: nspec, spec
    use species, only: ion_species, electron_species, alpha_species
    use species, only: beam_species
    real, dimension(:,:), intent(inout) :: epoints
    real, dimension(:,:), intent(inout) :: wgts
    real, intent(in) :: vcut_in
    logical,intent(in):: genquad_flag_in
    integer :: is
    genquad_flag = genquad_flag_in
    negrid = size(epoints(:,1))
    vcut = vcut_in
    call allocate_arrays
    egrid = epoints
    weights = wgts
    egrid_maxwell(:) = epoints(:,1)
    weights_maxwell(:) = wgts(:,1)
    do is = 1,nspec
      select case (spec(is)%type)
      case (ion_species)
        call calculate_f0_arrays_maxwellian(is)
      case (electron_species)
        call calculate_f0_arrays_maxwellian(is)
      case (alpha_species)
        select case (alpha_f0_switch)
        case (alpha_f0_maxwellian)
          call calculate_f0_arrays_maxwellian(is)
        case (alpha_f0_analytic)
          call check_electromagnetic
          call calculate_f0_arrays_analytic(is)
        case (alpha_f0_kappa)
          call check_electromagnetic
          call calculate_f0_arrays_kappa(is)
        case (alpha_f0_semianalytic)
          call check_electromagnetic
          call calculate_f0_arrays_semianalytic(is)
        case (alpha_f0_split)
          call check_electromagnetic
          call calculate_f0_arrays_split(is)
        case (alpha_f0_equivmaxw)
          call calculate_f0_arrays_equivmaxw(is)
        case (alpha_f0_external)
          call check_electromagnetic
          call calculate_f0_arrays_external(is)
        end select
      case (beam_species)
        select case (beam_f0_switch)
        case (beam_f0_maxwellian)
          call calculate_f0_arrays_maxwellian(is)
        end select
      end select
    end do
    epoints = egrid
    wgts = weights
    
  end subroutine calculate_f0_arrays

  subroutine check_electromagnetic
    use mp, only: mp_abort
    !use run_parameters, only: fapar, fbpar
    !if (abs(fapar) > epsilon(0.0) .or. abs(fbpar) > epsilon(0.0)) then
      !write(*,*) &
      !'Non-maxwellian species not implemented for electromagnetic runs'
      !call mp_abort('')
    !end if
  end subroutine check_electromagnetic


  function general_f0_unit_test_calculate_f0_arrays(epoints, wgts, vcut_in, rslts, err)
    use unit_tests
    real, dimension(:,:), intent(inout) :: epoints
    real, dimension(:,:), intent(inout) :: wgts
    real, dimension(:,:,:), intent(in) :: rslts
    real, intent(in) :: err
    real, intent(in) :: vcut_in
    logical :: general_f0_unit_test_calculate_f0_arrays
    logical :: tr ! test results
    tr = .true.
    call calculate_f0_arrays(epoints, wgts, vcut_in,genquad_flag)

    call announce_check('ion energy grid')
    tr = tr .and. agrees_with(egrid(:,1), epoints(:,1), err)
    call process_check(tr, 'ion energy grid')

    call announce_check('electron energy grid')
    tr = tr .and. agrees_with(egrid(:,2), epoints(:,1), err)
    call process_check(tr, 'electron energy grid')

    call announce_check('ion f0')
    tr = tr .and. agrees_with(f0_values(:,1), rslts(:,1,1), err)
    call process_check(tr, 'ion f0')
    call announce_check('electron f0')
    tr = tr .and. agrees_with(f0_values(:,2), rslts(:,2,1), err)
    call process_check(tr, 'electron f0')
    call announce_check('alpha f0')
    tr = tr .and. agrees_with(f0_values(:,3), rslts(:,3,1), err)
    call process_check(tr, 'alpha f0')

    call announce_check('ion gentemp')
    tr = tr .and. agrees_with(generalised_temperature(:,1), rslts(:,1,2), err)
    call process_check(tr,' ion gentemp')
    call announce_check('electron gentemp')
    tr = tr .and. agrees_with(generalised_temperature(:,2), rslts(:,2,2), err)
    call process_check(tr,' electron gentemp')
    call announce_check('alpha gentemp')
    tr = tr .and. agrees_with(generalised_temperature(:,3), rslts(:,3,2), err)
    call process_check(tr,' alpha gentemp')

    !write (*,*) 'f0prim(:,3),', f0prim(:,3)
    call announce_check('ion f0prim')
    !tr = tr .and. agrees_with(f0prim(:,1), rslts(:,1,3), err)
    call process_check(tr,' ion f0prim')
    call announce_check('electron f0prim')
    !tr = tr .and. agrees_with(f0prim(:,2), rslts(:,2,3), err)
    call process_check(tr,' electron f0prim')
    call announce_check('alpha f0prim')
    tr = tr .and. agrees_with(f0prim(:,3), rslts(:,3,3), err*10.0)
    call process_check(tr,' alpha f0prim')
    !write (*,*) 'f0prim(:,3),', f0prim(:,3)


    general_f0_unit_test_calculate_f0_arrays = tr

  end function general_f0_unit_test_calculate_f0_arrays

  subroutine allocate_arrays
    use species, only: nspec

    allocate(egrid(negrid,nspec))
    allocate(egrid_maxwell(negrid))
    allocate(weights(negrid,nspec))
    allocate(weights_maxwell(negrid))
    allocate(f0_values(negrid,nspec))
    allocate(generalised_temperature(negrid,nspec))
    allocate(stm(negrid,nspec))
    allocate(zstm(negrid,nspec))
    allocate(zogtemp(negrid,nspec))
    allocate(gtempoz(negrid,nspec))
    allocate(smz(negrid,nspec))
    allocate(f0prim(negrid,nspec))

  end subroutine allocate_arrays

  subroutine broadcast_arrays
    use mp, only: proc0, broadcast
    
    call broadcast(negrid)
    if (.not. proc0) call allocate_arrays
    call broadcast(egrid)
    call broadcast(egrid_maxwell)
    call broadcast(weights)
    call broadcast(weights_maxwell)
    call broadcast(f0_values)
    call broadcast(generalised_temperature)
    call broadcast(stm)
    call broadcast(zstm)
    call broadcast(zogtemp)
    call broadcast(gtempoz)
    call broadcast(smz)
    call broadcast(f0prim)
  end subroutine broadcast_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maxwellian
!! 
!! This subsection contains routines for handling 
!! Maxwellian distributions. 


  subroutine calculate_f0_arrays_maxwellian(is)
    use constants, only: pi
    use species, only: spec
    integer, intent(in) :: is
    integer :: ie
    ! Energy grid now set for every species in module egrid
    !egrid(:,is) = egrid_maxwell(:)
    f0_values(:, is) = exp(-egrid(:,is))/(pi**1.5)
    !weights(:,is) = weights_maxwell(:) * f0_values(:,is)
    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
    do ie = 1,negrid
       generalised_temperature(:,is) = spec(is)%temp
       if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), & 
                                   generalised_temperature(ie,is)
    end do
    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
    
    f0prim(:,is) = -( spec(is)%fprim + (egrid(:,is) - 1.5)*spec(is)%tprim)
  end subroutine calculate_f0_arrays_maxwellian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Semianalytic
!! 
!! This subsection contains routines for handling 
!! distributions given by the semianalytical form,
!! ammended such that all alphas, including ash,
!! are combined in one species.

  subroutine calculate_f0_arrays_semianalytic(is)
    use semianalytic_falpha, only: semianalytic_falpha_parameters_type
    use semianalytic_falpha, only: calculate_arrays
    use constants, only: pi
    use species, only: spec 
    use species, only: nspec
    use species, only: electron_species
    use species, only: ion_species
    use species, only: has_electron_species
    use mp, only: mp_abort, proc0
    type(semianalytic_falpha_parameters_type) :: parameters
    integer, intent(in) :: is
    integer :: ie
    integer :: electron_spec
    integer :: i
    real :: shift, scal ! For modifying the energy grid for alphas

    !> Determine which species is the electrons
    !main_ion_species = -1
    do i = 1,nspec
      if (spec(i)%type .eq. electron_species) electron_spec = i
      if (main_ion_species < 1 .and. spec(i)%type .eq. ion_species) &
        main_ion_species = i
    end do

    write(*,*) 'main_ion_species', main_ion_species, 'electron_spec', electron_spec

    !if (main_ion_species < 1)  call mp_abort(&
      !'main_ion_species < 1: please set main_ion_species in general_f0_parameters') 

    parameters%source_prim = spec(is)%sprim
    
    parameters%alpha_density = spec(is)%dens      

    parameters%alpha_is = is

    parameters%alpha_ion_collision_rate = spec(is)%gamma_ai
    parameters%alpha_electron_collision_rate = spec(is)%gamma_ae

    parameters%alpha_vth       = spec(is)%stm

    parameters%source          = spec(is)%source
    parameters%ash_fraction    = spec(is)%ash_fraction
    parameters%source_prim     = spec(is)%sprim
    !write (*,*) 'Source prim is', parameters%source_prim

    parameters%alpha_injection_energy = spec(is)%temp

    parameters%alpha_charge    = spec(is)%z

    parameters%alpha_mass      = spec(is)%mass

    parameters%negrid = negrid
    if (has_electron_species(spec)) then 
      parameters%electron_fprim   = spec(electron_spec)%fprim
      parameters%electron_tprim   = spec(electron_spec)%tprim
      parameters%electron_mass   = spec(electron_spec)%mass
      parameters%electron_charge = spec(electron_spec)%z
      parameters%electron_vth    = spec(electron_species)%stm
      parameters%electron_temp   = spec(electron_species)%temp
      if (main_ion_species>0) then 
        parameters%ash_temp        = spec(main_ion_species)%temp
        parameters%ion_mass        = spec(main_ion_species)%mass
        parameters%ion_tprim        = spec(main_ion_species)%tprim
        parameters%ion_fprim        = spec(main_ion_species)%fprim
        parameters%ion_charge      = spec(main_ion_species)%z
        parameters%ion_vth         = spec(main_ion_species)%stm
        parameters%ion_temp        = spec(main_ion_species)%temp
      else
        ! Assume main ions are deuterium with ti = te
        parameters%ion_mass       = spec(electron_spec)%mass * 2.0 * 1836.0
        parameters%ion_tprim      = 0.0
        parameters%ion_fprim      = spec(electron_spec)%fprim * 0.0
        parameters%ion_charge     = -spec(electron_spec)%z
        parameters%ion_temp       = spec(electron_spec)%temp
        parameters%ion_vth        = sqrt(parameters%ion_temp/parameters%ion_mass)
      end if


    else
      if (main_ion_species>0) then 

        parameters%electron_fprim   =  0.0
        parameters%electron_tprim   =  0.0
        ! Assume main ions are deuterium with ti = te
        parameters%electron_mass   = spec(main_ion_species)%mass/1836.0/2.0
        parameters%electron_charge = -spec(main_ion_species)%z
        parameters%electron_temp   = spec(main_ion_species)%temp
        parameters%electron_vth    = sqrt(parameters%electron_temp/parameters%electron_mass)
      else
        write (*,*) 'You have no electrons and no ions: not sure how to set'
        write (*,*) 'alpha parameters!'
        call mp_abort('')
      end if
    end if

    !if (parameters%source .eq. 0.0) then 
      !write (*,*) 'You have source = 0.0 for alphas!!!'
      !call mp_abort(' ')
    !end if

    ! This bit of code shifts the grid so that vcut for alphas is 1
    ! and the lowest point is energy min
    !scal =  (energy_min - 1.0) / (egrid(1,1) - vcut)
    !shift = energy_min - egrid(1,1) * scal
    !egrid(:,is) = egrid_maxwell(:) *(1.0 - energy_min)/vcut + energy_min
    !egrid(:,is) = egrid_maxwell(:) *scal + shift

    !write (*,*) 'parameters', parameters
    !write(*,*) 'calling calculate_arrays'
    !write (*,*) 'egrid', egrid(:,is)
    call calculate_arrays(parameters,&
                          egrid, &
                          f0_values, &
                          generalised_temperature, &
                          f0prim)
    spec(is)%source = parameters%source
    !write (*,*) 'f0prim', ',is', f0prim(:,is), is
    !write (*,*) 'f0prim(:,3),', f0prim(:,3)
    !write (*,*) 'f0_values(:,3),', f0_values(:,3)
    !weights(:,is) =  weights_maxwell(:)*scal * f0_values(:,is)
    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)

    if (proc0) write(*,*) "source = ", parameters%source
    if (proc0) write(*,*) "ash_fraction = ", parameters%ash_fraction
!    do ie = 1,negrid
!       if (proc0) write(*,*) sqrt(egrid(ie,is)), f0_values(ie,is), generalised_temperature(ie,is)
!    end do
    
  end subroutine calculate_f0_arrays_semianalytic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Split
!! 
!! This subsection contains routines for handling 
!! distributions given by the semianalytical form
!! found with fast alphas as a seperate species 
!! from ash.

  subroutine calculate_f0_arrays_split(is)
    use split_falpha, only: split_falpha_parameters_type
    use split_falpha, only: calculate_arrays
    use constants, only: pi
    use species, only: spec 
    use species, only: nspec
    use species, only: electron_species
    use species, only: ion_species
    use species, only: has_electron_species
    use mp, only: mp_abort
    type(split_falpha_parameters_type) :: parameters
    integer, intent(in) :: is
    integer :: ie
    integer :: electron_spec
    integer :: i
    real :: shift, scal ! For modifying the energy grid for alphas

    !> Determine which species is the electrons
    !main_ion_species = -1
    do i = 1,nspec
      if (spec(i)%type .eq. electron_species) electron_spec = i
      if (main_ion_species < 1 .and. spec(i)%type .eq. ion_species) &
        main_ion_species = i
    end do

    write(*,*) 'main_ion_species', main_ion_species, 'electron_spec', electron_spec

    !if (main_ion_species < 1)  call mp_abort(&
      !'main_ion_species < 1: please set main_ion_species in general_f0_parameters') 

    parameters%source_prim = spec(is)%sprim
    
    parameters%alpha_density = spec(is)%dens

    parameters%alpha_is = is

    parameters%alpha_ion_collision_rate = spec(is)%gamma_ai
    parameters%alpha_electron_collision_rate = spec(is)%gamma_ae

    parameters%alpha_vth       = spec(is)%stm

    parameters%source          = spec(is)%source
    parameters%source_prim     = spec(is)%sprim
    !write (*,*) 'Source prim is', parameters%source_prim

    parameters%alpha_injection_energy = spec(is)%temp

    parameters%alpha_charge    = spec(is)%z

    parameters%alpha_mass      = spec(is)%mass

    parameters%negrid = negrid
    if (has_electron_species(spec)) then 
      parameters%electron_fprim   = spec(electron_spec)%fprim
      parameters%electron_tprim   = spec(electron_spec)%tprim
      parameters%electron_mass   = spec(electron_spec)%mass
      parameters%electron_charge = spec(electron_spec)%z
      parameters%electron_vth    = spec(electron_species)%stm
      parameters%electron_temp   = spec(electron_species)%temp
      if (main_ion_species>0) then 
        parameters%ion_mass        = spec(main_ion_species)%mass
        parameters%ion_tprim        = spec(main_ion_species)%tprim
        parameters%ion_fprim        = spec(main_ion_species)%fprim
        parameters%ion_charge      = spec(main_ion_species)%z
        parameters%ion_vth         = spec(main_ion_species)%stm
        parameters%ion_temp        = spec(main_ion_species)%temp
      else
        ! Assume main ions are deuterium with ti = te
        parameters%ion_mass       = spec(electron_spec)%mass * 2.0 * 1836.0
        parameters%ion_tprim      = 0.0
        parameters%ion_fprim      = spec(electron_spec)%fprim * 0.0
        parameters%ion_charge     = -spec(electron_spec)%z
        parameters%ion_temp       = spec(electron_spec)%temp
        parameters%ion_vth        = sqrt(parameters%ion_temp/parameters%ion_mass)
      end if


    else
      if (main_ion_species>0) then 

        parameters%electron_fprim   =  0.0
        parameters%electron_tprim   =  0.0
        ! Assume main ions are deuterium with ti = te
        parameters%electron_mass   = spec(main_ion_species)%mass/1836.0/2.0
        parameters%electron_charge = -spec(main_ion_species)%z
        parameters%electron_temp   = spec(main_ion_species)%temp
        parameters%electron_vth    = sqrt(parameters%electron_temp/parameters%electron_mass)
      else
        write (*,*) 'You have no electrons and no ions: not sure how to set'
        write (*,*) 'alpha parameters!'
        call mp_abort('')
      end if
    end if

    !if (parameters%source .eq. 0.0) then 
      !write (*,*) 'You have source = 0.0 for alphas!!!'
      !call mp_abort(' ')
    !end if

    ! This bit of code shifts the grid so that vcut for alphas is 1
    ! and the lowest point is energy min
    !scal =  (energy_min - 1.0) / (egrid(1,1) - vcut)
    !shift = energy_min - egrid(1,1) * scal
    !egrid(:,is) = egrid_maxwell(:) *(1.0 - energy_min)/vcut + energy_min
    !egrid(:,is) = egrid_maxwell(:) *scal + shift

    !write (*,*) 'parameters', parameters
    !write(*,*) 'calling calculate_arrays'
    !write (*,*) 'egrid', egrid(:,is)
    call calculate_arrays(parameters,&
                          egrid, &
                          f0_values, &
                          generalised_temperature, &
                          f0prim)
    spec(is)%source = parameters%source
    !write (*,*) 'f0prim', ',is', f0prim(:,is), is
    !write (*,*) 'f0prim(:,3),', f0prim(:,3)
    !write (*,*) 'f0_values(:,3),', f0_values(:,3)
    !weights(:,is) =  weights_maxwell(:)*scal * f0_values(:,is)
    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
    
  end subroutine calculate_f0_arrays_split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! External
!! 
!! This subsection contains routines for handling 
!! distributions that are given by an external
!! input file.

!> This subroutine calculates f0 on the grid from an external
!! input file. The grid on the input file can differ from that 
!! of gs2. A cubic-spline is used to interpolate between the two.
!! This is a logarithmic interpolation on the splines to ensure positive-
!! definiteness of F0. The user can either specify df0/dE or it 
!! can be estimated internally.
  subroutine calculate_f0_arrays_external(is)
    !
    ! The egrid is already given by get_legendre_grids. The input values
    ! are cubic-spline interpolated to the energy grid of gs2.
    !
    ! The input file can take two forms, as controlled by the
    ! control paramter num_cols, the first integer 
    ! read from the file
    !  - Single-Column Mode: Only F0(E) is given. To calculate
    !      generalized temperature, a cubic spline is used to 
    !      interpolate and the slope is taken from that.
    !  - Two-Column Mode: First column is F0(E), the second is
    !      dF0/dE.
    ! 
    use mp, only: broadcast
    use constants, only: pi
    use species, only: spec, nspec
    use file_utils, only: run_name
    use splines, only: fitp_curvd, fitp_curv1, fitp_curv2
    implicit none
    integer, intent(in) :: is
    integer:: f0in_unit = 21, num_cols, ie, ierr, numdat, il, it
    real:: df0dE, test, moment0, moment2
    real:: pick_spec(nspec)
    real, dimension(:), allocatable:: f0_values_dat, df0dE_dat, egrid_dat, &
                                      f0_values_dat_log, df0dE_dat_log, yp, temp
    
    ! Open file and read column option
    open(unit=f0in_unit,file=trim(run_name)//'.f0in',status='old',action='read')
    read(f0in_unit,*) num_cols
    read(f0in_unit,*) numdat

    allocate(f0_values_dat(numdat))
    allocate(df0dE_dat(numdat))
    allocate(egrid_dat(numdat))
    allocate(f0_values_dat_log(numdat))
    allocate(df0dE_dat_log(numdat))
    allocate(yp(numdat))
    allocate(temp(numdat))

    if (num_cols .EQ. 2) then
       ! Read f0 values
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie)
          ! Interpolate the *logarithm* of f0 to ensure positive-definiteness
          f0_values_dat_log(ie) = log(f0_values_dat(ie))
       end do
       close(f0in_unit)

       ! Perform cubic spline to get F0 and its slope at grid points

       ! Generate spline parameters
       call fitp_curv1(numdat,egrid_dat,f0_values_dat_log,0.0,0.0,3,yp,temp,1.0,ierr)
       if (ierr .NE. 0) then
          write(*,*) "fitp_curv1 returned error code ", ierr
          stop 1
       end if

       do ie = 1,negrid
          ! Interpolate to get f0 at grid points
          f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat_log,yp,1.0)

          ! Recover F0 from its logarithm
          f0_values(ie,is) = exp(f0_values(ie,is))

          ! Calculate d/dE lnF0 to get generalised temperature
          df0dE = fitp_curvd(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat_log,yp,1.0)
          generalised_temperature(ie,is) = - spec(is)%temp/df0dE

          ! Diagnostic output
          if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), df0dE, & 
                                   generalised_temperature(ie,is)
       end do

    else if (num_cols .EQ. 3) then
       ! Read both f0 and df0/dE
       do ie = 1,numdat
          read(f0in_unit,*) egrid_dat(ie), f0_values_dat(ie), df0dE_dat(ie)
          f0_values_dat_log(ie) = log(f0_values_dat(ie))
          df0dE_dat_log(ie) = log(abs(df0dE_dat(ie)))
       end do
       close(f0in_unit)

       ! Generate spline parameters for f0
       call fitp_curv1(numdat,egrid_dat,f0_values_dat_log,0.0,0.0,3,yp,temp,1.0,ierr)

       if (ierr .NE. 0) then
          write(*,*) "fitp_curv1 returned error code ", ierr
          stop 1
       end if

       do ie = 1,negrid

          ! Interpolate to get F0 at grid points
          f0_values(ie,is) = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             f0_values_dat_log,yp,1.0)
       end do

       ! Recover F0 from its logarithm
       f0_values = exp(f0_values)

       ! Generate spline parameters for df0/dE
       call fitp_curv1(numdat,egrid_dat,df0dE_dat_log,0.0,0.0,3,yp,temp,1.0,ierr)
       if (ierr .NE. 0) then
          write(*,*) "fitp_curv1 returned error code ", ierr
          stop 1
       end if

       do ie = 1,negrid
          ! Interpolate to get f0 at grid points
          df0dE            = fitp_curv2(egrid(ie,is),numdat,egrid_dat, &
                             df0dE_dat_log,yp,1.0)
 
          ! Recover df0/dE from its logarithm (maintaining whatever sign it had before)
          df0dE = sign(1.0,df0dE_dat(ie))* exp(df0dE)

          generalised_temperature(ie,is) = -spec(is)%temp*f0_values(ie,is)/df0dE

          ! Diagnostic output
          if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), df0dE, & 
                                   generalised_temperature(ie,is)

       end do
    else
       write(*,*) "ERROR. First line in f0 input file should be num_cols: " 
       write(*,*) " num_cols = 1 if only f0 is to be input, "
       write(*,*) " num_cols=2 if f0 and df0/dE are input."
       stop 1
    end if

    ! Now calculate moments of F0 from trapezoidal rule.
    ! Should probably do this by quadrature, 
    ! but le_grids is not defined yet

    moment0 = 0.0
    do ie = 1,negrid-1
       moment0  = moment0 + 0.5*(egrid(ie+1,is)-egrid(ie,is)) * &
                           ( sqrt(egrid(ie,is))*f0_values(ie,is) + &
                             sqrt(egrid(ie+1,is))*f0_values(ie+1,is) )
    end do
!    do ie = 1,numdat-1
!       moment0  = moment0 + 0.5*(egrid_dat(ie+1)-egrid_dat(ie)) * &
!                           ( sqrt(egrid_dat(ie))*f0_values_dat(ie) + &
!                             sqrt(egrid_dat(ie+1))*f0_values_dat(ie+1) )
!    end do
    moment0 = moment0*4.0*pi

    ! Input parameter rescale_f0 determines the priority between the input species
    ! density or the input F0.
    if (rescale_f0) then
       write(*,*) "rescale_f0 = T: Rescaling magnitude of F0 to agree with input dens"
       write(*,*) "rescale factor is ", spec(is)%dens / moment0
       f0_values = f0_values * spec(is)%dens / moment0
    else
       write(*,*) "rescale_f0 = F: Overriding species dens to agree with input F0"
       write(*,*) "New dens = ", moment0
       spec(is)%dens = moment0
    end if

    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
    
    f0prim(:,is) = -( spec(is)%fprim + (egrid(:,is) - 1.5)*spec(is)%tprim)

    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)

    deallocate(egrid_dat,f0_values_dat,df0dE_dat,f0_values_dat_log,df0dE_dat_log,yp,temp)

  end subroutine calculate_f0_arrays_external

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Analytic
!! 
!! This subsection contains routines for handling 
!! the analytic approximation to the fast-particle 
!! slowing-down distribution 
!!
!! f0(v) = A/(vc^3 + v^3)                                           (if v < vstar)
!!       = B*exp(-Ealpha*v^2/Ti)                                    (if v > vstar)
!! 
!! vc    = critical speed, input
!! A     = normalization chosen so that f0 integrates to unity
!! B     = constant chosen to enforce continuity at v = vstar
!!
!! To do: 
!! - Calculate vc based on parameters of other species.
!! - Implement an appropriate generalization of temperature gradient, perhaps based on other species
!!
  subroutine calculate_f0_arrays_analytic(is)
    use constants, only: pi
    use species, only: spec, nspec
    use mp, only: mp_abort, proc0
    use species, only: ion_species, electron_species, alpha_species
    implicit none
    integer,intent(in):: is
    real:: A, Ti, Ealpha, vstar, v, df0dv, df0dE, B
    integer:: ie, i, electron_spec

!    if (vcrit .LE. 0.0) then
!       write(*,*) "ERROR: in general_f0. If alpha_f0='analytic' is chosen, vcrit must also be specified."
!       call mp_abort('')
!    end if

    do ie = 1,negrid 
       v = sqrt(egrid(ie,is))
       call eval_f0_analytic(is,v,f0_values(ie,is),generalised_temperature(ie,is),f0prim(ie,is))
    end do

    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)

    if (simple_gradient) f0prim(:,is) = - spec(is)%fprim 

    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)

  end subroutine calculate_f0_arrays_analytic

  subroutine eval_f0_analytic(is,v,f0,gentemp,f0prim)
    use species, only: ion_species, electron_species, alpha_species, spec, nspec, ZI_fac, vte
    use constants, only: pi
    use mp, only: mp_abort, proc0
    use spfunc, only: erf => erf_ext
    implicit none
    integer,intent(in):: is
    real,intent(in):: v
    real,intent(inout):: f0,gentemp, f0prim
    real:: A, Ti, Ealpha, vstar, df0dv, df0dE, vta, vti, Zi,ni,ne, ni_prim,ne_prim,Ti_prim,Te_prim
    integer:: ie, i, electron_spec

    electron_spec = -1
    do i = 1,nspec
      if (spec(i)%type .eq. electron_species) electron_spec = i
      if (main_ion_species < 1 .and. spec(i)%type .eq. ion_species) &
        main_ion_species = i
    end do

    Ealpha = spec(is)%temp                    !< temp for alpha species is interpreted as normalized injection energy
    vta = sqrt(Ealpha/spec(is)%mass)

       if (electron_spec .GT. 0) then
          ne = spec(electron_spec)%dens
          ne_prim = spec(electron_spec)%fprim
          Te_prim = spec(electron_spec)%tprim
       else
          ne = sum(spec(:)%z * spec(:)%dens )
          ne_prim = sum(spec(:)%z * spec(:)%dens * spec(:)%fprim )/ne
          Te_prim = spec(main_ion_species)%tprim
       end if

    if (vcrit .LE. 0.0) then
       !> Note that the ratio of me/ma is hardcoded. 
       vcrit = (0.25*3.0*sqrt(pi)*ZI_fac*1.371e-4)**(1.0/3.0)*vte/vta
    end if

    A = (4.0*pi/3.0)*log( (vcrit**3 + 1.0)/(vcrit**3))
    A = 1.0/A

    if (v .LE. 1.0) then
       f0 = A/(vcrit**3 + v**3)
       df0dv = -A*3.0*v**2/(vcrit**3 + v**3)**2
       df0dE = (0.5/v)*df0dv
       gentemp = -spec(is)%temp*f0/df0dE
       if (vcrit_opt .EQ. 1) then
          f0prim = -spec(is)%fprim - (-(vcrit**3/(vcrit**3 + v**3)) + (1.0/( log(1.0 + vcrit**(-3)) * (1.0 + vcrit**(3)))))*(ni_prim - ne_prim + Ti_prim + 0.5*Te_prim)
       else if (vcrit_opt .EQ. 2) then
          f0prim = -spec(is)%fprim - (-(vcrit**3/(vcrit**3 + v**3)) + (1.0/( log(1.0 + vcrit**(-3)) * (1.0 + vcrit**(3)))))*(1.5*spec(is)%tprim)
       end if
     
    else
       f0 = exp(-Ealpha*(v**2-1.0)/Ti) * A / (vcrit**3 + 1.0)
       df0dE = -Ealpha*f0/Ti
       gentemp = -spec(is)%temp*f0/df0dE
       f0prim = -spec(is)%fprim - (-(vcrit**3/(vcrit**3 + 1.0)) + (1.0/( log(1.0 + vcrit**(-3)) * (1.0 + vcrit**(3)))))*(ni_prim - ne_prim + Ti_prim + 0.5*Te_prim)
       f0prim = f0prim - (Ealpha/Ti)*(v**2-1.0)*Ti_prim
    end if

  end subroutine eval_f0_analytic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Equivalent Maxwellian
!! 
!! This option takes input parameters and constructs a Maxwellian alpha species with the same 
!! first two moments as the slowing-down distribution that would have been constructed by them.
!! A seperate option because the radial derivatives are treated differently than a regular
!! Maxwellian species. A possible source of confusion: input parameter "temp" is still the 
!! alpha injection energy! It gets changed here to the equivalent temperature.

  subroutine calculate_f0_arrays_equivmaxw(is)
    use constants, only: pi
    use mp, only: mp_abort, proc0, broadcast
    use species, only: spec,nspec
    implicit none
    integer, intent(in):: is
    integer:: electron_spec, i, ie
    real:: vti,ni,Zi,Ti,Ti_prim,ni_prim,vte,ne,ne_prim,Te_prim,vcva,Ealpha,vta
    real:: veff2va2, dveff2dvc, temp, LTa, Teff, mass, z

    f0_values(:,is) = exp(-egrid(:,is)) / (pi**1.5)

    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)

    do ie = 1,negrid
       generalised_temperature(:,is) = spec(is)%temp
       if (print_egrid) write(*,*) ie, egrid(ie,is), f0_values(ie,is), & 
                                   generalised_temperature(ie,is)
    end do

    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)
 
    f0prim(:,is) = -( spec(is)%fprim + (egrid(:,is) - 1.5)*spec(is)%tprim)


!    if (proc0) then
!       write(*,*) "Using equivalent-Maxwellian option for species ", is
!       write(*,*) "Parameters have been changed to the following values:"
!       write(*,*) "  vc/va = ", vcva
!       write(*,*) "  temp = ", temp
!       write(*,*) "  tprim = ", LTa
!    end if


  end subroutine calculate_f0_arrays_equivmaxw

  subroutine eval_f0_kappa(is,v,f0,gentemp,f0prim)
    use species, only: ion_species, electron_species, alpha_species, spec, nspec
    use constants, only: pi
    use mp, only: mp_abort, proc0
    implicit none
    integer,intent(in):: is
    real,intent(in):: v
    real,intent(inout):: f0,gentemp, f0prim
    real:: gammakmh,gammak,A,df0dv,df0dE
   
    gammakmh = exp(lgamma(kappa-0.5))  
    gammak = exp(lgamma(kappa))

    A = 1.0/(pi**1.5*gammakmh/gammak)

    f0 = A/(1.0 + v**2/kappa)**(kappa+1.0)

    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)
 
    df0dv = - 2.0*v*A*(kappa+1.0)/( kappa* (1.0 + v**2/kappa))
    df0dE = (0.5/v)*df0dv
    gentemp = -spec(is)%temp*f0/df0dE
 
    f0prim = 0.0
     
  end subroutine eval_f0_kappa

  subroutine calculate_f0_arrays_kappa(is)
    use constants, only: pi
    use species, only: spec, nspec
    use mp, only: mp_abort, proc0
    implicit none
    integer,intent(in):: is
    real:: A, Ti, Ealpha, vstar, v, df0dv, df0dE, B
    integer:: ie, i, electron_spec


    do ie = 1,negrid 
       v = sqrt(egrid(ie,is))
       call eval_f0_kappa(is,v,f0_values(ie,is),generalised_temperature(ie,is),f0prim(ie,is))
    end do

    gtempoz(:,is) = generalised_temperature(:,is) / spec(is)%z
    zogtemp(:,is) = spec(is)%z / generalised_temperature(:,is)

    if (.NOT. genquad_flag) weights(:,is) = weights(:,is) * f0_values(:,is)

  end subroutine calculate_f0_arrays_kappa


  function lgamma(x)

    use constants, only: pi

    real :: lgamma
    real, intent(in) :: x
    integer :: k
    real :: w, v, t, y

    real, parameter :: a(0:21) = (/ &
         &  0.00009967270908702825, -0.00019831672170162227, &
         & -0.00117085315349625822,  0.00722012810948319552, &
         & -0.00962213009367802970, -0.04219772092994235254, &
         &  0.16653861065243609743, -0.04200263501129018037, &
         & -0.65587807152061930091,  0.57721566490153514421, &
         &  0.99999999999999999764,  0.00004672097259011420, &
         & -0.00006812300803992063, -0.00132531159076610073, &
         &  0.00733521178107202770, -0.00968095666383935949, &
         & -0.04217642811873540280,  0.16653313644244428256, &
         & -0.04200165481709274859, -0.65587818792782740945, &
         &  0.57721567315209190522,  0.99999999973565236061 /)
    real, parameter :: b(0:97) = (/ &
         & -0.00000000004587497028,  0.00000000019023633960, &
         &  0.00000000086377323367,  0.00000001155136788610, &
         & -0.00000002556403058605, -0.00000015236723372486, &
         & -0.00000316805106385740,  0.00000122903704923381, &
         &  0.00002334372474572637,  0.00111544038088797696, &
         &  0.00344717051723468982,  0.03198287045148788384, &
         & -0.32705333652955399526,  0.40120442440953927615, &
         & -0.00000000005184290387, -0.00000000083355121068, &
         & -0.00000000256167239813,  0.00000001455875381397, &
         &  0.00000013512178394703,  0.00000029898826810905, &
         & -0.00000358107254612779, -0.00002445260816156224, &
         & -0.00004417127762011821,  0.00112859455189416567, &
         &  0.00804694454346728197,  0.04919775747126691372, &
         & -0.24818372840948854178,  0.11071780856646862561, &
         &  0.00000000030279161576,  0.00000000160742167357, &
         & -0.00000000405596009522, -0.00000005089259920266, &
         & -0.00000002029496209743,  0.00000135130272477793, &
         &  0.00000391430041115376, -0.00002871505678061895, &
         & -0.00023052137536922035,  0.00045534656385400747, &
         &  0.01153444585593040046,  0.07924014651650476036, &
         & -0.12152192626936502982, -0.07916438300260539592, &
         & -0.00000000050919149580, -0.00000000115274986907, &
         &  0.00000001237873512188,  0.00000002937383549209, &
         & -0.00000030621450667958, -0.00000077409414949954, &
         &  0.00000816753874325579,  0.00002412433382517375, &
         & -0.00026061217606063700, -0.00091000087658659231, &
         &  0.01068093850598380797,  0.11395654404408482305, &
         &  0.07209569059984075595, -0.10971041451764266684, &
         &  0.00000000040119897187, -0.00000000013224526679, &
         & -0.00000001002723190355,  0.00000002569249716518, &
         &  0.00000020336011868466, -0.00000118097682726060, &
         & -0.00000300660303810663,  0.00004402212897757763, &
         & -0.00001462405876235375, -0.00164873795596001280, &
         &  0.00513927520866443706,  0.13843580753590579416, &
         &  0.32730190978254056722,  0.08588339725978624973, &
         & -0.00000000015413428348,  0.00000000064905779353, &
         &  0.00000000160702811151, -0.00000002655645793815, &
         &  0.00000007619544277956,  0.00000047604380765353, &
         & -0.00000490748870866195,  0.00000821513040821212, &
         &  0.00014804944070262948, -0.00122152255762163238, &
         & -0.00087425289205498532,  0.14438703699657968310, &
         &  0.61315889733595543766,  0.55513708159976477557, &
         &  0.00000000001049740243, -0.00000000025832017855, &
         &  0.00000000139591845075, -0.00000000021177278325, &
         & -0.00000005082950464905,  0.00000037801785193343, &
         & -0.00000073982266659145, -0.00001088918441519888, &
         &  0.00012491810452478905, -0.00049171790705139895, &
         & -0.00425707089448266460,  0.13595080378472757216, &
         &  0.89518356003149514744,  1.31073912535196238583 /)
    real, parameter :: c(0:64) = (/ &
         &  0.0000000116333640008,  -0.0000000833156123568, &
         &  0.0000003832869977018,  -0.0000015814047847688, &
         &  0.0000065010672324100,  -0.0000274514060128677, &
         &  0.0001209015360925566,  -0.0005666333178228163, &
         &  0.0029294103665559733,  -0.0180340086069185819, &
         &  0.1651788780501166204,   1.1031566406452431944, &
         &  1.2009736023470742248,   0.0000000013842760642, &
         & -0.0000000069417501176,   0.0000000342976459827, &
         & -0.0000001785317236779,   0.0000009525947257118, &
         & -0.0000052483007560905,   0.0000302364659535708, &
         & -0.0001858396115473822,   0.0012634378559425382, &
         & -0.0102594702201954322,   0.1243625515195050218, &
         &  1.3888709263595291174,   2.4537365708424422209, &
         &  0.0000000001298977078,  -0.0000000008029574890, &
         &  0.0000000049454846150,  -0.0000000317563534834, &
         &  0.0000002092136698089,  -0.0000014252023958462, &
         &  0.0000101652510114008,  -0.0000774550502862323, &
         &  0.0006537746948291078,  -0.0066014912535521830, &
         &  0.0996711934948138193,   1.6110931485817511402, &
         &  3.9578139676187162939,   0.0000000000183995642, &
         & -0.0000000001353537034,   0.0000000009984676809, &
         & -0.0000000076346363974,   0.0000000599311464148, &
         & -0.0000004868554120177,   0.0000041441957716669, &
         & -0.0000377160856623282,   0.0003805693126824884, &
         & -0.0045979851178130194,   0.0831422678749791178, &
         &  1.7929113303999329439,   5.6625620598571415285, &
         &  0.0000000000034858778,  -0.0000000000297587783, &
         &  0.0000000002557677575,  -0.0000000022705728282, &
         &  0.0000000207024992450,  -0.0000001954426390917, &
         &  0.0000019343161886722,  -0.0000204790249102570, &
         &  0.0002405181940241215,  -0.0033842087561074799, &
         &  0.0713079483483518997,   1.9467574842460867884, &
         &  7.5343642367587329552 /)
    real, parameter :: d(0:6) = (/ &
         & -0.00163312359200500807,  0.00083644533703385956, &
         & -0.00059518947575728181,  0.00079365057505415415, &
         & -0.00277777777735463043,  0.08333333333333309869, &
         &  0.91893853320467274178 /)

    w = x
    if (x < 0.) w = 1. - x
    if (w < 0.5) then
       k = 0
       if (w >= 0.25) k = 11
       y = ((((((((((a(k) * w + a(k + 1)) * w + &
            & a(k + 2)) * w + a(k + 3)) * w + a(k + 4)) * w + &
            & a(k + 5)) * w + a(k + 6)) * w + a(k + 7)) * w + &
            & a(k + 8)) * w + a(k + 9)) * w + a(k + 10)) * w
       y = -log(y)
    else if (w < 3.5) then
       t = w - 4.5 / (w + 0.5)
       k = int(t + 4)
       t = t - (k - 3.5)
       k = k * 14
       y = ((((((((((((b(k) * t + b(k + 1)) * t + &
            & b(k +  2)) * t + b(k +  3)) * t + b(k +  4)) * t + &
            & b(k +  5)) * t + b(k +  6)) * t + b(k +  7)) * t + &
            & b(k +  8)) * t + b(k +  9)) * t + b(k + 10)) * t + &
            & b(k + 11)) * t + b(k + 12)) * t + b(k + 13)
    else if (w < 8.) then
       k = (int(w)) - 3
       t = w - (k + 3.5)
       k = k * 13
       y = (((((((((((c(k) * t + c(k + 1)) * t + &
            & c(k +  2)) * t + c(k +  3)) * t + c(k +  4)) * t + &
            & c(k +  5)) * t + c(k +  6)) * t + c(k +  7)) * t + &
            & c(k +  8)) * t + c(k +  9)) * t + c(k + 10)) * t + &
            & c(k + 11)) * t + c(k + 12)
    else
       v = 1. / w
       t = v * v
       y = (((((d(0) * t + d(1)) * t + d(2)) * t + &
            & d(3)) * t + d(4)) * t + d(5)) * v + d(6)
       y = y + ((w - 0.5) * log(w) - w)
    end if
    if (x < 0.) y = log(pi / sin(pi * x)) - y
    lgamma = y
  end function lgamma

end module general_f0
