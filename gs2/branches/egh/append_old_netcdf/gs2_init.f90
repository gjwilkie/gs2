! DO NOT EDIT THIS FILE
! This file is automatically generated by 
! gs2_init





!> This module is analogous to the init() function
!! in Linux-based operating systems: it initialises
!! gs2 to a certain init_level. At a given init level,
!! certain modules are initialised and certain are not.
!!
!! The gs2_init module is used by gs2_main to initialise modules. A typical
!! additional use case for this module is when it is desired
!! to override a given parameter (as in the override_* functions
!! in gs2_main). Gs2 must be taken down to the appropriate 
!! init_level, where all modules which contain any of those
!! parameters are uninitialized. The override is then set
!! and gs2 is brought back up to the highest init_level.
!! 
!! As in Linux, this module cannot be used until a certain
!! basic initialization has happened (think loading the kernel).
!! This basic initialization occurs in gs2_initialize in gs2_main,
!! and set the init_level to gs2_initialized.
!!
!! This is free software released under the MIT licence.
!! Written by:
!!            Edmund Highcock (edmundhighcock@users.sourceforge.net)
module gs2_init
  use overrides, only: miller_geometry_overrides_type
  use overrides, only: profiles_overrides_type
  use overrides, only: initial_values_overrides_type
  use overrides, only: optimisations_overrides_type
  public :: init_type
  !> A list of possible intialization levels.
  public :: init_level_list
  !> Bring gs2 to the target initialization
  !! level.
  public :: init
  !> Reads the gs2_init namelist
  public :: init_gs2_init
  !> Finalize the module
  public :: finish_gs2_init
  !> Save the current state of the fields and 
  !! distribution function, either to a file 
  !! or to a temporary array, depending on the
  !! value of in_memory. (NB if there is 
  !! sufficient memory to allocate a temporary
  !! array, in_memory will be overriden to 
  !! false).
  !public :: save_fields_and_dist_fn

  !> If initval_ov%in_memory is set to true (the default)
  !! write the current values of the dist fn into the 
  !! initval_ov object. If initval_ov%force_maxwell_reinit is
  !! set to false, also write the current field values into
  !! the object. 
  !!
  !! If initval_ov%in_memory is false, this will cause the current
  !! field and dist fn values to be written to file.
  !!
  !! You must have initialized the initval_ov
  !! object using init_initial_values_overrides before calling
  !! this function.
  !public :: set_initial_values_overrides_to_current_values

  public :: in_memory


  !> A type for labelling the different init
  !! levels available in gs2.
  type init_level_list_type
    !> The init_level reaches basic
    !! when initialize_gs2 has been called.
    integer :: basic = 1
    integer :: gs2_layouts = 2
    integer :: normalisations = 3
    integer :: theta_grid_params = 4
    integer :: gs2_save = 5
    integer :: init_g = 6
    integer :: species = 7
    integer :: override_optimisations = 8
    integer :: override_miller_geometry = 9
    integer :: theta_grid = 10
    integer :: kt_grids = 11
    integer :: le_grids = 12
    integer :: run_parameters = 13
    integer :: hyper = 14
    integer :: dist_fn_parameters = 15
    integer :: dist_fn_layouts = 16
    integer :: nonlinear_terms = 17
    integer :: dist_fn_arrays = 18
    integer :: dist_fn_level_1 = 19
    integer :: override_profiles = 20
    integer :: antenna = 21
    integer :: dist_fn_level_2 = 22
    integer :: override_timestep = 23
    integer :: dist_fn_level_3 = 24
    integer :: collisions = 25
    integer :: fields = 26
    integer :: override_initial_values = 27
    integer :: set_initial_values = 28
    integer :: full = 29
  end type init_level_list_type

  type(init_level_list_type) :: init_level_list

  !> A type for storing the current initialization
  !! status, as well as all the overrides.
  type init_type
    !> The current init level
    integer :: level = 0
    !> Whether or not diagnostics have been initialized
    logical :: diagnostics_initialized = .false.
    !> An object for overriding all or selected
    !! Miller geometry parameters. You must call
    !! gs2_main::prepare_miller_geometry_overrides 
    !! before setting these overrides. See 
    !! documentation for the overrides::miller_geometry_overrides_type
    !! for more information.
    type(miller_geometry_overrides_type) :: mgeo_ov
    !> An object for overriding all or selected
    !! profile parameters such as species temperature, density, and gradients
    !! as well as the flow gradient and mach number. You must call
    !! gs2_main::prepare_profiles_overrides 
    !! before setting these overrides. See 
    !! documentation for the overrides::profiles_overrides_type
    !! for more information.
    type(profiles_overrides_type) :: prof_ov
    !> An object for overriding the initial values of 
    !! the fields and distribution function. You must call
    !! gs2_main::prepare_initial_values_overrides 
    !! before setting these overrides. This override
    !! is very complicated. See 
    !! documentation for the overrides::initial_values_overrides_type
    !! for more information.
    type(initial_values_overrides_type) :: initval_ov

    !> An object for overriding non physics parameters which
    !! may alter run time and efficiency. You must call
    !! gs2_main::prepare_optimisations_overrides 
    !! before setting these overrides. 
    type(optimisations_overrides_type) :: opt_ov


    !> A list of possible init levels
    !type(init_level_list_type) :: levels
  end type init_type

  private

  !complex, dimension(:,:,:), allocatable :: phi_tmp, apar_tmp, bpar_tmp
  !logical :: fields_and_dist_fn_saved = .false.
  logical :: in_memory = .false.  
contains
  !> Initialize gs2 to the level of target_level.
  !! The init_type current contains info
  !! about the current initialization level. At the end
  !! of the subroutine, current%level is set to target_level
  subroutine init(current, target_level)
    use fields, only: init_fields
    implicit none
    type(init_type), intent(inout) :: current
    integer, intent(in) :: target_level
    integer :: i
    !logical :: up, down
    if (current%level .lt. init_level_list%basic) then
      write (*,*) "gs2_init::init cannot be called before &
       & initialize_gs2 in gs2 main"
      stop 1
    end if 


    if (current%level .eq. target_level) then
      return
    else
      if (up()) then 
        if (up() .and. current%level .lt. init_level_list%gs2_layouts) call gs2_layouts_subroutine
        if (up() .and. current%level .lt. init_level_list%normalisations) call normalisations_subroutine
        if (up() .and. current%level .lt. init_level_list%theta_grid_params) call theta_grid_params_subroutine
        if (up() .and. current%level .lt. init_level_list%gs2_save) call gs2_save_subroutine
        if (up() .and. current%level .lt. init_level_list%init_g) call init_g_subroutine
        if (up() .and. current%level .lt. init_level_list%species) call species_subroutine
        if (up() .and. current%level .lt. init_level_list%override_optimisations) call override_optimisations_subroutine
        if (up() .and. current%level .lt. init_level_list%override_miller_geometry) call override_miller_geometry_subroutine
        if (up() .and. current%level .lt. init_level_list%theta_grid) call theta_grid_subroutine
        if (up() .and. current%level .lt. init_level_list%kt_grids) call kt_grids_subroutine
        if (up() .and. current%level .lt. init_level_list%le_grids) call le_grids_subroutine
        if (up() .and. current%level .lt. init_level_list%run_parameters) call run_parameters_subroutine
        if (up() .and. current%level .lt. init_level_list%hyper) call hyper_subroutine
        if (up() .and. current%level .lt. init_level_list%dist_fn_parameters) call dist_fn_parameters_subroutine
        if (up() .and. current%level .lt. init_level_list%dist_fn_layouts) call dist_fn_layouts_subroutine
        if (up() .and. current%level .lt. init_level_list%nonlinear_terms) call nonlinear_terms_subroutine
        if (up() .and. current%level .lt. init_level_list%dist_fn_arrays) call dist_fn_arrays_subroutine
        if (up() .and. current%level .lt. init_level_list%dist_fn_level_1) call dist_fn_level_1_subroutine
        if (up() .and. current%level .lt. init_level_list%override_profiles) call override_profiles_subroutine
        if (up() .and. current%level .lt. init_level_list%antenna) call antenna_subroutine
        if (up() .and. current%level .lt. init_level_list%dist_fn_level_2) call dist_fn_level_2_subroutine
        if (up() .and. current%level .lt. init_level_list%override_timestep) call override_timestep_subroutine
        if (up() .and. current%level .lt. init_level_list%dist_fn_level_3) call dist_fn_level_3_subroutine
        if (up() .and. current%level .lt. init_level_list%collisions) call collisions_subroutine
        if (up() .and. current%level .lt. init_level_list%fields) call fields_subroutine
        if (up() .and. current%level .lt. init_level_list%override_initial_values) call override_initial_values_subroutine
        if (up() .and. current%level .lt. init_level_list%set_initial_values) call set_initial_values_subroutine
        if (up() .and. current%level .lt. init_level_list%full) call full_subroutine
      else if (down()) then
        if (down () .and. current%level .le. init_level_list%full) call full_subroutine
        if (down () .and. current%level .le. init_level_list%set_initial_values) call set_initial_values_subroutine
        if (down () .and. current%level .le. init_level_list%override_initial_values) call override_initial_values_subroutine
        if (down () .and. current%level .le. init_level_list%fields) call fields_subroutine
        if (down () .and. current%level .le. init_level_list%collisions) call collisions_subroutine
        if (down () .and. current%level .le. init_level_list%dist_fn_level_3) call dist_fn_level_3_subroutine
        if (down () .and. current%level .le. init_level_list%override_timestep) call override_timestep_subroutine
        if (down () .and. current%level .le. init_level_list%dist_fn_level_2) call dist_fn_level_2_subroutine
        if (down () .and. current%level .le. init_level_list%antenna) call antenna_subroutine
        if (down () .and. current%level .le. init_level_list%override_profiles) call override_profiles_subroutine
        if (down () .and. current%level .le. init_level_list%dist_fn_level_1) call dist_fn_level_1_subroutine
        if (down () .and. current%level .le. init_level_list%dist_fn_arrays) call dist_fn_arrays_subroutine
        if (down () .and. current%level .le. init_level_list%nonlinear_terms) call nonlinear_terms_subroutine
        if (down () .and. current%level .le. init_level_list%dist_fn_layouts) call dist_fn_layouts_subroutine
        if (down () .and. current%level .le. init_level_list%dist_fn_parameters) call dist_fn_parameters_subroutine
        if (down () .and. current%level .le. init_level_list%hyper) call hyper_subroutine
        if (down () .and. current%level .le. init_level_list%run_parameters) call run_parameters_subroutine
        if (down () .and. current%level .le. init_level_list%le_grids) call le_grids_subroutine
        if (down () .and. current%level .le. init_level_list%kt_grids) call kt_grids_subroutine
        if (down () .and. current%level .le. init_level_list%theta_grid) call theta_grid_subroutine
        if (down () .and. current%level .le. init_level_list%override_miller_geometry) call override_miller_geometry_subroutine
        if (down () .and. current%level .le. init_level_list%override_optimisations) call override_optimisations_subroutine
        if (down () .and. current%level .le. init_level_list%species) call species_subroutine
        if (down () .and. current%level .le. init_level_list%init_g) call init_g_subroutine
        if (down () .and. current%level .le. init_level_list%gs2_save) call gs2_save_subroutine
        if (down () .and. current%level .le. init_level_list%theta_grid_params) call theta_grid_params_subroutine
        if (down () .and. current%level .le. init_level_list%normalisations) call normalisations_subroutine
        if (down () .and. current%level .le. init_level_list%gs2_layouts) call gs2_layouts_subroutine
      end if
    end if
  contains
    function up()
      logical :: up
      up = (target_level>current%level) 
    end function up
    function down()
      logical :: down
      down = (target_level<current%level) 
    end function down
      subroutine gs2_layouts_subroutine
        use unit_tests, only: debug_message
        use gs2_layouts, only: init_gs2_layouts
        use gs2_layouts, only: finish_gs2_layouts
        if (up()) call init_gs2_layouts
        if (down()) call finish_gs2_layouts
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... gs2_layouts   ')
          current%level = 2
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... gs2_layouts   ')
          current%level = 2 - 1
        end if
      end subroutine gs2_layouts_subroutine

      subroutine normalisations_subroutine
        use unit_tests, only: debug_message
        use normalisations, only: init_normalisations
        use normalisations, only: finish_normalisations
        if (up()) call init_normalisations
        if (down()) call finish_normalisations
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... normalisations   ')
          current%level = 3
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... normalisations   ')
          current%level = 3 - 1
        end if
      end subroutine normalisations_subroutine

      subroutine theta_grid_params_subroutine
        use unit_tests, only: debug_message
        use theta_grid_params, only: init_theta_grid_params
        use theta_grid_params, only: finish_theta_grid_params
        if (up()) call init_theta_grid_params
        if (down()) call finish_theta_grid_params
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... theta_grid_params   ')
          current%level = 4
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... theta_grid_params   ')
          current%level = 4 - 1
        end if
      end subroutine theta_grid_params_subroutine

      subroutine gs2_save_subroutine
        use unit_tests, only: debug_message
        use gs2_save, only: init_gs2_save
        use gs2_save, only: finish_gs2_save
        if (up()) call init_gs2_save
        if (down()) call finish_gs2_save
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... gs2_save   ')
          current%level = 5
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... gs2_save   ')
          current%level = 5 - 1
        end if
      end subroutine gs2_save_subroutine

      subroutine init_g_subroutine
        use unit_tests, only: debug_message
        use init_g, only: init_init_g
        use init_g, only: finish_init_g
        if (up()) call init_init_g
        if (down()) call finish_init_g
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... init_g   ')
          current%level = 6
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... init_g   ')
          current%level = 6 - 1
        end if
      end subroutine init_g_subroutine

      subroutine species_subroutine
        use unit_tests, only: debug_message
        use species, only: init_species
        use species, only: finish_species
        if (up()) call init_species
        if (down()) call finish_species
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... species   ')
          current%level = 7
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... species   ')
          current%level = 7 - 1
        end if
      end subroutine species_subroutine

      subroutine override_optimisations_subroutine
        use unit_tests, only: debug_message
          use gs2_layouts, only: lso=>set_overrides
          if (up() .and. current%opt_ov%init) call lso(current%opt_ov)

       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... override_optimisations   ')
          current%level = 8
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... override_optimisations   ')
          current%level = 8 - 1
        end if
      end subroutine override_optimisations_subroutine

      subroutine override_miller_geometry_subroutine
        use unit_tests, only: debug_message
          use theta_grid_params, only: tgpso=>set_overrides
          if (up() .and. current%mgeo_ov%init) call tgpso(current%mgeo_ov)

       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... override_miller_geometry   ')
          current%level = 9
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... override_miller_geometry   ')
          current%level = 9 - 1
        end if
      end subroutine override_miller_geometry_subroutine

      subroutine theta_grid_subroutine
        use unit_tests, only: debug_message
        use theta_grid, only: init_theta_grid
        use theta_grid, only: finish_theta_grid
        if (up()) call init_theta_grid
        if (down()) call finish_theta_grid
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... theta_grid   ')
          current%level = 10
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... theta_grid   ')
          current%level = 10 - 1
        end if
      end subroutine theta_grid_subroutine

      subroutine kt_grids_subroutine
        use unit_tests, only: debug_message
        use kt_grids, only: init_kt_grids
        use kt_grids, only: finish_kt_grids
        if (up()) call init_kt_grids
        if (down()) call finish_kt_grids
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... kt_grids   ')
          current%level = 11
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... kt_grids   ')
          current%level = 11 - 1
        end if
      end subroutine kt_grids_subroutine

      subroutine le_grids_subroutine
        use unit_tests, only: debug_message
        use le_grids, only: init_le_grids
        use le_grids, only: finish_le_grids
        logical :: dummy1, dummy2
        if (up()) call init_le_grids(dummy1, dummy2)
        if (down()) call finish_le_grids
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... le_grids   ')
          current%level = 12
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... le_grids   ')
          current%level = 12 - 1
        end if
      end subroutine le_grids_subroutine

      subroutine run_parameters_subroutine
        use unit_tests, only: debug_message
        use run_parameters, only: init_run_parameters
        use run_parameters, only: finish_run_parameters
        if (up()) call init_run_parameters
        if (down()) call finish_run_parameters
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... run_parameters   ')
          current%level = 13
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... run_parameters   ')
          current%level = 13 - 1
        end if
      end subroutine run_parameters_subroutine

      subroutine hyper_subroutine
        use unit_tests, only: debug_message
        use hyper, only: init_hyper
        use hyper, only: finish_hyper
        if (up()) call init_hyper
        if (down()) call finish_hyper
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... hyper   ')
          current%level = 14
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... hyper   ')
          current%level = 14 - 1
        end if
      end subroutine hyper_subroutine

      subroutine dist_fn_parameters_subroutine
        use unit_tests, only: debug_message
        use dist_fn, only: init_dist_fn_parameters
        use dist_fn, only: finish_dist_fn_parameters
        if (up()) call init_dist_fn_parameters
        if (down()) call finish_dist_fn_parameters
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... dist_fn_parameters   ')
          current%level = 15
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... dist_fn_parameters   ')
          current%level = 15 - 1
        end if
      end subroutine dist_fn_parameters_subroutine

      subroutine dist_fn_layouts_subroutine
        use unit_tests, only: debug_message
        use gs2_layouts, only: init_dist_fn_layouts
        use gs2_layouts, only: finish_dist_fn_layouts
        use kt_grids, only: naky, ntheta0
        use le_grids, only: nlambda, negrid
        use species, only: nspec
        if (up()) call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec) 
        if (down()) call finish_dist_fn_layouts
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... dist_fn_layouts   ')
          current%level = 16
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... dist_fn_layouts   ')
          current%level = 16 - 1
        end if
      end subroutine dist_fn_layouts_subroutine

      subroutine nonlinear_terms_subroutine
        use unit_tests, only: debug_message
        use nonlinear_terms, only: init_nonlinear_terms
        use nonlinear_terms, only: finish_nonlinear_terms
        if (up()) call init_nonlinear_terms
        if (down()) call finish_nonlinear_terms
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... nonlinear_terms   ')
          current%level = 17
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... nonlinear_terms   ')
          current%level = 17 - 1
        end if
      end subroutine nonlinear_terms_subroutine

      subroutine dist_fn_arrays_subroutine
        use unit_tests, only: debug_message
        use dist_fn, only: init_dist_fn_arrays
        use dist_fn, only: finish_dist_fn_arrays
        if (up()) call init_dist_fn_arrays
        if (down()) call finish_dist_fn_arrays
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... dist_fn_arrays   ')
          current%level = 18
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... dist_fn_arrays   ')
          current%level = 18 - 1
        end if
      end subroutine dist_fn_arrays_subroutine

      subroutine dist_fn_level_1_subroutine
        use unit_tests, only: debug_message
        use dist_fn, only: init_dist_fn_level_1
        use dist_fn, only: finish_dist_fn_level_1
        if (up()) call init_dist_fn_level_1
        if (down()) call finish_dist_fn_level_1
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... dist_fn_level_1   ')
          current%level = 19
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... dist_fn_level_1   ')
          current%level = 19 - 1
        end if
      end subroutine dist_fn_level_1_subroutine

      subroutine override_profiles_subroutine
        use unit_tests, only: debug_message
          use dist_fn, only: dfso=>set_overrides
          use species, only: sso=>set_overrides
          if (up() .and. current%prof_ov%init) call dfso(current%prof_ov)
          if (up() .and. current%prof_ov%init) call sso(current%prof_ov)

       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... override_profiles   ')
          current%level = 20
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... override_profiles   ')
          current%level = 20 - 1
        end if
      end subroutine override_profiles_subroutine

      subroutine antenna_subroutine
        use unit_tests, only: debug_message
        use antenna, only: init_antenna
        use antenna, only: finish_antenna
        if (up()) call init_antenna
        if (down()) call finish_antenna
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... antenna   ')
          current%level = 21
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... antenna   ')
          current%level = 21 - 1
        end if
      end subroutine antenna_subroutine

      subroutine dist_fn_level_2_subroutine
        use unit_tests, only: debug_message
        use dist_fn, only: init_dist_fn_level_2
        use dist_fn, only: finish_dist_fn_level_2
        if (up()) call init_dist_fn_level_2
        if (down()) call finish_dist_fn_level_2
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... dist_fn_level_2   ')
          current%level = 22
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... dist_fn_level_2   ')
          current%level = 22 - 1
        end if
      end subroutine dist_fn_level_2_subroutine

      subroutine override_timestep_subroutine
        use unit_tests, only: debug_message

       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... override_timestep   ')
          current%level = 23
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... override_timestep   ')
          current%level = 23 - 1
        end if
      end subroutine override_timestep_subroutine

      subroutine dist_fn_level_3_subroutine
        use unit_tests, only: debug_message
        use dist_fn, only: init_dist_fn_level_3
        use dist_fn, only: finish_dist_fn_level_3
        if (up()) call init_dist_fn_level_3
        if (down()) call finish_dist_fn_level_3
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... dist_fn_level_3   ')
          current%level = 24
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... dist_fn_level_3   ')
          current%level = 24 - 1
        end if
      end subroutine dist_fn_level_3_subroutine

      subroutine collisions_subroutine
        use unit_tests, only: debug_message
        use collisions, only: init_collisions
        use collisions, only: finish_collisions
        if (up()) call init_collisions
        if (down()) call finish_collisions
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... collisions   ')
          current%level = 25
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... collisions   ')
          current%level = 25 - 1
        end if
      end subroutine collisions_subroutine

      subroutine fields_subroutine
        use unit_tests, only: debug_message
        use fields, only: init_fields, fields_pre_init
        use fields, only: finish_fields
        if (up()) then 
          call fields_pre_init
          !write (*,*) 'called fields_pre_init'
          call init_fields
        end if
        if (down()) call finish_fields
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... fields   ')
          current%level = 26
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... fields   ')
          current%level = 26 - 1
        end if
      end subroutine fields_subroutine

      subroutine override_initial_values_subroutine
        use unit_tests, only: debug_message

       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... override_initial_values   ')
          current%level = 27
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... override_initial_values   ')
          current%level = 27 - 1
        end if
      end subroutine override_initial_values_subroutine

      subroutine set_initial_values_subroutine
        use unit_tests, only: debug_message
          if (up()) call set_initial_field_and_dist_fn_values(current)
       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... set_initial_values   ')
          current%level = 28
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... set_initial_values   ')
          current%level = 28 - 1
        end if
      end subroutine set_initial_values_subroutine

      subroutine full_subroutine
        use unit_tests, only: debug_message

       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... full   ')
          current%level = 29
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... full   ')
          current%level = 29 - 1
        end if
      end subroutine full_subroutine



  end subroutine init

  subroutine init_gs2_init(current)
    use file_utils, only: input_unit, input_unit_exist
    use mp, only: proc0
    namelist /init_knobs/ in_memory
    type(init_type), intent(inout) :: current
    integer :: in_file
    logical :: exist

    in_memory = .false.
    if (proc0) then
      in_file = input_unit_exist("init_knobs",exist)
      if(exist) read (unit=in_file, nml=init_knobs)
    endif
  end subroutine init_gs2_init

  subroutine finish_gs2_init(current)
    type(init_type), intent(inout) :: current
    !if (fields_and_dist_fn_saved .and. in_memory) then
    !call deallocate_saved_arrays
    !end if
  end subroutine finish_gs2_init

  !subroutine save_fields_and_dist_fn
    !use dist_fn_arrays, only: gnew, g_restart_tmp
    !use gs2_save, only: gs2_save_for_restart
    !use mp, only: proc0, broadcast, mp_abort
    !use collisions, only: vnmult
    !use gs2_layouts, only: g_lo
    !use theta_grid, only: ntgrid
    !use file_utils, only: error_unit
    !use kt_grids, only: ntheta0, naky
    !use run_parameters, only: fphi, fapar, fbpar
    !use fields, only:  force_maxwell_reinit
    !use fields_arrays, only: phinew, aparnew, bparnew
    !use gs2_time, only: user_time, user_dt
    !use antenna, only: dump_ant_amp
    !implicit none
    !integer :: iostat
    !integer :: istatus

    !!write (*,*) 'save_fields_and_dist_fn in_memory', in_memory

    !if (fields_and_dist_fn_saved) then 
      !call mp_abort("ERROR: In save_fields_and_dist_fn &
        !& fields_and_dist_fn_saved is .true. This means that &
        !& save_fields_and_dist_fn has been called twice without &
        !& reinitalising.", .true.)
    !end if



    !!If we want to do restarts in memory then try to allocate storage
    !if(in_memory)then
      !!Try to allocate storage to hold g
      !allocate(g_restart_tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc),stat=iostat)

      !!If allocate failed
      !if(iostat.ne.0)then
        !!Disable in_memory flag
        !in_memory=.false.
        !!Print error message
        !if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for g --> Reverting to file based restart"
      !else
        !!Copy into temporary
        !g_restart_tmp=gnew
      !endif

      !!!!!
      !!! NOW WE MAKE COPIES OF THE FIELDS
      !!! --> Don't bother if force_maxwell_reinit as we're going to recalculate
      !!!!!

      !!Try to allocate storage to hold phi
      !if(fphi.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
        !allocate(phi_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

        !!If allocate failed
        !if(iostat.ne.0)then
          !!Disable in_memory flag
          !in_memory=.false.
          !!Print error message
          !if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for phi --> Reverting to file based restart"
        !else
          !!Copy into temporary
          !phi_tmp=phinew
        !endif
      !endif

      !!Try to allocate storage to hold apar
      !if(fapar.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
        !allocate(apar_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

        !!If allocate failed
        !if(iostat.ne.0)then
          !!Disable in_memory flag
          !in_memory=.false.
          !!Print error message
          !if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for apar --> Reverting to file based restart"
        !else
          !!Copy into temporary
          !apar_tmp=aparnew
        !endif
      !endif

      !!Try to allocate storage to hold bpar
      !if(fbpar.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
        !allocate(bpar_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

        !!If allocate failed
        !if(iostat.ne.0)then
          !!Disable in_memory flag
          !in_memory=.false.
          !!Print error message
          !if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for bpar --> Reverting to file based restart"
        !else
          !!Copy into temporary
          !bpar_tmp=bparnew
        !endif
      !endif

    !endif

    !if(.not.in_memory)then
      !!Should really do this with in_memory=.true. as well but
      !!not sure that we really need to as we never read in the dumped data.
      !if (proc0) call dump_ant_amp

      !call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    !endif
    !fields_and_dist_fn_saved = .true.

  !end subroutine save_fields_and_dist_fn
  subroutine set_initial_field_and_dist_fn_values(current)
    use dist_fn_arrays, only: g, gnew
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields, only: force_maxwell_reinit
    use fields, only: set_init_fields
    use file_utils, only: error_unit
    use init_g, only: ginit
    use init_g, only: ginitopt_restart_memory
    use init_g, only: ginitopt_restart_many
    use run_parameters, only: fphi, fapar, fbpar
    use mp, only: proc0
    use unit_tests, only: job_id
    implicit none
    type (init_type), intent(in) :: current
    logical :: restarted

    !write (*,*) 'set_init_field in_memory', in_memory, fields_and_dist_fn_saved

    if (.not. current%initval_ov%override) then 
      ! This is the usual initial setup 
      call ginit (restarted)
    else
      if (current%initval_ov%in_memory) then 
        g = current%initval_ov%g 
        gnew = g
      else
        call ginit(restarted, ginitopt_restart_many)
      end if
    end if

    if (current%initval_ov%force_maxwell_reinit)then
      call set_init_fields
      !if (.not. proc0) write (*,*) 'field value jjjj kkkk', phinew(-5, 3, 2), job_id
      !if (.not. proc0) write (*,*) 'field value sum jjjj kkkk', sum(real(conjg(phinew)*phinew)), job_id
    else
      write(error_unit(), *) "INFO: You have disabled &
        & force_maxwell_reinit which causes the fields to be &
        & recalculated self-consistently from the the dist fn. &
        & You are on your own and here be dragons."
      if(current%initval_ov%override .and. current%initval_ov%in_memory) then 
        if(fphi.gt.0) &
          phinew=current%initval_ov%phi
        if(fapar.gt.0) &
          aparnew=current%initval_ov%apar
        if(fbpar.gt.0) &
          bparnew=current%initval_ov%bpar
      else
        ! No need to do anything: fields read from file.
      end if
    end if
  end subroutine set_initial_field_and_dist_fn_values
  !subroutine deallocate_saved_arrays 
    !use dist_fn_arrays, only: g_restart_tmp
    !if(allocated(g_restart_tmp)) deallocate(g_restart_tmp)
    !if(allocated(phi_tmp)) deallocate(phi_tmp)
    !if(allocated(apar_tmp)) deallocate(apar_tmp)
    !if(allocated(bpar_tmp)) deallocate(bpar_tmp)
  !end subroutine deallocate_saved_arrays
end module gs2_init

