! This file is used by gs2_init.rb as the template for gs2_init.f90.
! gs2_init.rb takes this template and replaces certain sections of it
! with automatically generated code. 
!
! Do not modify any lines containing the word TEMPLATE
! TEMPLATE_END_HEADER



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
    ! TEMPLATE_GS2_INIT_LEVEL
    ! TEMPLATE_LEVEL_DECLARATIONS
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

    ! TEMPLATE_CHECK_GS2_INITIALISED

    if (current%level .eq. target_level) then
      return
    else
      if (up()) then 
        ! TEMPLATE_UP
      else if (down()) then
        ! TEMPLATE_DOWN
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
    ! TEMPLATE_SUBROUTINES


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

