# This file generates the file gs2_init.f90.
# For more information see the help in gs2_init.f90
# To run:
#    $ ruby gs2_init.rb gs2_init.f90
#
# This is free software released under the MIT licence.
# Written by:
#           Edmund Highcock (edmundhighcock@users.sourceforge.net)

class GenerateInit

  # Initialization levels should be listed in order of cost
  # with highest cost at the top. Overrides should all be at
  # the bottome to make sure they appear as high as possible
  # in the level list.
  DEPENDENCIES = [
                  ['fields' , ['collisions', 'antenna', 'dist_fn_level_3']],
                  ['dist_fn_level_3' , ['dist_fn_level_2', 'hyper', 'override_timestep']],
                  ['dist_fn_level_1' , ['dist_fn_arrays']], 
                  ['collisions' ,
                    ['species', 'kt_grids', 'gs2_layouts', 'theta_grid', 'le_grids',
                      'dist_fn_layouts', 'run_parameters', 'dist_fn_level_3']],
                  ['nonlinear_terms' ,
                    ['species', 'kt_grids', 'gs2_layouts', 'theta_grid', 'le_grids',
                      'dist_fn_layouts']],
                  ['gs2_layouts' , []],
                  ['le_grids' ,
                    ['species', 'kt_grids', 'gs2_layouts', 'theta_grid']],
                  ['antenna' , ['species', 'run_parameters']],
                  ['theta_grid' , []],
                  ['kt_grids' , []],
                  ['gs2_save' , []],
                  ['run_parameters' , ['kt_grids']],
                  ['hyper' , ['kt_grids', 'gs2_layouts']],
                  ['init_g' , ['gs2_layouts']],
                  ['species' , []],
                  ['dist_fn_parameters' , 
                    ['gs2_layouts', 'species', 'theta_grid', 'kt_grids', 'le_grids'   ]],
                  ['dist_fn_arrays' , 
                    ['dist_fn_parameters', 'run_parameters', 'nonlinear_terms',
                      'dist_fn_layouts','nonlinear_terms' ]],
                  ['dist_fn_layouts' ,
                    ['species', 'kt_grids', 'gs2_layouts', 'theta_grid']],
                  ['dist_fn_level_2' , ['dist_fn_level_1']], 
                  ['override_timestep' , ['run_parameters']],
  ]
                  

                  


  # A list of levels in ascending order of dependence, i.e. the leftmost
  # in the array must be initialized first.
  LEVELS = []
            #['theta_grid', 'kt_grids', 'gs2_layouts', 'gs2_save'],
            #['run_parameters', 'init_g', 'species', 'le_grids'],
            #['dist_fn_layouts', 'antenna'],
            #['nonlinear_terms', 'collisions'],
            #['fields'],
           #].flatten

  # Here we make a hash of {module => [dependencies]}
  # We can't make DEPENDENCIES itself a hash because 
  # we care about the order in which the modules
  # are listed and older versions of Ruby may not 
  # preserve the order of keys in a hash
  deps = {}
  DEPENDENCIES.each{|mod,dependencies| deps[mod] = dependencies}
  modules_remaining = DEPENDENCIES.map{|mod,dependencies| mod}

  p deps, 'deps'
  #exit

  while modules_remaining.size > 0
    p 'modules_remaining', modules_remaining, 'LEVELS', LEVELS
    #deps.keys.each do |k|
      #if deps[k] - LEVELS == []
        #LEVELS.push k
        #deps.delete k
      #end 
    #end
    modules_remaining.each do |mod|
      if deps[mod] - LEVELS == [] # i.e. all dependencies already in LEVELS
        LEVELS.push mod
        modules_remaining.delete mod
        break
      end 
    end
  end 
  
  p LEVELS, 'LEVELS'

  GS2_LEVEL = 1

  @@level_counter = GS2_LEVEL+1
  
  def initialize(level)
    @level_name = level
    @module_name = case level
                   when 'dist_fn_layouts' 
                     'gs2_layouts' 
                   when /^dist_fn_*/
                     'dist_fn'
                   else
                     @level_name
                   end
    @level_number = @@level_counter
    @@level_counter += 1
  end

  def level_declaration
    "integer :: #@level_name = #@level_number"
  end 

  def up
    "if (current%level .lt. init_level_list%#@level_name) call #@level_name"
  end  
  def down
    "if (current%level .le. init_level_list%#@level_name) call #@level_name"
  end  
  def subroutine
    return <<EOF
      subroutine #@level_name
        use unit_tests, only: debug_message
#{
        case @level_name
        when /override.*/
          # Nothing needs to be done for the overrides,
          # they are just placeholders
          str = "\n"
        when 'dist_fn_layouts'
          str = <<EOF2
        use gs2_layouts, only: init_#@level_name
        use gs2_layouts, only: finish_#@level_name
        use kt_grids, only: naky, ntheta0
        use le_grids, only: nlambda, negrid
        use species, only: nspec
        if (up) call init_#@level_name(naky, ntheta0, nlambda, negrid, nspec) 
        if (down) call finish_#@level_name
EOF2
        when 'le_grids'
          str = <<EOF2
        use #@module_name, only: init_#@level_name
        use #@module_name, only: finish_#@level_name
        logical :: dummy1, dummy2
        if (up) call init_#@level_name(dummy1, dummy2)
        if (down) call finish_#@level_name
EOF2
        else
          str = <<EOF2
        use #@module_name, only: init_#@level_name
        use #@module_name, only: finish_#@level_name
        if (up) call init_#@level_name
        if (down) call finish_#@level_name
EOF2
        end
        str
        
}       
        if (up) then
          call debug_message(2, 'gs2_init::init reached init level... #@level_name   ')
          current%level = #@level_number
        else  ! (down)
          call debug_message(2, 'gs2_init::init left init level... #@level_name   ')
          current%level = #@level_number - 1
        end if
      end subroutine #@level_name
EOF
  end



end 

generators = GenerateInit::LEVELS.map{|l| GenerateInit.new(l)}

string = <<EOF
! DO NOT EDIT THIS FILE
! This file is automatically generated by 
! gs2_init.rb


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
  public :: init_level_type
  public :: init_level_list
  public :: init
  !> Save the current state of the fields and 
  !! distribution function, either to a file 
  !! or to a temporary array, depending on the
  !! value of in_memory. (NB if there is 
  !! sufficient memory to allocate a temporary
  !! array, in_memory will be overriden to 
  !! false).
  public :: save_fields_and_dist_fn
  public :: load_saved_field_values

  public :: in_memory

  private

  !> A type for labelling the different init
  !! levels available in gs2.
  type init_level_list_type
    !> The init_level reaches gs2
    !! when initialize_gs2 has been called.
    integer :: gs2 = #{GenerateInit::GS2_LEVEL}
    #{generators.map{|g| g.level_declaration}.join("\n    ")}
  end type init_level_list_type

  type(init_level_list_type) :: init_level_list
  
  !> A type for storing the init_level of gs2.
  type init_level_type
    !> The current init level
    integer :: level = 0
    !> Whether or not diagnostics have been initialized
    logical :: diagnostics_initialized = .false.
    !> A list of possible init levels
    !type(init_level_list_type) :: levels
  end type init_level_type

  complex, dimension(:,:,:), allocatable :: phi_tmp, apar_tmp, bpar_tmp
  logical :: fields_and_dist_fn_saved = .false.
  logical :: in_memory = .false.  
contains
  !> Initialize gs2 to the level of target_level.
  !! The init_level_type current contains info
  !! about the current initialization level. At the end
  !! of the subroutine, current%level is set to target_level
  subroutine init(current, target_level)
    use fields, only: init_fields
    implicit none
    type(init_level_type), intent(inout) :: current
    integer, intent(in) :: target_level
    integer :: i
    logical :: up, down
    down = (target_level<current%level) 
    up = (target_level>current%level) 

    if (current%level .lt. init_level_list%gs2) then
      write (*,*) "gs2_init::init cannot be called before &
       & initialize_gs2 in gs2 main"
      stop 1
    end if 


    if (current%level .eq. target_level) then
      return
    else
      if (up) then 
        #{generators.map{|g| g.up}.join("\n        ")}
      else
        #{generators.reverse.map{|g| g.down}.join("\n        ")}
      end if
    end if
    contains
#{generators.map{|g| g.subroutine}.join("\n")}
  end subroutine init

  subroutine save_fields_and_dist_fn
    use dist_fn_arrays, only: gnew, g_restart_tmp
    use gs2_save, only: gs2_save_for_restart
    use mp, only: proc0, broadcast
    use collisions, only: vnmult
    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use file_utils, only: error_unit
    use kt_grids, only: ntheta0, naky
    use run_parameters, only: fphi, fapar, fbpar
    use fields, only:  force_maxwell_reinit
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_time, only: user_time, user_dt
    use antenna, only: dump_ant_amp
    use file_utils, only: input_unit, input_unit_exist
    integer :: iostat, istatus
    integer :: in_file
    logical :: exist


    namelist /init_knobs/ in_memory

    if (proc0) then
       in_file = input_unit_exist("init_knobs",exist)
       if(exist) read (unit=in_file, nml=init_knobs)
    endif

    !If we want to do restarts in memory then try to allocate storage
    if(in_memory)then
       !Try to allocate storage to hold g
       allocate(g_restart_tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc),stat=iostat)

       !If allocate failed
       if(iostat.ne.0)then
          !Disable in_memory flag
          in_memory=.false.
          !Print error message
          if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for g --> Reverting to file based restart"
       else
          !Copy into temporary
          g_restart_tmp=gnew
       endif

       !!!!
       !! NOW WE MAKE COPIES OF THE FIELDS
       !! --> Don't bother if force_maxwell_reinit as we're going to recalculate
       !!!!

       !Try to allocate storage to hold phi
       if(fphi.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
          allocate(phi_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

          !If allocate failed
          if(iostat.ne.0)then
             !Disable in_memory flag
             in_memory=.false.
             !Print error message
             if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for phi --> Reverting to file based restart"
          else
             !Copy into temporary
             phi_tmp=phinew
          endif
       endif

       !Try to allocate storage to hold apar
       if(fapar.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
          allocate(apar_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

          !If allocate failed
          if(iostat.ne.0)then
             !Disable in_memory flag
             in_memory=.false.
             !Print error message
             if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for apar --> Reverting to file based restart"
          else
             !Copy into temporary
             apar_tmp=aparnew
          endif
       endif

       !Try to allocate storage to hold bpar
       if(fbpar.gt.0.and.(.not.force_maxwell_reinit).and.in_memory)then
          allocate(bpar_tmp(-ntgrid:ntgrid,ntheta0,naky),stat=iostat)

          !If allocate failed
          if(iostat.ne.0)then
             !Disable in_memory flag
             in_memory=.false.
             !Print error message
             if (proc0) write(error_unit(), *) "Couldn't allocate temporary storage for bpar --> Reverting to file based restart"
          else
             !Copy into temporary
             bpar_tmp=bparnew
          endif
       endif

    endif

    if(.not.in_memory)then
       !Should really do this with in_memory=.true. as well but
       !not sure that we really need to as we never read in the dumped data.
       if (proc0) call dump_ant_amp

       call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    endif
    fields_and_dist_fn_saved = .true.

  end subroutine save_fields_and_dist_fn
  subroutine load_saved_field_values
    use fields_arrays, only: phinew, aparnew, bparnew
    use fields, only: force_maxwell_reinit
    use dist_fn_arrays, only: g_restart_tmp
    use run_parameters, only: fphi, fapar, fbpar
    if(in_memory.and.(.not.force_maxwell_reinit))then
       if(fphi.gt.0) phinew=phi_tmp
       if(fapar.gt.0) aparnew=apar_tmp
       if(fbpar.gt.0) bparnew=bpar_tmp
    endif
    !Deallocate tmp memory
    if(allocated(g_restart_tmp)) deallocate(g_restart_tmp)
    if(allocated(phi_tmp)) deallocate(phi_tmp)
    if(allocated(apar_tmp)) deallocate(apar_tmp)
    if(allocated(bpar_tmp)) deallocate(bpar_tmp)
    fields_and_dist_fn_saved = .false.
  end subroutine load_saved_field_values
end module gs2_init

EOF

File.open(ARGV[-1], 'w'){|f| f.puts string}
