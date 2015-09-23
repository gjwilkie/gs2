# This file generates the file gs2_init.f90.
# For more information see the help in gs2_init.f90
# If you want to add a new module to the initialization
# system, see the help for the DEPENDENCIES constant
# just below. 
#
# To run:
#    $ ruby gs2_init.rb gs2_init.f90
#
# This is free software released under the MIT licence.
# Written by:
#           Edmund Highcock (edmundhighcock@users.sourceforge.net)

class GenerateInit

  # If you want to add a new module to this list, what you
  # have to do is include it in the following form:
  #
  # ['new_module_name', ['dependency1', 'dependency2']],
  #
  # where dependency1 etc are modules that need to be 
  # initialized before your new module.
  #
  # Importantly, you also have to make sure that your
  # module has the public subroutines init_new_module_name
  # and finish_new_module_name. This module will call
  # the first subroutine when initializing your module,
  # and the second when finishing it.
  #
  # More sophisticated initialisation of your module is
  # possible, (e.g. passing parameters to the init
  # or finish functions, but then you have to add your own case 
  # to the subroutine function below. See dist_fn_layouts
  # as an example.
  #
  # Initialization levels should be listed in order of cost
  # with highest cost at the top. Overrides should all be at
  # the bottome to make sure they appear as high as possible
  # in the level list.
  DEPENDENCIES = [ 
    # full is a special level to signify complete
    # initialiation (excluding diagnostics, which are not handled
    # by this module)
    ['full', ['set_initial_values', 'normalisations']],
    ['fields_level_2' , ['collisions', 'antenna_level_2', 'dist_fn_level_3', 'fields_level_1']],
    ['dist_fn_level_3' , ['dist_fn_level_2', 'hyper', 'override_timestep']],
    ['dist_fn_level_1' , ['dist_fn_arrays']], 
    ['fields_level_1' , ['kt_grids', 'antenna_level_1', 'gs2_layouts', 'fields_parameters']],
    ['collisions' ,
      ['species', 'kt_grids', 'gs2_layouts', 'theta_grid', 'le_grids',
        'dist_fn_layouts', 'run_parameters', 'dist_fn_level_3']],
    ['nonlinear_terms' ,
      ['species', 'kt_grids', 'gs2_layouts', 'theta_grid', 'le_grids',
        'dist_fn_layouts', 'override_optimisations']],
    ['gs2_layouts' , []],
    ['set_initial_values', ['fields_level_2', 'init_g', 'override_initial_values']],
    ['le_grids' ,
      ['species', 'kt_grids', 'gs2_layouts', 'theta_grid']],
    ['antenna_level_2' , ['species', 'run_parameters', 'override_profiles']],
    ['theta_grid' , ['theta_grid_params', 'override_miller_geometry']],
    ['normalisations', []],
    ['theta_grid_params' , []],
    ['fields_parameters' , []],
    ['kt_grids_parameters' , ['theta_grid']],
    ['kt_grids' , ['theta_grid', 'kt_grids_parameters', 'override_kt_grids']],
    ['gs2_save' , []],
    ['run_parameters' , ['kt_grids']],
    ['hyper' , ['kt_grids', 'gs2_layouts']],
    ['init_g' , ['gs2_layouts']],
    ['species' , ['kt_grids']],
    ['dist_fn_parameters' , 
      ['gs2_layouts', 'species', 'theta_grid', 'kt_grids', 'le_grids'   ]],
    ['dist_fn_arrays' , 
      ['dist_fn_parameters', 'run_parameters', 'nonlinear_terms',
        'dist_fn_layouts','nonlinear_terms' ]],
    ['dist_fn_layouts' ,
      ['species', 'kt_grids', 'gs2_layouts', 'theta_grid']],

    ['dist_fn_level_2' , ['dist_fn_level_1', 'override_profiles']], 
    ['antenna_level_1', ['species']],
    
    ['override_kt_grids' , ['kt_grids_parameters']],
    ['override_optimisations' , ['gs2_layouts', 'fields_parameters']],
    ['override_miller_geometry' , ['theta_grid_params']],
    # Override tprim, fprim, vnewk, temp and dens in species
    ['override_profiles' , ['species']],
    # Override the timestep set in run_parameters
    ['override_timestep' , ['run_parameters']],
    ['override_initial_values', ['fields_level_2', 'init_g']],
  ]
                  

                  


  # A list of levels in ascending order of dependence, i.e. the leftmost
  # in the array must be initialized first.
  LEVELS = []
  # Here we make a hash of {module => [dependencies]}
  # We can't make DEPENDENCIES itself a hash because 
  # we care about the order in which the modules
  # are listed and older versions of Ruby may not 
  # preserve the order of keys in a hash
  deps = {}
  DEPENDENCIES.each{|mod,dependencies| deps[mod] = dependencies}
  modules_remaining = DEPENDENCIES.map{|mod,dependencies| mod}

  # Now we make the LEVELS list. This block determines
  # the order the levels have to be reached, i.e.
  # the order in which modules are initialized.
  while modules_remaining.size > 0
    #p modules_remaining 
    modules_remaining.each do |mod|
      if deps[mod] - LEVELS == [] # i.e. all dependencies already in LEVELS
        LEVELS.push mod
        modules_remaining.delete mod
        break
      end 
    end
  end 
  

  # This is the basic level that must be reached
  # before calling the init subroutine from this
  # module.
  GS2_LEVEL = 1

  @@level_counter = GS2_LEVEL+1
  
  # The constructor for each level object.
  def initialize(level)
    @level_name = level
		@level_sub_name = level + '_subroutine'
    @module_name = case level
                   when 'dist_fn_layouts' 
                     'gs2_layouts' 
                   when /^dist_fn_*/
                     'dist_fn'
                   when /^fields_level_[12]/
                     'fields'
                   when /^antenna_level_[12]/
                     'antenna'
                   when /^kt_grids*/
                     'kt_grids'
                   when /^fields_parameters/
                     'fields'
                   else
                     @level_name
                   end
    @level_number = @@level_counter
    @@level_counter += 1
  end

  # Declaration of the level within the 
  # level_list_type object.
  def level_declaration
    "integer :: #@level_name = #@level_number"
  end 

  # Code for determining if the module should
  # be initialized.
  def up
    "if (up() .and. current%level == init_level_list%#@level_name-1) call #@level_sub_name"
  end  
  def down
    "if (down () .and. current%level == init_level_list%#@level_name) call #@level_sub_name"
  end  
  def subroutine
    return <<EOF
      subroutine #@level_sub_name
        use unit_tests, only: debug_message
#{
        case @level_name
        when 'full', 'override_initial_values'
          str = "\n"
        when /override_kt_grids/
          str = <<EOF2
          use kt_grids, only: ktso=>set_overrides
          if (up() .and. current%kt_ov%init) call ktso(current%kt_ov)

EOF2
        when /override_optimisations/
          str = <<EOF2
          use gs2_layouts, only: lso=>set_overrides
          use fields, only: fso=>set_overrides
          use dist_fn, only: dso=>set_overrides
          if (up() .and. current%opt_ov%init) call lso(current%opt_ov)
          if (up() .and. current%opt_ov%init) call fso(current%opt_ov)
          if (up() .and. current%opt_ov%init) call dso(current%opt_ov)

EOF2
        when /override_miller_geometry/
          str = <<EOF2
          use theta_grid_params, only: tgpso=>set_overrides
          if (up() .and. current%mgeo_ov%init) call tgpso(current%mgeo_ov)

EOF2
			  when /override_timestep/
          str = <<EOF2
          use run_parameters, only: rso=>set_overrides
          if (up() .and. current%tstep_ov%init) call rso(current%tstep_ov)

EOF2
        when /override_profiles/
          str = <<EOF2
          use dist_fn, only: dfso=>set_overrides
          use species, only: sso=>set_overrides
          if (up() .and. current%prof_ov%init) call dfso(current%prof_ov)
          if (up() .and. current%prof_ov%init) call sso(current%prof_ov)

EOF2
        when 'set_initial_values'
          str = <<EOF2
          if (up()) call set_initial_field_and_dist_fn_values(current)
EOF2
        when 'dist_fn_layouts'
          str = <<EOF2
        use gs2_layouts, only: init_#@level_name
        use gs2_layouts, only: finish_#@level_name
        use kt_grids, only: naky, ntheta0
        use le_grids, only: nlambda, negrid
        use species, only: nspec
        if (up()) call init_#@level_name(naky, ntheta0, nlambda, negrid, nspec) 
        if (down()) call finish_#@level_name
EOF2
        when 'le_grids'
          str = <<EOF2
        use #@module_name, only: init_#@level_name
        use #@module_name, only: finish_#@level_name
        logical :: dummy1, dummy2
        if (up()) call init_#@level_name(dummy1, dummy2)
        if (down()) call finish_#@level_name
EOF2
        when 'fields'
          str = <<EOF2
        use #@module_name, only: init_#@level_name, fields_pre_init
        use #@module_name, only: finish_#@level_name
        if (up()) then 
          call fields_pre_init
          !write (*,*) 'called fields_pre_init'
          call init_#@level_name
        end if
        if (down()) call finish_#@level_name
EOF2
        else
          str = <<EOF2
        use #@module_name, only: init_#@level_name
        use #@module_name, only: finish_#@level_name
        if (up()) call init_#@level_name
        if (down()) call finish_#@level_name
EOF2
        end
        str
        
}       
        if (up()) then
          call debug_message(1, 'gs2_init::init reached init level... #@level_name   ')
          current%level = #@level_number
        else if (down()) then  ! (down)
          call debug_message(1, 'gs2_init::init left init level... #@level_name   ')
          current%level = #@level_number - 1
        end if
      end subroutine #@level_sub_name
EOF
  end



end 

generators = GenerateInit::LEVELS.map{|l| GenerateInit.new(l)}

dir = File.dirname(File.expand_path(__FILE__)) # i.e. the directory this file is in

string = File.read(dir + '/templates/gs2_init_template.f90').sub(
  /\A.*!\s*TEMPLATE_END_HEADER/m, ''
).sub(
  /!\s*TEMPLATE_GS2_INIT_LEVEL.*$/, 
  "integer :: basic = #{GenerateInit::GS2_LEVEL}"
).sub(
  /!\s*TEMPLATE_LEVEL_DECLARATIONS.*$/,
  "#{generators.map{|g| g.level_declaration}.join("\n    ")}"
).sub(
  /!\s*TEMPLATE_UP.*$/,
  "#{generators.map{|g| g.up}.join("\n        ")}"
).sub(
  /!\s*TEMPLATE_DOWN.*$/,
  "#{generators.reverse.map{|g| g.down}.join("\n        ")}"
).sub(
  /^\s+!\s*TEMPLATE_SUBROUTINES.*$/,
  "#{generators.map{|g| g.subroutine}.join("\n")}"
).sub(
  /^\s+!\s*TEMPLATE_CHECK_GS2_INITIALISED.*$/,
  <<EOF
    if (current%level .lt. init_level_list%basic) then
      write (*,*) "gs2_init::init cannot be called before &
       & initialize_gs2 in gs2 main"
      stop 1
    end if 
EOF
)

warning = <<EOF
! DO NOT EDIT THIS FILE
! This file is automatically generated by 
! gs2_init

EOF

string = warning + string 



File.open(ARGV[-1], 'w'){|f| f.puts string}
