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
    # full is a special level to signify complete
    # initialiation
    ['full', ['set_initial_values']],
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
    ['set_initial_values', ['fields', 'init_g']],
    ['le_grids' ,
      ['species', 'kt_grids', 'gs2_layouts', 'theta_grid']],
    ['antenna' , ['species', 'run_parameters', 'override_profiles']],
    ['theta_grid' , ['theta_grid_params', 'override_miller_geometry']],
    ['theta_grid_params' , []],
    ['kt_grids' , ['theta_grid']],
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

    ['dist_fn_level_2' , ['dist_fn_level_1', 'override_profiles']], 
    
    ['override_miller_geometry' , ['theta_grid_params']],
    # Override tprim, fprim, vnewk, temp and dens in species
    ['override_profiles' , ['species']],
    # Override the timestep set in run_parameters
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

  #p deps, 'deps'
  #exit

  while modules_remaining.size > 0
    #p 'modules_remaining', modules_remaining, 'LEVELS', LEVELS
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
  
  #p LEVELS, 'LEVELS'

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
    "if (up() .and. current%level .lt. init_level_list%#@level_name) call #@level_name"
  end  
  def down
    "if (down () .and. current%level .le. init_level_list%#@level_name) call #@level_name"
  end  
  def subroutine
    return <<EOF
      subroutine #@level_name
        use unit_tests, only: debug_message
#{
        case @level_name
        when 'full', 'override_timestep'
          str = "\n"
        when /override_miller_geometry/
          str = <<EOF2
          use theta_grid_params, only: tgpso=>set_overrides
          if (up() .and. current%mgeo_ov%set) call tgpso(current%mgeo_ov)

EOF2
        when /override_profiles/
          str = <<EOF2
          use dist_fn, only: dfso=>set_overrides
          use species, only: sso=>set_overrides
          if (up() .and. current%prof_ov%set) call dfso(current%prof_ov)
          if (up() .and. current%prof_ov%set) call sso(current%prof_ov)

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
      end subroutine #@level_name
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
