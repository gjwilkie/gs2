# A file for generating overrides.f90
# To run:
#   ruby overrides.rb overrides.f90
#
# This is free software released under the MIT license
# Written by:
#           Edmund Highcock (edmundhighcock@users.sourceforge.net)
#

class Generator
  def self.generate_type(name, parameter_list)
    return <<EOF
  type #{name}_overrides_type
    logical :: set = .false.
    #{parameter_list.map{|p| p.switch}.join("\n    ")}
    #{parameter_list.map{|p| p.value}.join("\n    ")}
  end type #{name}_overrides_type

EOF
  end
  def self.generate_initialize(name, parameter_list)
    return <<EOF
  #{ 
  if countpar = parameter_list.find{|p| p.count}
   "subroutine init_#{name}_overrides(overrides, #{countpar.count})
    integer, intent(in) :: #{countpar.count}"
  else
   "subroutine init_#{name}_overrides(overrides)"
  end}
    type(#{name}_overrides_type), intent(inout) :: overrides
    overrides%set = .true.
    #{parameter_list.map{|p| p.init}.join("\n    ")}
  end subroutine init_#{name}_overrides

EOF
  end
  def self.generate_finish(name, parameter_list)
    return <<EOF
  subroutine finish_#{name}_overrides(overrides)
    type(#{name}_overrides_type), intent(inout) :: overrides
    overrides%set = .false.
    #{parameter_list.map{|p| p.finish}.join("\n    ")}
  end subroutine finish_#{name}_overrides

EOF
  end

  def switch
    if @count
      "logical, dimension(:), pointer :: override_#@name"
    else
      "logical :: override_#@name"
    end
  end

  def value
    if @count
      "#@type, dimension(:), pointer :: #@name"
    else
      "#@type :: #@name"
    end
  end

  def init
    str = "overrides%override_#@name = .false."
    str = "allocate(overrides%override_#@name(#@count), overrides%#@name(#@count))\n    " + str if @count 
    return str
  end

  def finish
    str = "overrides%override_#@name = .false."
    str = "deallocate(overrides%override_#@name, overrides%#@name)\n    " + str if @count 
    return str
  end

  attr_reader :count

  def initialize(p)
    @type = p[0]
    @name = p[1]
    @count = p[2]
  end
end 



parameter_list_geo = [
['real', 'rhoc'],
['real', 'qinp'],
['real', 'shat'],
['real', 'rgeo_lcfs'],
['real', 'rgeo_local'],
['real', 'akappa'],
['real', 'akappri'],
['real', 'tri'],
['real', 'tripri'],
['real', 'shift'],
['real', 'betaprim'],
].compact.map{|p| Generator.new(p)}

parameter_list_profs = [
['real', 'dens', 'nspec'],
['real', 'temp', 'nspec'],
['real', 'tprim', 'nspec'],
['real', 'fprim', 'nspec'],
['real', 'vnewk', 'nspec'],
['real', 'g_exb'],
['real', 'mach'],
].compact.map{|p| Generator.new(p)}


string = <<EOF
! DO NOT EDIT THIS FILE
! This file is automatically generated by overrides.rb

module overrides
#{Generator.generate_type('miller_geometry', parameter_list_geo)}
#{Generator.generate_type('profiles', parameter_list_profs)}

contains
#{Generator.generate_initialize('miller_geometry', parameter_list_geo)}
#{Generator.generate_finish('miller_geometry', parameter_list_geo)}
#{Generator.generate_initialize('profiles', parameter_list_profs)}
#{Generator.generate_finish('profiles', parameter_list_profs)}

end module overrides
EOF

File.open(ARGV[-1], 'w'){|f| f.puts string}


