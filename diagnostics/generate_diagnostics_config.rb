# This Ruby script automatically generates the module diagnostics_config. It is compatible with any version of Ruby, which means that it should work out of the box on any system, since even the most decrepit systems usually have Ruby 1.8.7. It is automatically invoked by the Makefile; thus, after editing it a simple make will trigger the generation of the file diagnostics_config.f90.
#

# This is a list of all the input parameters to the namelist
# diagnostics_config, with their default values.
# Leave third element of array empty to use the default default
#
# This list is used in generating diagnostics_config.f90
input_variables_for_diagnostics_config = [
	['integer', 'nwrite', '10'],
	['logical', 'write_any', '.true.'],

	# Parameters for writing out fields
	['logical', 'write_fields', '.true.'],
	['logical', 'write_phi_over_time'],
	['logical', 'write_apar_over_time'],
	['logical', 'write_bpar_over_time'],

	# Parameters for writing out fluxes
	['logical', 'write_fluxes', '.true.'],
	['logical', 'write_fluxes_by_mode'],

	# Parameters for writing out growth rates and frequencies
	['logical', 'write_omega', '.true.'],
	['integer', 'navg', '10'],
	['integer', 'igomega', '0'],
	['real', 'omegatinst', '1.0e6'],
	['real', 'omegatol', '-0.001'],
	['logical', 'exit_when_converged', '.true.'],
]


class Generator
	attr_accessor :name, :type
	def initialize(type, name, default)
		@type = type
		@name = name
		@default = default
	end
	def declaration
		"#@type :: #@name"
	end
	def set_default
		"#@name = #{default}"
	end
	def default
		@default || case type
								when /logical/
									'.false.'
								when /integer/
									'1'
								end
	end
	def parameters_type_value
		"gnostics%#@name"
	end
	def set_parameters_type_value
		"#{parameters_type_value} = #@name"
	end
	def broadcast
		"call broadcast (#{parameters_type_value})"
	end
end

generators = input_variables_for_diagnostics_config.map{|type, name, default| Generator.new(type,name,default)}


string = <<EOF

!> DO NOT EDIT THIS FILE
!! This file has been automatically generated using generate_diagnostics_config.rb

!> A module for handling the configuration of the diagnostics
!! module via the namelist diagnostics_config.
module diagnostics_config

  use simpledataio, only: sdatio_file

  public :: init_diagnostics_config
  public :: finish_diagnostics_config
  public ::diagnostics_type

  type diagnostics_type
   type(sdatio_file) :: sfile
   !> Integer below gives the sdatio type 
   !! which corresponds to a gs2 real
   integer :: rtype
   integer :: istep
   logical :: create
   logical :: wryte
   logical :: distributed
   logical :: parallel
   logical :: exit
   #{generators.map{|g| g.declaration}.join("\n   ") }
  end type diagnostics_type


  private

contains
  subroutine init_diagnostics_config(gnostics)
    use file_utils, only: open_output_file
    type(diagnostics_type), intent(out) :: gnostics
    call read_parameters(gnostics)
  end subroutine init_diagnostics_config

  subroutine finish_diagnostics_config
  end subroutine finish_diagnostics_config

  subroutine read_parameters(gnostics)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    implicit none
    type(diagnostics_type), intent(out) :: gnostics
    #{generators.map{|g| g.declaration}.join("\n    ") }
    namelist /diagnostics_config/ &
      #{generators.map{|g| g.name}.join(", &\n      ")}

    integer :: in_file
    logical :: exist

    if (proc0) then
      #{generators.map{|g| g.set_default}.join("\n      ")}

      in_file = input_unit_exist ("diagnostics_config", exist)
      if (exist) read (unit=in_file, nml=diagnostics_config)

      #{generators.map{|g| g.set_parameters_type_value}.join("\n      ")}

    end if

    #{generators.map{|g| g.broadcast}.join("\n    ")}



  end subroutine read_parameters
end module diagnostics_config



EOF

File.open(ARGV[-1], 'w'){|f| f.puts string}