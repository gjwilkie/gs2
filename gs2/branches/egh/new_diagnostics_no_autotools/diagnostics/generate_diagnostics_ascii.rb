
# This Ruby script automatically generates the module diagnostics_ascii. It is compatible with any version of Ruby, which means that it should work out of the box on any system, since even the most decrepit systems usually have Ruby 1.8.7. It is automatically invoked by the Makefile; thus, after editing it a simple make will trigger the generation of the file diagnostics_ascii.f90.
#

# This is a list of all the 
# ascii output files
ascii_files = [
  'fields',
  'heat',
  'heat2',
  'vres',
  'lpc',
  'vres2',
  'parity',
  'jext'
]


class Generator
	attr_accessor :name
	def initialize(name)
		@name = name
	end
	def declaration
		"integer :: #@name"
	end
	def switch
		"logical :: write_to_#@name = .false."
	end
  def open_output_file
    "if (ascii_files%write_to_#@name) call open_output_file(ascii_files%#@name, '.new.#@name')"
  end
  def close_output_file
    "if (ascii_files%write_to_#@name) call close_output_file(ascii_files%#@name)"
  end

end

generators = ascii_files.map{|name| Generator.new(name)}


string = <<EOF

!> DO NOT EDIT THIS FILE
!! This file has been automatically generated using generate_diagnostics_ascii.rb

!> A module for managing text-based output files 
!! for the new diagnostics module
module diagnostics_ascii


  !> Holds integers corresponding to the
  !! ascii output units
  type diagnostics_ascii_type
   #{generators.map{|g| g.declaration}.join("\n   ") }
   #{generators.map{|g| g.switch}.join("\n   ") }
  end type diagnostics_ascii_type



contains
  subroutine init_diagnostics_ascii(ascii_files)
    use file_utils, only: open_output_file
    type(diagnostics_ascii_type), intent(out) :: ascii_files
    #{generators.map{|g| g.open_output_file}.join("\n    ") }
  end subroutine init_diagnostics_ascii

  subroutine finish_diagnostics_ascii(ascii_files)
    use file_utils, only: close_output_file
    type(diagnostics_ascii_type), intent(out) :: ascii_files
    #{generators.map{|g| g.close_output_file}.join("\n    ") }
  end subroutine finish_diagnostics_ascii

end module diagnostics_ascii



EOF

File.open(ARGV[-1], 'w'){|f| f.puts string}
