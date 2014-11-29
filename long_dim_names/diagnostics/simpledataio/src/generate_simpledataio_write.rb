class Generator
  def initialize(type, dimsize, offset)
    @dimsize = dimsize
    @offset = offset
    @type = type
    if dimsize==0
      @dimension = ""
    else
      @dimension = ", dimension(#{dimension_spec}) "
    end
  end
  def dimension_spec
    ([":"]*@dimsize).join(",")
  end
  def procedure_name
    if @offset
      "write_variable_with_offset_#{@type.gsub(/[* ]/, '_')}_#{@dimsize}"
    else
      "write_variable_#{@type.gsub(/[* ]/, '_')}_#{@dimsize}"
    end
  end
  def interface_name
    #unless @type=~/^(real|complex)$/
    "  module procedure " + procedure_name
    #else
    #"#ifdef SINGLE_PRECISION\n  module procedure #{procedure_name}\n#endif"
    #end
  end
  def get_n2
    if @offset
      <<EOF
call number_of_unlimited_dimensions(sfile, variable_name, n)
   n2 = n+#@dimsize#{complex ? "+1" : nil}
   ndims = number_of_dimensions(sfile, variable_name)
   if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to write_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use write_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
   end if
EOF
    else
      "n2 = number_of_dimensions(sfile, variable_name)"
    end
  end
  def val_get
    if @offset
      if complex  
        "realval(1:,#{@dimsize.times.map{|i| "(offsets(#{i+2})):"}.join(",&\n                 ")})" 
      else
        "val(#{@dimsize.times.map{|i| "(offsets(#{i+1})):"}.join(",&   \n")})"
      end
    else
      if complex  
        "realval(1:,#{dimension_spec})" 
      else
        "val(#{dimension_spec})"
      end
    end
  end
  def val_get_0
    complex ?  "realval" : "(/val/)"
  end
  def complex
    @type =~ /complex/i
  end
  def real_type_from_complex
    (@type =~ /16/ ? 'double precision' : 'real')   
  end
  def realval_declaration
    real_type_from_complex  + ', dimension(1' + @dimsize.times.map{|i| ", size(val, #{i+1})"}.join('') + ') :: realval'
  end
  def realimag
    ['real', 'aimag']
  end
  def netcdf_inputs
    "call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)"
  end
  def put_variable
    "status =  nf90_put_var(fileid, varid+1, &
        #{@dimsize == 0 ? val_get_0 : val_get}, start=int(starts), count=int(counts))"
  end
  def check_error
    "if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', &
        variable_name, ', ',  nf90_strerror(status)"
  end
  def function_string
    _string = <<EOF
 subroutine #{procedure_name}(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
   use netcdf, only: nf90_put_var, nf90_strerror
#endif
   use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
   use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   #{@type}, intent(in)#{@dimension} :: val
   integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
   integer :: fileid, varid, status, n2
   #{@offset ? "integer :: n, ndims" : ""}
   #{complex ? realval_declaration : nil}

#ifdef FORTRAN_NETCDF
   #{get_n2}
   allocate(starts(n2))
   allocate(counts(n2))
   allocate(offsets(n2))
   #{
if complex 
   2.times.map{|i| "
   call set_count(sfile, variable_name, \"r\", 1)
   call set_start(sfile, variable_name, \"r\", #{i+1})
   realval(1#{([",:"]*@dimsize).join('')}) = #{realimag[i]}(#{@dimsize==0 ? "val" : "val(#{dimension_spec})"})
   #{netcdf_inputs} 
   #{put_variable} 
   #{check_error}"
  }.join("\n\n")        
else "  
   #{netcdf_inputs}
   #{put_variable}
   #{check_error}"
end
   }

   deallocate(starts, offsets, counts)
#endif
 end subroutine #{procedure_name}
EOF
  end

end

begin
  4.times.map{|i|}
rescue
  puts "You appear to be running ruby 1.8.6 or lower... suggest you upgrade your ruby version!"
  class Integer
    def times(&block)
      if block
        (0...self).to_a.each{|i| yield(i)}
      else
        return  (0...self).to_a
      end
    end
  end
end

generators = []
generators_no_offset = []
['real', 'double precision', 'integer', 'character', 'complex', 'complex*16'].each do |type|
  (0..6).each do |dimsize|
    generators.push Generator.new(type, dimsize, true)
    generators_no_offset.push Generator.new(type, dimsize, false)
  end
end

string = <<EOF
module simpledataio_write
  implicit none
  private
  public :: write_variable_with_offset, write_variable

  interface write_variable_with_offset
#{generators.map{|g| "  "+g.interface_name}.join("\n")}
  end interface write_variable_with_offset

  interface write_variable
#{generators_no_offset.map{|g| "  "+g.interface_name}.join("\n")}
  end interface write_variable

contains

#{generators.map{|g| g.function_string}.join("\n")}

#{generators_no_offset.map{|g| g.function_string}.join("\n")}

end module simpledataio_write

EOF


File.open(ARGV[-1], 'w'){|file| file.puts string}
