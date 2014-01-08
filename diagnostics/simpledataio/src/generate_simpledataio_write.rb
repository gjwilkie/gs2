
class Generator
	def initialize(type, dimsize)
		@dimsize = dimsize
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
		"write_variable_#{@type.gsub(/[* ]/, '_')}_#{@dimsize}"
	end
	def val_get
		complex ? 
		"realval(1:,#{@dimsize.times.map{|i| "starts(#{i+2}):"}.join(",")})"
		:
		"val(#{@dimsize.times.map{|i| "starts(#{i+1}):"}.join(",")})"
	end
	def val_get_0
		complex ?
    "realval"
			:
    "(/val/)"
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
	def function_string
		string = <<EOF
 subroutine #{procedure_name}(sfile, variable_name, val)
   use netcdf
	 use simpledataio
   type(sdatio_file), intent(in) :: sfile
   character(*), intent(in) :: variable_name
   #{@type}, intent(in)#{@dimension} :: val
   integer, dimension(:), allocatable :: starts, counts
   integer :: fileid, varid, status, n
	 #{complex ? realval_declaration : nil}

   call number_of_unlimited_dimensions(sfile, variable_name, n)
   allocate(starts(n+#@dimsize#{complex ? "+1" : nil}))
   allocate(counts(n+#@dimsize#{complex ? "+1" : nil}))
	 #{complex ?
 	
	2.times.map{|i| "
	 call set_count(sfile, variable_name, \"r\", 1)
	 call set_start(sfile, variable_name, \"r\", #{i+1})
	 realval(1#{([",:"]*@dimsize).join('')}) = #{realimag[i]}(#{@dimsize==0 ? "val" : "val(#{dimension_spec})"})
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   #{@dimsize == 0 ? val_get_0 : val_get}, start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)"
	}.join("\n\n")

	:

  " 
   call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts)
   status =  nf90_put_var(fileid, varid+1, &
		   #{@dimsize == 0 ? "(/val/)" : val_get}, start=starts, count=counts)
   if (.not. status .eq. 0) write (*,*) 'Error writing variable: ', variable_name, ', ',  nf90_strerror(status)"
	 }


 end subroutine #{procedure_name}
EOF
  end

end

generators = []
['real', 'double precision', 'integer', 'character', 'complex', 'complex*16'].each do |type|
	(0..6).each do |dimsize|
		generators.push Generator.new(type, dimsize)
	end
end

string = <<EOF

module simpledataio_write

interface write_variable
#{generators.map{|g| "  module procedure " + g.procedure_name}.join("\n")}
end interface write_variable

contains

#{generators.map{|g| g.function_string}.join("\n")}

end module simpledataio_write

EOF


puts string
