class Generator
  def initialize(type, dimsize)
    @dimsize = dimsize
    @type = type
    if dimsize==0
      @dimension = ""
    else
      @dimension = ", dimension(#{([":"]*dimsize).join(",")}) "
    end
  end
  def procedure_name
    "create_and_write_variable_#{@type.gsub(' ', '_')}_#{@dimsize}"
  end
  def val_get
    "val(#{@dimsize.times.map{|i| "starts(#{i+1}):"}.join(",")})"
  end
  def function_string
    _string = <<EOF
  subroutine #{procedure_name}(gnostics, variable_type, variable_name, dimension_list, variable_description, variable_units, val)
    use simpledataio, only: create_variable
    use simpledataio_write, only: write_variable
    use diagnostics_config, only: diagnostics_type
    type(diagnostics_type), intent(in) :: gnostics
    integer, intent(in) :: variable_type
    character(*), intent(in) :: variable_name
    character(*), intent(in) :: dimension_list
    character(*), intent(in) :: variable_description
    character(*), intent(in) :: variable_units
    #{@type.sub(/_/, '*')}, intent(in)#{@dimension} :: val
 
    if (gnostics%create) then 
       call create_variable(gnostics%sfile, variable_type, variable_name, trim(dimension_list), variable_description, variable_units)
    end if

    if (gnostics%create .or. .not. gnostics%wryte) return

    call write_variable(gnostics%sfile, variable_name, val)
  end subroutine #{procedure_name}
EOF
  end

end

class GeneratorDistributed
  def initialize(type, dimsize)
    @dimsize = dimsize
    @type = type
    if dimsize==0
      @dimension = ""
    else
      @dimension = ", dimension(#{([":"]*dimsize).join(",")}) "
    end
  end
  def procedure_name
    "create_and_write_dstrb_fieldlike_variable_#{@type.gsub(' ', '_')}_#{@dimsize}"
  end
  def val_get
    "val(#{@dimsize.times.map{|i| "starts(#{i+1}):"}.join(",")})"
  end
  def function_string
    _string = <<EOF
  subroutine #{procedure_name}(gnostics, variable_type, variable_name, dimension_list, &
    variable_description, variable_units, val)
    use simpledataio, only: create_variable, set_start, set_count
    use simpledataio, only: set_independent, set_collective
    use simpledataio_write, only: write_variable, write_variable_with_offset
    use diagnostics_config, only: diagnostics_type
    use mp, only: mp_abort,barrier
    use file_utils, only: error_unit
    use fields_parallelization, only: field_k_local
    use kt_grids, only: naky, ntheta0
    type(diagnostics_type), intent(in) :: gnostics
    integer, intent(in) :: variable_type
    character(*), intent(in) :: variable_name
    character(*), intent(in) :: dimension_list
    character(*), intent(in) :: variable_description
    character(*), intent(in) :: variable_units
    integer :: xdim
    integer :: id, it, ik
    integer :: i1 !, i2, i3
    #{@type.sub(/_/, '*')}, intent(in)#{@dimension} :: val
    #{@type.sub(/_/, '*')} :: dummy
   
    !return 
    ! Find location of the x dimension
!<DD>This may need changing if we change the string used to identify kx and ky.
    xdim = index(dimension_list, "XY")

    if (xdim .eq. 0) then
       write(error_unit(), *) "The function create_and_write_dstrb_field_like_variable should &
            & only be called for arrays whose dimension list contains XY in that order"
       call mp_abort("")
    end if
 
    if (gnostics%create) then 
       call create_variable(gnostics%sfile, variable_type, variable_name, trim(dimension_list), variable_description, variable_units)
       if (gnostics%distributed) then
       end if
    end if


    if (gnostics%wryte) then
       if (.not.  gnostics%distributed) then
          call write_variable(gnostics%sfile, variable_name, val)
       else
          ! For some reason every process has to make at least
          ! one write to a variable with an infinite dimension.
          ! Here we make some dummy writes to satisfy that
!<DD>This will need changing if we allow multi-character dimension names
          do id = 1,len(trim(dimension_list))
             if (dimension_list(id:id) .eq. 't') cycle
             call set_count(gnostics%sfile, variable_name, dimension_list(id:id), 1)
             !call set_start(gnostics%sfile, variable_name, dimension_list(id:id), 1)
          end do
          call write_variable(gnostics%sfile, variable_name, dummy)
          do id = 1,len(trim(dimension_list))
             !! Reset the starts and counts
             if (dimension_list(id:id) .eq. 't') cycle
             call set_count(gnostics%sfile, variable_name, dimension_list(id:id), -1)
             !call set_start(gnostics%sfile, variable_name, dimension_list(id:id), -1)
          end do
          call barrier
          call set_count(gnostics%sfile, variable_name, "X", 1)
          call set_count(gnostics%sfile, variable_name, "Y", 1)
          call set_independent(gnostics%sfile, variable_name)
          do ik = 1,naky
             do it = 1,ntheta0
                if (field_k_local(it,ik)) then
                   call set_start(gnostics%sfile, variable_name, "X", it)
                   call set_start(gnostics%sfile, variable_name, "Y", ik)
                   ! Now we treat cases where X and Y are not the two most
                   ! slowly varying indices
                   if (xdim < #@dimsize - 1) then
                      call set_count(gnostics%sfile, variable_name, dimension_list(xdim+2:xdim+2), 1)
                      ! This loop will normally be over species
                      do i1 = 1,size(val, xdim+2)
                         call set_start(gnostics%sfile, variable_name, dimension_list(xdim+2:xdim+2), i1)
                         if (xdim < #@dimsize - 2) then
                            write (*,*) "Case with two dimensions to the right of X and Y not implemented"
                            !<DD>Should this be an mp_abort?
                            stop 1
                         else
                            call write_variable_with_offset(gnostics%sfile, variable_name, val)
                         end if
                      end do 
                   else
                      call write_variable_with_offset(gnostics%sfile, variable_name, val)
                   end if
                end if
             end do
          end do
          call set_collective(gnostics%sfile, variable_name)
       end if
    end if

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
distributed_generators = []
['real', 'integer', 'character', 'double precision', 'complex', 'complex_16'].each do |type| # 
  (0..6).each do |dimsize|
    generators.push Generator.new(type, dimsize)
    next if dimsize < 2
    distributed_generators.push GeneratorDistributed.new(type, dimsize)
  end
end

string = <<EOF

! DO NOT EDIT THIS FILE
! This file is automatically generated by generate_diagnostics_create_and_write

!> A module which contains two high level interfaces for writing 
!! to the netcdf file, one for variables which are local to all
!! processors, and one for variables which may be distributed in 
!! the same manner as the fields.
module diagnostics_create_and_write

  implicit none

  private

  !> Create and/or write the given variable depending on the values
  !! of the flags gnostics%create and gnostics%wryte
  public :: create_and_write_variable

  !> These are a set of subroutines for writing variables which 
  !! have the dimensions X and Y (for example, fields, moments or 
  !! fluxes) which may be distributed, ie. different XY combinations
  !! may be on different processors. The locality is determined through
  !! the function field_k_local.
  public :: create_and_write_distributed_fieldlike_variable

  interface create_and_write_variable
#{generators.map{|g| "     "+"module procedure " + g.procedure_name}.join("\n")}
  end interface create_and_write_variable

  interface create_and_write_distributed_fieldlike_variable
#{distributed_generators.map{|g| "     "+"module procedure " + g.procedure_name}.join("\n")}
  end interface create_and_write_distributed_fieldlike_variable

contains

#{generators.map{|g| g.function_string}.join("\n")}

#{distributed_generators.map{|g| g.function_string}.join("\n")}

end module diagnostics_create_and_write

EOF


File.open(ARGV[-1], 'w'){|f| f.puts string}
