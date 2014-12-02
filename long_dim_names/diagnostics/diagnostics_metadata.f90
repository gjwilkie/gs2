!> A module for for writing help, metadata and the input file 
!! to the netcdf file. 
!! 
!! To extract the input file use this command
!!    ncdump <outputfile> -v input_file| \
!!      grep 'input_file = ' | sed s/.\ \;\$// | sed s/input_file.....// | sed \
!!       s/\\\\\\\"/\\\"/g | sed s/\\\(\\\\\\\\n\ \\\)*\\\\\\\\n/\\n/g 
module diagnostics_metadata
  use diagnostics_config, only: diagnostics_type
  implicit none
  character, dimension(100000) :: inputfile_array
  integer :: inputfile_length

contains
  subroutine write_metadata(gnostics)
    use simpledataio, only: add_metadata
    use simpledataio, only: add_standard_metadata
    use runtime_tests, only: build_indentifier
    use runtime_tests, only: get_svn_rev
    type(diagnostics_type), intent(in) :: gnostics
    character (20) :: datestamp, timestamp, timezone
    character (31) :: strtime
    !integer :: inttime
    integer, dimension(8) :: values

    call add_metadata(gnostics%sfile, "title", &
      new_line("")//&
      "=================================="&
      //new_line("")//&
      "          GS2 Output File"&
      //new_line("")//&
      "==================================")
    call add_metadata(gnostics%sfile, "file_description", &
      new_line("")//"This file contains output from the gyrokinetic &
        &flux tube code GS2. "//new_line("")//&
      " GS2 is a code which evolves the gyrokinetic equation, "//new_line("")//&
      " which describes the behaviour of a strongly magnetized plasma. "//new_line("")// &
      " Among other information this file may contain: "//new_line("")//&
      " (*) Information about the geometry of the magnetic field. "//new_line("")//&
      " (*) Values of the electromagnetic fields. "//new_line("")//&
      " (*) Values of linear growth rates. "//new_line("")//&
      " (*) Values of turbulent fluxes. "//new_line("")//&
      " (*) Values of moments such as density, parallel flow and temperature. "//new_line("")//&
      " ."//new_line("")//&
      " The details of what is contained are controlled by the namelist &
        diagnostics_config. "//new_line(""))
    call add_metadata(gnostics%sfile, "data_description", &
      new_line("")//"Most data consists of physical quantities given as "//new_line("")//&
      " a function of one more of the five dimensions in the GK eqn: "//new_line("")//&
      " theta, kx, ky, energy and lamda. "//new_line("")// &
      " Data is typically double precision real numbers.  "//new_line("")//&
      " Complex numbers are handled by having an extra dimension ri. "//new_line(""))
    call add_metadata(gnostics%sfile, "normalization_description", &
      new_line("")//"Quantities in this file are given in "//new_line("")//&
      " dimensionless form.  Dimensional quantities can  "//new_line("")//&
      " be reconstructed from these dimensionless "//new_line("")//&
      " quantities by multiplying by the expression "//new_line("")// &
      " given in the units attribute. This expression"//new_line("")//&
      " is constructed from one or more of the following "//new_line("")//&
      " normalising quantities: "//new_line("")//&
      " (*) a_ref: the normalising length "//new_line("")//&
      "       (this is half-diameter of "//new_line("")//&
      "       the LCFS for numerical equilibria "//new_line("")//&
      "       but has no physical meaning  "//new_line("")//&
      "       in Miller, circular and slab geometries) "//new_line("")//&
      " (*) B_ref: the normalising field "//new_line("")//&
      "       (this is the field on the magnetic "//new_line("")//&
      "       axis for numerical equilibria "//new_line("")//&
      "       but has no physical meaning in "//new_line("")//&
      "       Miller, circular and slab geometries) "//new_line("")//&
      " (*) Properties of the reference species (a "//new_line("")//&
      "       hypothetical species whose properties  "//new_line("")//&
      "       may be equal to one or none of the species "//new_line("")//&
      "       in the simulation). "//new_line("")//&
      "    (*) T_ref: the reference temperature. "//new_line("")//&
      "    (*) n_ref: the reference density. "//new_line("")//&
      "    (*) n_ref: the reference mass. "//new_line("")//&
      "    (*) Z_ref: the reference charge (which is always equal to the proton charge). "//new_line("")//&
      " (*) vth_ref = sqrt(2 T_ref/m_ref). "//new_line("")//&
      " (*) rho_ref = vth_ref / (Z_ref B_ref / m_ref c) . "//new_line(""))

    call add_metadata(gnostics%sfile, "gs2_help", &
      new_line("")//"At the time of writing you can obtain "//new_line("")//&
      "help for using GS2 in the following"//new_line("")//&
      "places: "//new_line("")//&
      " (*) http://gyrokinetics.sourceforge.net/wiki"//new_line("")//&
      "       (User help: installing, running) "//new_line("")//&
      " (*) http://gyrokinetics.sourceforge.net/gs2_documentation/ "//new_line("")//& 
      "       (Doxygen documentation) "//new_line("")//&
      " (*) http://sourceforge.net/projects/gyrokinetics "//new_line("")//& 
      "       (Downloads, bug tracker) "//new_line(""))
    call add_metadata(gnostics%sfile, "input_file_extraction", new_line("")//&
      "A bash command for extracting the input file from this "//new_line("")//&
      "file is printed on the wiki and in the sourcefile "//new_line("")//&
      "diagnostics_metadata.f90")
    call add_metadata(gnostics%sfile, "svn_revision", &
      trim(get_svn_rev()))
    call add_metadata(gnostics%sfile, "build_indentifier", &
      trim(build_indentifier()))

    call add_standard_metadata(gnostics%sfile)
  end subroutine write_metadata

  subroutine read_input_file_and_get_size(gnostics)
    use file_utils, only: get_input_unit
    use mp, only: proc0, broadcast
    type(diagnostics_type), intent(in) :: gnostics
    character(len=1000) :: line
    integer :: inputunit
    integer :: ios, i

    ! Find out the size of the input file and read it in
    ! This is unnecessarily done twice, as this function 
    ! is called once to create the variable and once to
    ! write, but the overhead is trivial
    if (proc0) then
      call get_input_unit (inputunit)
      rewind (unit=inputunit)
      ios = 0
      inputfile_length = 0
      do while (ios .eq. 0)
        line = ''
        read(unit=inputunit, iostat=ios, fmt=('(A,$)')) &
          line

        !if (len(trim(line)) .gt. 1) then
          do i = 1, len(trim(line))
            inputfile_length = inputfile_length + 1
            inputfile_array(inputfile_length) = line(i:i)
          end do
          inputfile_length = inputfile_length + 1
          inputfile_array(inputfile_length) = " "
          inputfile_length = inputfile_length + 1
          inputfile_array(inputfile_length) = "\"
          inputfile_length = inputfile_length + 1
          inputfile_array(inputfile_length) = "n"
        !end if 
      end do
    end if
    call broadcast(inputfile_length)
    call broadcast(inputfile_array)
  end subroutine read_input_file_and_get_size

  !> Write the input file.
 
  subroutine write_input_file(gnostics)
    use file_utils, only: get_input_unit
    use diagnostics_create_and_write, only: create_and_write_variable
    use simpledataio, only: SDATIO_CHAR
    use diagnostics_dimensions, only: dim_string
    type(diagnostics_type), intent(in) :: gnostics

    call create_and_write_variable(gnostics, SDATIO_CHAR, "input_file", &
      trim(dim_string(gnostics%dims%input_file_dim)), &
       "Full text of the input file", "text", inputfile_array) 



  end subroutine write_input_file
end module diagnostics_metadata
