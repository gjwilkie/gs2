!> A module for for writing help, metadata and the input file 
!! to the netcdf file. 
!! 
!! To extract the input file use the bash utility scripts/extract_input_file
module diagnostics_metadata
  use diagnostics_config, only: diagnostics_type
  implicit none
  character, dimension(100000) :: inputfile_array
  integer :: inputfile_length

contains
  subroutine write_metadata(gnostics)
    use simpledataio, only: add_metadata
    use simpledataio, only: add_standard_metadata
    use runtime_tests, only: build_identifier
    use runtime_tests, only: get_svn_rev
    use run_parameters, only: user_comments
    type(diagnostics_type), intent(in) :: gnostics
    character (20) :: datestamp, timestamp, timezone
    character (31) :: strtime
    !integer :: inttime
    integer, dimension(8) :: values

    call add_metadata(gnostics%sfile, "title", &
      jline("")//&
      jline("=================================================")//&
      jline("                GS2 Output File")//&
      jline("=================================================")//&
      jline(""))
    call add_metadata(gnostics%sfile, "file_description", &
      jline("")//&
      jline(" This file contains output from the gyrokinetic")//&
      jline(" flux tube code GS2. GS2 is a code")//&
      jline(" which evolves the gyrokinetic equation, ")//&
      jline(" which describes the behaviour of a strongly")//& 
      jline(" magnetized plasma. ")// &
      jline(" Among other information this file may contain: ")//&
      jline(" (*) Values of geometric coefficients. ")//&
      jline(" (*) Values of the electromagnetic fields. ")//&
      jline(" (*) Values of linear growth rates. ")//&
      jline(" (*) Values of turbulent fluxes. ")//&
      jline(" (*) Values of moments such as density, ")//&
      jline("         parallel flow and temperature. ")//&
      jline(" The details of what is contained are controlled")//&
      jline(" by the namelist diagnostics_config. ")//&
      jline(""))
    call add_metadata(gnostics%sfile, "data_description", &
      jline("")//&
      jline(" Most data consists of physical quantities given")//&
      jline(" as a function of one more of the five dimensions")//&
      jline(" in the GK eqn: theta, kx, ky, energy and lamda. ")// &
      jline(" Data is typically double precision real numbers.  ")//&
      jline(" Complex numbers are handled by having an extra")//&
      jline(" dimension ri. ")//&
      jline(""))
    call add_metadata(gnostics%sfile, "normalization_description", &
      jline("")//&
      jline(" Quantities in this file are given in ")//&
      jline(" dimensionless form.  Dimensional quantities can")//&
      jline(" be reconstructed from these dimensionless ")//&
      jline(" quantities by multiplying by the expression ")// &
      jline(" given in the units attribute. This expression")//&
      jline(" is constructed from one or more of the following")//&
      jline(" normalising quantities: ")//&
      jline(" (*) a_ref: the normalising length ")//&
      jline("       (this is half-diameter of ")//&
      jline("       the LCFS for numerical equilibria ")//&
      jline("       but has no physical meaning  ")//&
      jline("       in Miller, circular and slab geometries) ")//&
      jline(" (*) B_ref: the normalising field ")//&
      jline("       (this is the field on the magnetic ")//&
      jline("       axis for numerical equilibria ")//&
      jline("       but has no physical meaning in ")//&
      jline("       Miller, circular and slab geometries) ")//&
      jline(" (*) Properties of the reference species (a ")//&
      jline("       hypothetical species whose properties  ")//&
      jline("       may be equal to one or none of the species ")//&
      jline("       in the simulation): ")//&
      jline("    (*) T_ref: the reference temperature ")//&
      jline("    (*) n_ref: the reference density ")//&
      jline("    (*) n_ref: the reference mass ")//&
      jline("    (*) Z_ref: the reference charge, which ")//&
      jline("          is always equal to the proton charge")//&
      jline(" (*) vth_ref = sqrt(2 T_ref/m_ref). ")//&
      jline(" (*) rho_ref = vth_ref / (Z_ref B_ref / m_ref c)")//&
      jline(""))

    call add_metadata(gnostics%sfile, "gs2_help", &
      jline("")//&
      jline(" At the time of writing you can obtain ")//&
      jline(" help for using GS2 in the following")//&
      jline(" places: ")//&
      jline(" (*) http://gyrokinetics.sourceforge.net/wiki")//&
      jline("       (User help: installing, running) ")//&
            new_line("")//& 
            " (*) http://gyrokinetics.sourceforge.net/gs2_documentation/ "//&
      jline("       (Doxygen documentation) ")//&
      jline(" (*) http://sourceforge.net/projects/gyrokinetics ")//& 
      jline("       (Downloads, bug tracker) ")//&
      jline(""))
    call add_metadata(gnostics%sfile, "input_file_extraction", &
      jline("")//&
      jline(" A bash utility for extracting the input file ")//&
      jline(" from this file is included in the scripts ")//&
      jline(" folder in the GS2 source. To invoke it:")//&
      jline(" $ ./extract_input_file <netcdf_file>")//&
      jline(""))
    call add_metadata(gnostics%sfile, "user_comments", &
      trim(user_comments))
    call add_metadata(gnostics%sfile, "svn_revision", &
      trim(get_svn_rev()))
    call add_metadata(gnostics%sfile, "build_identifier", &
      trim(build_identifier()))

    call add_standard_metadata(gnostics%sfile)
  end subroutine write_metadata

  function jline(inputline)
    use mp, only: mp_abort
    character(*), intent(in) :: inputline
    character(len=52) :: jline
    integer :: lenin
    lenin = len(trim(inputline))
    if (lenin .gt. 50) then
      call mp_abort("Long line: "//inputline, .true.)
    end if

    jline = ''

    jline(2:lenin+1) = inputline(1:lenin)
    jline(1:1) = new_line("")
    
  end function jline

  subroutine read_input_file_and_get_size(gnostics)
    use file_utils, only: get_input_unit
    use mp, only: proc0, broadcast
    type(diagnostics_type), intent(in) :: gnostics
    character(len=1000) :: line
    integer :: inputunit
    integer :: ios, i
    !Note : The input file used here is not directly the input file passed to GS2
    !as during init_input_unit (file_utils) the passed input file is processed to 
    !a) Strip out any comments
    !b) Deal with any included files.
    !It is this processed file (written to .<run_name>.in) which is stored.
    !We may wish to add the option to retain comments as these may be useful
    !when revisting old files to remind the user why the file is setup as it is.

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
          inputfile_array(inputfile_length) = "\" !" !These comment characters are for editors which think \ is escaping "
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
