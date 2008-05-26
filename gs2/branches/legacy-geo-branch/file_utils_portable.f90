module file_utils
  implicit none
  private

  public :: init_file_utils
  ! subroutine init_file_utils (list, input, error, name)
  ! logical, intent (out) :: list
  ! logical, intent (in), optional :: input, error
  ! character(*), intent (in), optional :: name
  !   default: INPUT=.true., ERROR=.true., NAME="unknown"
  !   Set up run_name(s) and list_name for output files
  !   Open input file and strip comments, unless disabled with INPUT=.false.
  !   Open error output file, unless disabled with ERROR=.false.

  public :: init_job_name
  ! subroutine ...

  public :: finish_file_utils
  ! subroutine finish_file_utils
  !   Clean up files opened in init

  public :: run_name
  ! character(500) :: run_name
  !    Label for the run, taken from the command line

  public :: list_name
  ! character(500) :: list_name
  !    Label for the list, taken from the command line

  public :: input_unit
  ! function input_unit (nml)
  ! character(*), intent (in) :: nml
  ! integer :: input_unit
  !    Rewind the input file to start of namelist NML,
  !    and return its unit number

  public :: input_unit_exist
  ! function input_unit_exist (nml,exist)
  ! character(*), intent (in) :: nml
  ! integer :: input_unit
  !    Rewind the input file to start of namelist NML,
  !    and return its unit number, setexist=.true.
  !    If the namelist NML isn't found, set exist=.false.

  public :: error_unit
  ! function error_unit ()
  ! integer :: error_unit
  !    Return the error unit number

  public :: get_input_unit

  public :: open_output_file
  ! subroutine open_output_file (unit, ext)
  ! integer, intent (out) :: unit
  ! character (*), intent (in) :: ext
  !    Open a file with name made from the run_name with the EXT appended
  !    and return its unit number in UNIT

  public :: close_output_file
  ! subroutine close_output_file (unit)
  ! integer, intent (in) :: unit
  !    Close the file associated with UNIT from open_output_file

  public :: flush_output_file
  ! subroutine flush_output_file (unit)
  ! integer, intent (in) :: unit
  !    Close/open-append the file associated with UNIT from open_output_file

  public :: get_unused_unit
  ! subroutine get_unused_unit (unit)
  ! integer, intent (out) :: unit
  !    Return a unit number not associated with any file

  public :: get_indexed_namelist_unit
  ! subroutine get_indexed_namelist_unit (unit, nml, index)
  ! integer, intent (out) :: unit
  ! character (*), intent (in) :: nml
  ! integer, intent (in) :: index
  !    Copy namelist, NML // '_' // INDEX, from the input file to
  !    namelist, NML, in a temporary file, UNIT

  character(500), pointer :: run_name
  character(500), target :: arun_name, job_name
  character(len=500) :: list_name
  integer, save :: input_unit_no, error_unit_no
  integer, save, public :: num_input_lines

contains

  subroutine init_file_utils (list, input, error, name)
    implicit none
    logical :: list
    logical, intent (in), optional :: input, error
    character(*), intent (in), optional :: name
    logical :: inp, err

    if (.not. present (input)) then
       inp = .true.
    else
       inp = input
    end if
    if (.not. present (error)) then
       err = .true.
    else
       err = error
    end if
    if (.not. present (name)) then
       arun_name = "unknown"
    else
       arun_name = trim(name)
    end if

    call run_type (list)

    if (list) then
       list_name = arun_name
    else
       call init_run_name
       call init_error_unit (err)
       call init_input_unit (inp)
    end if

  end subroutine init_file_utils

  subroutine run_type (list)
    use command_line
    implicit none
    logical :: list
    integer :: l, ierr

    if (iargc() /= 0) then
       call cl_getarg (1, arun_name, 500, ierr)
       if (ierr /= 0) then
          print *, "Error getting run name."
       end if
    end if

    l = len_trim (arun_name)
    if (l>5) then
       list = arun_name(l-4:l) == ".list"
    else
       list = .false.
    end if

  end subroutine run_type

  subroutine init_run_name

    implicit none
    integer :: l

    l = len_trim (arun_name)
    if (l > 3 .and. arun_name(l-2:l) == ".in") then
       arun_name = arun_name(1:l-3)
    end if
    run_name => arun_name

  end subroutine init_run_name

!  subroutine init_job_name (njobs, group0, job_list)
  subroutine init_job_name (jobname)
!    use command_line
!    use mp
    implicit none
!    integer, intent (in) :: njobs
!    integer, intent (in), dimension(0:) :: group0
!    character (len=500), dimension(0:) :: job_list
    character (len=500), intent(in) :: jobname
    logical :: err = .true., inp = .true.
!    integer :: i

!    call scope (subprocs)
!    job_name = trim(job_list(job))
    job_name = trim(jobname)
    run_name => job_name

    call init_error_unit (err)
    call init_input_unit (inp)

  end subroutine init_job_name

  subroutine get_unused_unit (unit)
    implicit none
    integer, intent (out) :: unit
    character(20) :: read, write, readwrite
    unit = 50
    do
       inquire (unit=unit, read=read, write=write, readwrite=readwrite)
       if (read == "UNKNOWN" .and. write == "UNKNOWN" &
            .and. readwrite == "UNKNOWN") &
       then
          return
       end if
       unit = unit + 1
    end do
  end subroutine get_unused_unit

  subroutine open_output_file (unit, ext)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: ext
    character (500) :: hack
    call get_unused_unit (unit)
    hack=trim(run_name)//ext
    open (unit=unit, file=trim(hack), status="replace", action="write")
  end subroutine open_output_file

  subroutine close_output_file (unit)
    implicit none
    integer, intent (in) :: unit
    close (unit=unit)
  end subroutine close_output_file

  subroutine flush_output_file (unit, ext)
    implicit none
    integer, intent (in) :: unit
    character (*), intent (in) :: ext
    character (500) :: hack
    hack=trim(run_name)//ext
    close(unit=unit)
    open (unit=unit, file=trim(hack), status="old", action="write", position="append")
  end subroutine flush_output_file

  subroutine init_error_unit (open_it)
    implicit none
    logical, intent (in) :: open_it
    error_unit_no = 6
    if (run_name /= "unknown" .and. open_it) then
       call open_output_file (error_unit_no, ".error")
    end if
  end subroutine init_error_unit

  subroutine strip_comments (line)
    implicit none
    character(*), intent (in out) :: line
    logical :: in_single_quotes, in_double_quotes
    integer :: i, length

    length = len_trim(line)
    i = 1
    in_single_quotes = .false.
    in_double_quotes = .false.
    loop: do
       if (in_single_quotes) then
          if (line(i:i) == "'") in_single_quotes = .false.
       else if (in_double_quotes) then
          if (line(i:i) == '"') in_double_quotes = .false.
       else
          select case (line(i:i))
          case ("'")
             in_single_quotes = .true.
          case ('"')
             in_double_quotes = .true.
          case ("!")
             i = i - 1
             exit loop
          end select
       end if
       if (i >= length) exit loop
       i = i + 1
    end do loop
    line = line(1:i)
  end subroutine strip_comments

  subroutine init_input_unit (open_it)
    implicit none
    logical, intent (in) :: open_it
    integer :: in_unit, out_unit, iostat
    character(500) :: line

    ! for includes
    integer, parameter :: stack_size = 10
    integer, dimension (stack_size) :: stack
    integer :: stack_ptr

    if (.not. open_it) then
       input_unit_no = -1
       return
    end if

    call get_unused_unit (in_unit)
    open (unit=in_unit, file=trim(run_name)//".in", status="old", &
         action="read", iostat=iostat)
    if (iostat /= 0) then
       print "(a)", "Couldn't open input file: "//trim(run_name)//".in"
    end if

    call get_unused_unit (out_unit)
!    open (unit=out_unit, status="scratch", action="readwrite")
    open (unit=out_unit, file="."//trim(run_name)//".in")

    iostat = 0
    stack_ptr = 0
    num_input_lines = 0
    do
       read (unit=in_unit, fmt="(a)", iostat=iostat) line
       if (iostat /= 0) then
          if (stack_ptr <= 0) exit
          close (unit=in_unit)
          iostat = 0
          in_unit = stack(stack_ptr)
          stack_ptr = stack_ptr - 1
          cycle
       end if
       if (line(1:9) == "!include ") then
          if (stack_ptr >= stack_size) then
             print "(a)", "!include ignored: nesting too deep: "//trim(line)
             cycle
          end if
          stack_ptr = stack_ptr + 1
          stack(stack_ptr) = in_unit
          call get_unused_unit (in_unit)
          open (unit=in_unit, file=trim(line(10:)), status="old", &
                action="read", iostat=iostat)
          if (iostat /= 0) then
             print "(a)", "!include ignored: file unreadable: "//trim(line)
             in_unit = stack(stack_ptr)
             stack_ptr = stack_ptr - 1
             cycle
          end if
          cycle
       end if
       call strip_comments (line)
       write (unit=out_unit, fmt="(a)") trim(line)
       num_input_lines = num_input_lines + 1
    end do
    close (unit=in_unit)

    input_unit_no = out_unit
  end subroutine init_input_unit

  subroutine finish_file_utils
    implicit none
    if (input_unit_no > 0) then
       close (unit=input_unit_no)
       input_unit_no = -1
    end if
    if (error_unit_no > 0 .and. error_unit_no /= 6) then
       close (unit=error_unit_no)
       error_unit_no = -1
    end if
  end subroutine finish_file_utils

  function input_unit (nml)
    implicit none
    character(*), intent (in) :: nml
    integer :: input_unit, iostat
    character(500) :: line
    intrinsic adjustl, trim
    input_unit = input_unit_no
    if (input_unit_no > 0) then
       rewind (unit=input_unit_no)
       do
          read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
          if (iostat /= 0) then
             rewind (unit=input_unit_no)
             exit
          end if
          if (trim(adjustl(line)) == "&"//nml) then
             backspace (unit=input_unit_no)
             return
          end if
       end do
    end if
    write (unit=error_unit_no, fmt="('Couldn''t find namelist: ',a)") nml
    write (unit=*, fmt="('Couldn''t find namelist: ',a)") nml
  end function input_unit

  function input_unit_exist (nml,exist)
    implicit none
    character(*), intent (in) :: nml
    logical, intent(out) :: exist
    integer :: input_unit_exist, iostat
    character(500) :: line
    intrinsic adjustl, trim
    input_unit_exist = input_unit_no
    exist = .true.
    if (input_unit_no > 0) then
       rewind (unit=input_unit_no)
       do
          read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
          if (iostat /= 0) then
             rewind (unit=input_unit_no)
             exit
          end if
          if (trim(adjustl(line)) == "&"//nml) then
             backspace (unit=input_unit_no)
             return
          end if
       end do
    end if
    exist = .false.
  end function input_unit_exist

  function error_unit ()
    implicit none
    integer :: error_unit
    error_unit = error_unit_no
  end function error_unit

  subroutine get_input_unit (unit)
    implicit none
    integer, intent (out) :: unit

    unit = input_unit_no

  end subroutine get_input_unit

  subroutine get_indexed_namelist_unit (unit, nml, index)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: nml
    integer, intent (in) :: index
    character(500) :: line
    integer :: iunit, iostat, in_file
    logical :: exist

    call get_unused_unit (unit)
!    open (unit=unit, status="scratch", action="readwrite")
    open (unit=unit, file="."//trim(run_name)//".scratch")

    write (line, *) index
    line = nml//"_"//trim(adjustl(line))
    in_file = input_unit_exist(trim(line), exist)
    if (exist) then
       iunit = input_unit(trim(line))
    else
       return
    end if

    read (unit=iunit, fmt="(a)") line
    write (unit=unit, fmt="('&',a)") nml

    do
       read (unit=iunit, fmt="(a)", iostat=iostat) line
       if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
       write (unit=unit, fmt="(a)") trim(line)
    end do
    write (unit=unit, fmt="('/')")
    rewind (unit=unit)
  end subroutine get_indexed_namelist_unit

end module file_utils
