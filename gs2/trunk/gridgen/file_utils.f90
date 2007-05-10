module file_utils
  implicit none
  private

  public :: init
  ! subroutine init (input, error)
  ! logical, intent (in), optional :: input, error
  !   default: INPUT=.true., ERROR=.true.
  !   Set up run_name for output files
  !   Open input file and strip comments, unless disabled with INPUT=.false.
  !   Open error output file, unless disabled with ERROR=.false.

  public :: finish
  ! subroutine finish
  !   Clean up files opened in init

  public :: run_name
  ! character(500), save :: run_name
  !    Label for the run, taken from the command line

  public :: input_unit
  ! function input_unit (nml)
  ! character(*), intent (in) :: nml
  ! integer :: input_unit
  !    Rewind the input file to start of namelist NML,
  !    and return its unit number

  public :: error_unit
  ! function error_unit ()
  ! integer :: error_unit
  !    Return the error unit number

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

  public :: unused_unit
  ! function unused_unit ()
  ! integer :: unused_unit
  !    Return a unit number not associated with any file

  character(500), save :: run_name
  integer, save :: input_unit_no, error_unit_no
contains
  subroutine init (input, error, name)
    implicit none
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
       run_name = "unknown"
    else
       run_name = name
    end if

    call init_run_name
    call init_error_unit (err)
    call init_input_unit (inp)
  end subroutine init

  subroutine init_run_name
    use command_line
    implicit none
    integer :: l, ierr
    integer, intrinsic :: len_trim
    if (iargc() /= 0) then
       call cl_getarg (1, run_name, 500, ierr)
       if (ierr /= 0) then
          print *, "Error getting run name."
          stop
       end if
    end if
    l = len_trim (run_name)
    if (l > 3 .and. run_name(l-2:l) == ".in") then
       run_name = run_name(1:l-3)
    end if
  end subroutine init_run_name

  subroutine init_error_unit (open_it)
    implicit none
    logical, intent (in) :: open_it
    error_unit_no = 6
    if (run_name /= "unknown" .and. open_it) then
       call open_output_file (error_unit_no, ".error")
    end if
  end subroutine init_error_unit

  subroutine init_input_unit (open_it)
    implicit none
    logical, intent (in) :: open_it
    integer :: in_unit, out_unit, iostat
    character(500) :: line
    if (.not. open_it) then
       input_unit_no = -1
       return
    end if
    in_unit = unused_unit()
    open (unit=in_unit, file=trim(run_name)//".in", status="old", &
         action="read", iostat=iostat)
    if (iostat /= 0) then
       print "(a)", "Couldn't open input file: "//trim(run_name)//".in"
       stop
    end if
    out_unit = unused_unit()
    open (unit=out_unit, file=trim(run_name)//"._in", status="replace", &
         action="write")
    iostat = 0
    do
       read (unit=in_unit, fmt="(a)", iostat=iostat) line
       if (iostat /= 0) exit
       write (unit=out_unit, fmt="(a)") trim(strip_comments(line))
    end do
    close (unit=in_unit)
    close (unit=out_unit)
    input_unit_no = unused_unit()
    open (unit=input_unit_no, file=trim(run_name)//"._in", status="old", &
         action="read")
  end subroutine init_input_unit

  function strip_comments (line)
    implicit none
    character(*), intent (in) :: line
    character(len=len(line)) :: strip_comments
    logical :: in_single_quotes, in_double_quotes
    integer :: i, length
    integer, intrinsic :: len_trim
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
    strip_comments = line(1:i)
  end function strip_comments

  subroutine finish
    implicit none
    if (input_unit_no > 0) then
       close (unit=input_unit_no)
       input_unit_no = -1
    end if
    if (error_unit_no > 0 .and. error_unit_no /= 6) then
       close (unit=error_unit_no)
       error_unit_no = -1
    end if
  end subroutine finish

  function input_unit (nml)
    implicit none
    character(*), intent (in) :: nml
    integer :: input_unit, iostat
    character(500) :: line
    intrinsic adjustl, trim
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
             exit
          end if
       end do
    end if
    input_unit = input_unit_no
  end function input_unit

  function error_unit ()
    implicit none
    integer :: error_unit
    error_unit = error_unit_no
  end function error_unit

  function unused_unit ()
    implicit none
    integer :: unused_unit
    integer :: unit
    character(20) :: read, write, readwrite
    unit = 10
    do
       inquire (unit=unit, read=read, write=write, readwrite=readwrite)
       if (read == "UNKNOWN" .and. write == "UNKNOWN" &
            .and. readwrite == "UNKNOWN") &
       then
          unused_unit = unit
          return
       end if
       unit = unit + 1
    end do
  end function unused_unit

  subroutine open_output_file (unit, ext)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: ext
    unit = unused_unit()
    open (unit=unit, file=trim(run_name)//ext, status="replace", &
         action="write")
  end subroutine open_output_file

  subroutine close_output_file (unit)
    implicit none
    integer, intent (in) :: unit
    close (unit=unit)
  end subroutine close_output_file
end module file_utils
