module mp
  implicit none
  private

  public :: init_mp, finish_mp
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0, job
  public :: send, receive
  public :: barrier
  public :: scope, allprocs, subprocs
  public :: init_jobs

  integer, parameter :: nproc = 1, iproc = 0
  logical, parameter :: proc0 = .true.

  integer, parameter :: job = 0
  integer, parameter :: allprocs = 0, subprocs = 1

  interface broadcast
     module procedure broadcast_integer 
     module procedure broadcast_integer_array 

     module procedure broadcast_real    
     module procedure broadcast_real_array    

     module procedure broadcast_complex 
     module procedure broadcast_complex_array

     module procedure broadcast_logical 
     module procedure broadcast_logical_array 

     module procedure bcastfrom_integer 
     module procedure bcastfrom_integer_array 

     module procedure bcastfrom_real    
     module procedure bcastfrom_real_array    

     module procedure bcastfrom_complex 
     module procedure bcastfrom_complex_array 

     module procedure bcastfrom_logical 
     module procedure bcastfrom_logical_array 

     module procedure broadcast_character
     module procedure bcastfrom_character
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer
     module procedure sum_reduce_integer_array

     module procedure sum_reduce_real
     module procedure sum_reduce_real_array

     module procedure sum_reduce_complex
     module procedure sum_reduce_complex_array
  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer
     module procedure sum_allreduce_integer_array

     module procedure sum_allreduce_real
     module procedure sum_allreduce_real_array

     module procedure sum_allreduce_complex
     module procedure sum_allreduce_complex_array
  end interface

  interface max_reduce
     module procedure max_reduce_integer
     module procedure max_reduce_integer_array

     module procedure max_reduce_real
     module procedure max_reduce_real_array
  end interface

  interface max_allreduce
     module procedure max_allreduce_integer
     module procedure max_allreduce_integer_array

     module procedure max_allreduce_real
     module procedure max_allreduce_real_array
  end interface

  interface min_reduce
     module procedure min_reduce_integer
     module procedure min_reduce_integer_array

     module procedure min_reduce_real
     module procedure min_reduce_real_array
  end interface

  interface min_allreduce
     module procedure min_allreduce_integer
     module procedure min_allreduce_integer_array

     module procedure min_allreduce_real
     module procedure min_allreduce_real_array
  end interface

  interface send
     module procedure send_integer
     module procedure send_integer_array

     module procedure send_real
     module procedure send_real_array

     module procedure send_complex
     module procedure send_complex_array

     module procedure send_logical
     module procedure send_logical_array
  end interface

  interface receive
     module procedure receive_integer
     module procedure receive_integer_array

     module procedure receive_real
     module procedure receive_real_array

     module procedure receive_complex
     module procedure receive_complex_array

     module procedure receive_logical
     module procedure receive_logical_array
  end interface

contains

  subroutine init_mp
  end subroutine init_mp

  subroutine scope (focus)
    integer, intent (in) :: focus
  end subroutine scope

  subroutine finish_mp
  end subroutine finish_mp

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
  end subroutine broadcast_integer_array

  subroutine broadcast_real (x)
    implicit none
    real, intent (in out) :: x
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real, dimension (:), intent (in out) :: x
  end subroutine broadcast_real_array

  subroutine broadcast_complex (z)
    implicit none
    complex, intent (in out) :: z
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
  end subroutine broadcast_complex_array

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
  end subroutine broadcast_logical_array

  subroutine broadcast_character (s)
    implicit none
    character(*), intent (in out) :: s
  end subroutine broadcast_character

  subroutine bcastfrom_integer (i, src)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_integer

  subroutine bcastfrom_integer_array (i, src)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_integer_array

  subroutine bcastfrom_real (x, src)
    implicit none
    real, intent (in out) :: x
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_real

  subroutine bcastfrom_real_array (x, src)
    implicit none
    real, dimension (:), intent (in out) :: x
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_real_array

  subroutine bcastfrom_complex (z, src)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_complex

  subroutine bcastfrom_complex_array (z, src)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_complex_array

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_logical

  subroutine bcastfrom_logical_array (f, src)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_logical_array

  subroutine bcastfrom_character (s, src)
    implicit none
    character(*), intent (in out) :: s
    integer, intent (in) :: src
    if (src /= 0) call error ("broadcast from")
  end subroutine bcastfrom_character

  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine sum_reduce_integer_array

  subroutine sum_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine sum_reduce_complex_array

  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_complex (z)
    implicit none
    complex, intent (in out) :: z
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
  end subroutine sum_allreduce_complex_array

  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine max_reduce_integer_array

  subroutine max_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine max_reduce_real_array

  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
  end subroutine max_allreduce_real_array

  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine min_reduce_integer

  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine min_reduce_integer_array

  subroutine min_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    if (dest /= 0) call error ("reduce to")
  end subroutine min_reduce_real_array

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
  end subroutine min_allreduce_real_array

  subroutine barrier
  end subroutine barrier

  subroutine send_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_integer

  subroutine send_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_integer_array

  subroutine send_real (a, dest, tag)
    implicit none
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_real

  subroutine send_real_array (a, dest, tag)
    implicit none
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_real_array

  subroutine send_complex (z, dest, tag)
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_complex

  subroutine send_complex_array (z, dest, tag)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_complex_array

  subroutine send_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_logical

  subroutine send_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_logical_array

  subroutine send_character (s, dest, tag)
    implicit none
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    call error ("send")
  end subroutine send_character

  subroutine receive_integer (i, src, tag)
    implicit none
    integer, intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_integer

  subroutine receive_integer_array (i, src, tag)
    implicit none
    integer, dimension (:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_integer_array

  subroutine receive_real (a, src, tag)
    implicit none
    real, intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_real

  subroutine receive_real_array (a, src, tag)
    implicit none
    real, dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_real_array

  subroutine receive_complex (z, src, tag)
    implicit none
    complex, intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_complex

  subroutine receive_complex_array (z, src, tag)
    implicit none
    complex, dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_complex_array

  subroutine receive_logical (f, src, tag)
    implicit none
    logical, intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_logical

  subroutine receive_logical_array (f, src, tag)
    implicit none
    logical, dimension (:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_logical_array

  subroutine receive_character (s, src, tag)
    implicit none
    character(*), intent (out) :: s
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    call error ("receive")
  end subroutine receive_character

  subroutine init_jobs(ncolumns, group0, ierr)
    integer, intent (in) :: ncolumns
    integer, dimension(0:) :: group0
    integer :: ierr

    if (ncolumns /= 1) call error ("jobs")

  end subroutine init_jobs

  subroutine error (msg)
    implicit none
    character(*), intent (in) :: msg

    print *, "mp error: "//msg
    call coredump
    stop
  end subroutine error

  subroutine coredump
    real, dimension (:), allocatable :: a
    deallocate (a)
  end subroutine coredump

end module mp
