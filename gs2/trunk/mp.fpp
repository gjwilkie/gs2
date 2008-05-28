# include "define.inc"

module mp
!
! Easier Fortran90 interface to the MPI Message Passing Library.
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
! Note: mp_mpi_r8.f90 is a version of mp_mpi.f90 to use when compiling 
! with -r8 (where the default real type is taken to be 8 bytes).  Just 
! replaced all occurances of MPI_REAL with MPI_DOUBLE_PRECISION and 
! MPI_COMPLEX with MPI_DOUBLE_COMPLEX.
!
# ifdef MPI
  use mpi
  ! TT: I experienced a problem on Ranger with mvapich-devel.
  ! TT: Compiler complained an inconsistent variable type for reorder below.
  ! TT: In this case, the problem was solved by commenting out the above line
  ! TT: and use mpif.h below.  (4/16/08)
# endif
  implicit none
  private

  public :: init_mp, finish_mp
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0, job
  public :: send, receive
  public :: barrier
! TT> no need to open them any more
!  public :: scope, allprocs, subprocs
! <TT
  public :: init_jobs
  public :: communicator

# ifdef MPI
!  include 'mpif.h'
  ! TT: I'll leave the above line as a comment.
  integer, pointer :: nproc
  integer, target :: ntot_proc, ngroup_proc

  integer, pointer :: iproc
  integer, target :: aproc, gproc

  logical, pointer :: proc0
  logical, target :: aproc0, gproc0

  integer, pointer :: communicator
  integer, target :: comm_all, comm_group

  integer :: job
  integer (kind(MPI_REAL)) :: mpireal, mpicmplx
# else
  integer, parameter :: nproc = 1, iproc = 0
  logical, parameter :: proc0 = .true.

  integer, parameter :: job = 0, communicator = -1
# endif
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
     
     module procedure send_character
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

     module procedure receive_character
  end interface

contains

  subroutine init_mp
# ifdef MPI
    use constants, only: pi, kind_rs, kind_rd
    implicit none
    integer :: ierror!, rank

    call mpi_init (ierror)
    call mpi_comm_size (mpi_comm_world, ntot_proc, ierror)
    call mpi_comm_rank (mpi_comm_world, aproc, ierror)
    comm_all = mpi_comm_world
    aproc0 = aproc == 0

    call scope (allprocs)

    if ( (kind(pi)==kind_rs) .and. (kind_rs/=kind_rd) ) then
       mpireal = MPI_REAL4
       mpicmplx = MPI_COMPLEX8
    else if (kind(pi)==kind_rd) then
       mpireal = MPI_REAL8
       mpicmplx = MPI_COMPLEX16
    else
       write (*,*) 'ERROR: precision mismatch in mpi'
    end if
# endif

  end subroutine init_mp

  subroutine scope (focus)

    integer, intent (in) :: focus

# ifdef MPI
    if (focus == allprocs) then
       communicator => comm_all
       nproc => ntot_proc
       iproc => aproc
       proc0 => aproc0
    else
       communicator => comm_group
       nproc => ngroup_proc
       iproc => gproc
       proc0 => gproc0
    end if
# endif

  end subroutine scope

  subroutine finish_mp
# ifdef MPI
    implicit none
    integer :: ierror

    call mpi_finalize (ierror)
# endif
  end subroutine finish_mp

! ************** broadcasts *****************************

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
# ifdef MPI
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, communicator, ierror)
# endif
  end subroutine broadcast_character

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, communicator, ierror)
# endif
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, communicator, ierror)
# endif
  end subroutine broadcast_integer_array

  subroutine broadcast_real (x)
    implicit none
    real, intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, 1, mpireal, 0, communicator, ierror)
# endif
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real, dimension (:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, communicator, ierror)
# endif
  end subroutine broadcast_real_array

  subroutine broadcast_complex (z)
    implicit none
    complex, intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, 1, mpicmplx, 0, communicator, ierror)
# endif
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, communicator, ierror)
# endif
  end subroutine broadcast_complex_array

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, communicator, ierror)
# endif
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, communicator, ierror)
# endif
  end subroutine broadcast_logical_array

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_logical

  subroutine bcastfrom_logical_array (f, src)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_logical_array

  subroutine bcastfrom_character (c, src)
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_character

  subroutine bcastfrom_integer (i, src)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_integer

  subroutine bcastfrom_integer_array (i, src)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_integer_array

  subroutine bcastfrom_real (x, src)
    implicit none
    real, intent (in out) :: x
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, 1, mpireal, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_real

  subroutine bcastfrom_real_array (x, src)
    implicit none
    real, dimension (:), intent (in out) :: x
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_real_array

  subroutine bcastfrom_complex (z, src)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, 1, mpicmplx, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_complex

  subroutine bcastfrom_complex_array (z, src)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, src, communicator, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_complex_array

! ************** reductions ***********************

  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_integer_array

  subroutine sum_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, mpireal, MPI_SUM, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), mpireal, MPI_SUM, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, 1, mpicmplx, MPI_SUM, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    complex, dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), mpicmplx, MPI_SUM, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_array

  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, communicator, ierror)
# endif
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, communicator, ierror)
# endif
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, mpireal, MPI_SUM, communicator, ierror)
# endif
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), mpireal, MPI_SUM, communicator, ierror)
# endif
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_complex (z)
    implicit none
    complex, intent (in out) :: z
# ifdef MPI
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, 1, mpicmplx, MPI_SUM, communicator, ierror)
# endif
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
# ifdef MPI
    complex, dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), mpicmplx, MPI_SUM, communicator, ierror)
# endif
  end subroutine sum_allreduce_complex_array

  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_integer_array

  subroutine max_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, mpireal, MPI_MAX, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), mpireal, MPI_MAX, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_real_array

  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, communicator, ierror)
# endif
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, communicator, ierror)
# endif
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, mpireal, MPI_MAX, communicator, ierror)
# endif
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), mpireal, MPI_MAX, communicator, ierror)
# endif
  end subroutine max_allreduce_real_array

  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_integer

  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_integer_array

  subroutine min_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, mpireal, MPI_MIN, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), mpireal, MPI_MIN, dest, communicator, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_real_array

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, communicator, ierror)
# endif
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, communicator, ierror)
# endif
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, mpireal, MPI_MIN, communicator, ierror)
# endif
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), mpireal, MPI_MIN, communicator, ierror)
# endif
  end subroutine min_allreduce_real_array

! ********************* barrier **********************

  subroutine barrier
# ifdef MPI
    implicit none
    integer :: ierror
    call mpi_barrier (communicator, ierror)
# endif
  end subroutine barrier

! ********************* sends **********************

  subroutine send_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_integer

  subroutine send_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_integer_array

  subroutine send_real (a, dest, tag)
    implicit none
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, mpireal, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_real

  subroutine send_real_array (a, dest, tag)
    implicit none
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), mpireal, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_real_array

  subroutine send_complex (z, dest, tag)
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, mpicmplx, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_complex

  subroutine send_complex_array (z, dest, tag)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), mpicmplx, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_complex_array

  subroutine send_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_logical

  subroutine send_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_logical_array

  subroutine send_character (s, dest, tag)
    implicit none
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send &
         (s, len(s), MPI_CHARACTER, dest, tagp, communicator, ierror)
# else
    call error ("send")
# endif
  end subroutine send_character

! ********************* receives  **********************

  subroutine receive_integer (i, src, tag)
    implicit none
# ifdef MPI
    integer, intent (out) :: i
# else
    integer :: i
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_integer

  subroutine receive_integer_array (i, src, tag)
    implicit none
# ifdef MPI
    integer, dimension (:), intent (out) :: i
# else
    integer, dimension (:) :: i
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_integer_array

  subroutine receive_real (a, src, tag)
    implicit none
# ifdef MPI
    real, intent (out) :: a
# else
    real :: a
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, mpireal, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_real

  subroutine receive_real_array (a, src, tag)
    implicit none
# ifdef MPI
    real, dimension (:), intent (out) :: a
# else
    real, dimension (:) :: a
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), mpireal, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_real_array

  subroutine receive_complex (z, src, tag)
    implicit none
# ifdef MPI
    complex, intent (out) :: z
# else
    complex :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, mpicmplx, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex

  subroutine receive_complex_array (z, src, tag)
    implicit none
# ifdef MPI
    complex, dimension (:), intent (out) :: z
# else
    complex, dimension (:) :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), mpicmplx, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex_array

  subroutine receive_logical (f, src, tag)
    implicit none
# ifdef MPI
    logical, intent (out) :: f
# else
    logical :: f
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_logical

  subroutine receive_logical_array (f, src, tag)
    implicit none
# ifdef MPI
    logical, dimension (:), intent (out) :: f
# else
    logical, dimension (:) :: f
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, communicator, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_logical_array

  subroutine receive_character (s, src, tag)
    implicit none
# ifdef MPI
    character(*), intent (out) :: s
# else
    character(*) :: s
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (s, len(s), MPI_CHARACTER, src, tagp, communicator, &
         status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_character

  subroutine init_jobs (ncolumns, group0, ierr)

    implicit none  
# ifdef MPI
!    integer, parameter :: reorder=1
    ! TT: I changed variable definition by assuming integer 1 corresponds to
    ! TT: logical .true. but I'm not sure if reorder is needed.
    ! TT: In any case this subroutine is only called when you use job fork.
    logical, parameter :: reorder=.true.
    integer :: ip, j, comm2d, id2d, ierr, nrows
# endif
    integer, intent(in) :: ncolumns
    integer, dimension(0:) :: group0
# ifndef MPI
    integer :: ierr
    if (ncolumns /= 1) call error ("jobs")
# else
    integer, parameter :: ndim=2
    integer, dimension(ndim) :: dims
    integer, dimension(0:ndim-1) :: coords1d, coords2d
    logical, dimension(0:ndim-1) :: belongs
    logical, dimension(ndim) :: period
    
    logical :: isroot
    
! calculate dimensions  mpi processor grid will have and check that 
! ncolumns*nrows = number of processes
    
    nrows = ntot_proc/ncolumns
    dims=(/ ncolumns, nrows /)     
    if(ntot_proc /= ncolumns*nrows) then
       ierr = 1
       if(aproc0) write(*,*) 'Number of processes must be divisible by number of groups'
       return
    endif
    ngroup_proc = nrows
    
    ! create 2d cartesian topology for processes
    
    period=(/ .false., .false. /)  !! no circular shift

    call mpi_cart_create(mpi_comm_world, ndim, dims, period, reorder, comm2d, ierr)
    call mpi_comm_rank(comm2d, id2d, ierr)
    call mpi_cart_coords(comm2d, id2d, ndim, coords2d, ierr)
    
! each processor knows which subgrid it is in from variable mpi_group
    job = coords2d(0)

! create 1d subgrids from 2d processor grid, variable belongs denotes
! whether processor grid is split by column or row

    belongs(1) = .true.    ! this dimension belongs to subgrid
    belongs(0) = .false.  

    call mpi_cart_sub(comm2d, belongs, comm_group, ierr)
    call mpi_comm_rank(comm_group, gproc, ierr)     
    call mpi_cart_coords(comm_group, gproc, 1, coords1d, ierr)
    gproc0 = (gproc == 0)
    
! find root process of each 1d subgrid and place in array group0 indexed 
! from 0 to subgrids-1
     
    j=1
    group0(0) = 0
    do ip = 1, ntot_proc-1
       if (proc0) then
          call receive (isroot, ip)
          if (isroot) then
             group0(j) = ip
             j=j+1
          end if
       else if (ip == aproc) then
          call send (gproc0, 0)
       end if
       call barrier
    end do

!let all processors have the group0 array
    call broadcast (group0)

! TT> brought down here from init_job_name in file_utils.fpp
    call scope (subprocs)
! <TT

# endif
  end subroutine init_jobs

# ifndef MPI
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
# endif

end module mp
