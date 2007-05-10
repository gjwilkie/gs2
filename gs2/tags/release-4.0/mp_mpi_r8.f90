module mp
!
! Easier Fortran90 interface to the MPI Message Passing Library.
!
! Note: mp_mpi_r8.f90 is a version of mp_mpi.f90 to use when compiling 
! with -r8 (where the default real type is taken to be 8 bytes).  Just 
! replaced all occurances of MPI_REAL with MPI_DOUBLE_PRECISION and 
! MPI_COMPLEX with MPI_DOUBLE_COMPLEX.
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
  implicit none
  private

  public :: init_mp, finish_mp
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0
  public :: send, receive
  public :: barrier

  integer :: nproc, iproc
  logical :: proc0

  interface broadcast
     module procedure broadcast_integer 
     module procedure broadcast_integer_array 
     module procedure broadcast_integer_array2
     module procedure broadcast_integer_array3
     module procedure broadcast_integer_array4

     module procedure broadcast_real    
     module procedure broadcast_real_array    
     module procedure broadcast_real_array2    
     module procedure broadcast_real_array3    
     module procedure broadcast_real_array4

     module procedure broadcast_complex 
     module procedure broadcast_complex_array
     module procedure broadcast_complex_array2
     module procedure broadcast_complex_array3
     module procedure broadcast_complex_array4

     module procedure broadcast_logical 
     module procedure broadcast_logical_array 
     module procedure broadcast_logical_array2
     module procedure broadcast_logical_array3
     module procedure broadcast_logical_array4

     module procedure bcastfrom_integer 
     module procedure bcastfrom_integer_array 
     module procedure bcastfrom_integer_array2 
     module procedure bcastfrom_integer_array3 
     module procedure bcastfrom_integer_array4

     module procedure bcastfrom_real    
     module procedure bcastfrom_real_array    
     module procedure bcastfrom_real_array2    
     module procedure bcastfrom_real_array3
     module procedure bcastfrom_real_array4

     module procedure bcastfrom_complex 
     module procedure bcastfrom_complex_array 
     module procedure bcastfrom_complex_array2 
     module procedure bcastfrom_complex_array3
     module procedure bcastfrom_complex_array4

     module procedure bcastfrom_logical 
     module procedure bcastfrom_logical_array 
     module procedure bcastfrom_logical_array2 
     module procedure bcastfrom_logical_array3
     module procedure bcastfrom_logical_array4

     module procedure broadcast_character
     module procedure bcastfrom_character
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer
     module procedure sum_reduce_integer_array
     module procedure sum_reduce_integer_array2
     module procedure sum_reduce_integer_array3
     module procedure sum_reduce_integer_array4

     module procedure sum_reduce_real
     module procedure sum_reduce_real_array
     module procedure sum_reduce_real_array2
     module procedure sum_reduce_real_array3
     module procedure sum_reduce_real_array4

     module procedure sum_reduce_complex
     module procedure sum_reduce_complex_array
     module procedure sum_reduce_complex_array2
     module procedure sum_reduce_complex_array3
     module procedure sum_reduce_complex_array4
  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer
     module procedure sum_allreduce_integer_array
     module procedure sum_allreduce_integer_array2
     module procedure sum_allreduce_integer_array3
     module procedure sum_allreduce_integer_array4

     module procedure sum_allreduce_real
     module procedure sum_allreduce_real_array
     module procedure sum_allreduce_real_array2
     module procedure sum_allreduce_real_array3
     module procedure sum_allreduce_real_array4

     module procedure sum_allreduce_complex
     module procedure sum_allreduce_complex_array
     module procedure sum_allreduce_complex_array2
     module procedure sum_allreduce_complex_array3
     module procedure sum_allreduce_complex_array4
  end interface

  interface max_reduce
     module procedure max_reduce_integer
     module procedure max_reduce_integer_array
     module procedure max_reduce_integer_array2
     module procedure max_reduce_integer_array3
     module procedure max_reduce_integer_array4

     module procedure max_reduce_real
     module procedure max_reduce_real_array
     module procedure max_reduce_real_array2
     module procedure max_reduce_real_array3
     module procedure max_reduce_real_array4
  end interface

  interface max_allreduce
     module procedure max_allreduce_integer
     module procedure max_allreduce_integer_array
     module procedure max_allreduce_integer_array2
     module procedure max_allreduce_integer_array3
     module procedure max_allreduce_integer_array4

     module procedure max_allreduce_real
     module procedure max_allreduce_real_array
     module procedure max_allreduce_real_array2
     module procedure max_allreduce_real_array3
     module procedure max_allreduce_real_array4
  end interface

  interface min_reduce
     module procedure min_reduce_integer
     module procedure min_reduce_integer_array
     module procedure min_reduce_integer_array2
     module procedure min_reduce_integer_array3
     module procedure min_reduce_integer_array4

     module procedure min_reduce_real
     module procedure min_reduce_real_array
     module procedure min_reduce_real_array2
     module procedure min_reduce_real_array3
     module procedure min_reduce_real_array4
  end interface

  interface min_allreduce
     module procedure min_allreduce_integer
     module procedure min_allreduce_integer_array
     module procedure min_allreduce_integer_array2
     module procedure min_allreduce_integer_array3
     module procedure min_allreduce_integer_array4

     module procedure min_allreduce_real
     module procedure min_allreduce_real_array
     module procedure min_allreduce_real_array2
     module procedure min_allreduce_real_array3
     module procedure min_allreduce_real_array4
  end interface

  interface send
     module procedure send_integer
     module procedure send_integer_array
     module procedure send_integer_array2
     module procedure send_integer_array3
     module procedure send_integer_array4

     module procedure send_real
     module procedure send_real_array
     module procedure send_real_array2
     module procedure send_real_array3
     module procedure send_real_array4

     module procedure send_complex
     module procedure send_complex_array
     module procedure send_complex_array2
     module procedure send_complex_array3
     module procedure send_complex_array4

     module procedure send_logical
     module procedure send_logical_array
     module procedure send_logical_array2
     module procedure send_logical_array3
     module procedure send_logical_array4
  end interface

  interface receive
     module procedure receive_integer
     module procedure receive_integer_array
     module procedure receive_integer_array2
     module procedure receive_integer_array3
     module procedure receive_integer_array4

     module procedure receive_real
     module procedure receive_real_array
     module procedure receive_real_array2
     module procedure receive_real_array3
     module procedure receive_real_array4

     module procedure receive_complex
     module procedure receive_complex_array
     module procedure receive_complex_array2
     module procedure receive_complex_array3
     module procedure receive_complex_array4

     module procedure receive_logical
     module procedure receive_logical_array
     module procedure receive_logical_array2
     module procedure receive_logical_array3
     module procedure receive_logical_array4
  end interface

contains

  subroutine init_mp
    implicit none
    include 'mpif.h'
    integer :: ierror, rank

    call mpi_init (ierror)
    call mpi_comm_size (mpi_comm_world, nproc, ierror)
    call mpi_comm_rank (mpi_comm_world, iproc, ierror)
    proc0 = iproc == 0
  end subroutine init_mp

  subroutine finish_mp
    implicit none
    include 'mpif.h'
    integer :: ierror

    call mpi_finalize (ierror)
  end subroutine finish_mp

! ************** broadcasts *****************************

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array

  subroutine broadcast_integer_array2 (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array2

  subroutine broadcast_integer_array3 (i)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array3

  subroutine broadcast_integer_array4 (i)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array4

  subroutine broadcast_real (x)
    implicit none
    real, intent (in out) :: x
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real, dimension (:), intent (in out) :: x
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array

  subroutine broadcast_real_array2 (x)
    implicit none
    real, dimension (:,:), intent (in out) :: x
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array2

  subroutine broadcast_real_array3 (x)
    implicit none
    real, dimension (:,:,:), intent (in out) :: x
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array3

  subroutine broadcast_real_array4 (x)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: x
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array4

  subroutine broadcast_complex (z)
    implicit none
    complex, intent (in out) :: z
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array

  subroutine broadcast_complex_array2 (z)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array2

  subroutine broadcast_complex_array3 (z)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array3

  subroutine broadcast_complex_array4 (z)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array4

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array

  subroutine broadcast_logical_array2 (f)
    implicit none
    logical, dimension (:,:), intent (in out) :: f
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array2

  subroutine broadcast_logical_array3 (f)
    implicit none
    logical, dimension (:,:,:), intent (in out) :: f
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array3

  subroutine broadcast_logical_array4 (f)
    implicit none
    logical, dimension (:,:,:,:), intent (in out) :: f
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array4

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical

  subroutine bcastfrom_logical_array (f, src)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array

  subroutine bcastfrom_logical_array2 (f, src)
    implicit none
    logical, dimension (:,:), intent (in out) :: f
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array2

  subroutine bcastfrom_logical_array3 (f, src)
    implicit none
    logical, dimension (:,:,:), intent (in out) :: f
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array3

  subroutine bcastfrom_logical_array4 (f, src)
    implicit none
    logical, dimension (:,:,:,:), intent (in out) :: f
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array4

  subroutine bcastfrom_character (c, src)
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_character

  subroutine bcastfrom_integer (i, src)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer

  subroutine bcastfrom_integer_array (i, src)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array

  subroutine bcastfrom_integer_array2 (i, src)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array2

  subroutine bcastfrom_integer_array3 (i, src)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array3

  subroutine bcastfrom_integer_array4 (i, src)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array4

  subroutine bcastfrom_real (x, src)
    implicit none
    real, intent (in out) :: x
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real

  subroutine bcastfrom_real_array (x, src)
    implicit none
    real, dimension (:), intent (in out) :: x
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array

  subroutine bcastfrom_real_array2 (x, src)
    implicit none
    real, dimension (:,:), intent (in out) :: x
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array2

  subroutine bcastfrom_real_array3 (x, src)
    implicit none
    real, dimension (:,:,:), intent (in out) :: x
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array3

  subroutine bcastfrom_real_array4 (x, src)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: x
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array4

  subroutine bcastfrom_complex (z, src)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex

  subroutine bcastfrom_complex_array (z, src)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array

  subroutine bcastfrom_complex_array2 (z, src)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array2

  subroutine bcastfrom_complex_array3 (z, src)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array3

  subroutine bcastfrom_complex_array4 (z, src)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent (in) :: src
    include 'mpif.h'
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array4

! ************** reductions ***********************

  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array

  subroutine sum_reduce_integer_array2 (i, dest)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array2

  subroutine sum_reduce_integer_array3 (i, dest)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array3

  subroutine sum_reduce_integer_array4 (i, dest)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3),size(i,4)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array4

  subroutine sum_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_real_array2 (a, dest)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array2

  subroutine sum_reduce_real_array3 (a, dest)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array3

  subroutine sum_reduce_real_array4 (a, dest)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3),size(a,4)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array4

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: dest
    include 'mpif.h'
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    include 'mpif.h'
    complex, dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array

  subroutine sum_reduce_complex_array2 (z, dest)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    integer, intent (in) :: dest
    include 'mpif.h'
    complex, dimension (size(z,1),size(z,2)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array2

  subroutine sum_reduce_complex_array3 (z, dest)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
    include 'mpif.h'
    complex, dimension (size(z,1),size(z,2),size(z,3)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array3

  subroutine sum_reduce_complex_array4 (z, dest)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
    include 'mpif.h'
    complex, dimension (size(z,1),size(z,2),size(z,3),size(z,4)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array4

  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    include 'mpif.h'
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_integer_array2 (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array2

  subroutine sum_allreduce_integer_array3 (i)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array3

  subroutine sum_allreduce_integer_array4 (i)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3),size(i,4)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array4

  subroutine sum_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
    include 'mpif.h'
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_real_array2 (a)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array2

  subroutine sum_allreduce_real_array3 (a)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array3

  subroutine sum_allreduce_real_array4 (a)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3),size(a,4)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array4

  subroutine sum_allreduce_complex (z)
    implicit none
    complex, intent (in out) :: z
    include 'mpif.h'
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
    include 'mpif.h'
    complex, dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array

  subroutine sum_allreduce_complex_array2 (z)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    include 'mpif.h'
    complex, dimension (size(z,1),size(z,2)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array2

  subroutine sum_allreduce_complex_array3 (z)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    include 'mpif.h'
    complex, dimension (size(z,1),size(z,2),size(z,3)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array3

  subroutine sum_allreduce_complex_array4 (z)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    include 'mpif.h'
    complex, dimension (size(z,1),size(z,2),size(z,3),size(z,4)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array4

  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array

  subroutine max_reduce_integer_array2 (i, dest)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array2

  subroutine max_reduce_integer_array3 (i, dest)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array3

  subroutine max_reduce_integer_array4 (i, dest)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3),size(i,4)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array4

  subroutine max_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array

  subroutine max_reduce_real_array2 (a, dest)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array2

  subroutine max_reduce_real_array3 (a, dest)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array3

  subroutine max_reduce_real_array4 (a, dest)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3),size(a,4)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array4

  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    include 'mpif.h'
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_integer_array2 (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array2

  subroutine max_allreduce_integer_array3 (i)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array3

  subroutine max_allreduce_integer_array4 (i)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3),size(i,4)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array4

  subroutine max_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
    include 'mpif.h'
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array

  subroutine max_allreduce_real_array2 (a)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array2

  subroutine max_allreduce_real_array3 (a)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array3

  subroutine max_allreduce_real_array4 (a)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3),size(a,4)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array4

  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer

  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array

  subroutine min_reduce_integer_array2 (i, dest)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array2

  subroutine min_reduce_integer_array3 (i, dest)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array3

  subroutine min_reduce_integer_array4 (i, dest)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    integer, intent (in) :: dest
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3),size(i,4)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array4

  subroutine min_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array

  subroutine min_reduce_real_array2 (a, dest)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array2

  subroutine min_reduce_real_array3 (a, dest)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array3

  subroutine min_reduce_real_array4 (a, dest)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3),size(a,4)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array4

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    include 'mpif.h'
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_integer_array2 (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array2

  subroutine min_allreduce_integer_array3 (i)
    implicit none
    integer, dimension (:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array3

  subroutine min_allreduce_integer_array4 (i)
    implicit none
    integer, dimension (:,:,:,:), intent (in out) :: i
    include 'mpif.h'
    integer, dimension (size(i,1),size(i,2),size(i,3),size(i,4)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array4

  subroutine min_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
    include 'mpif.h'
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array

  subroutine min_allreduce_real_array2 (a)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array2

  subroutine min_allreduce_real_array3 (a)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array3

  subroutine min_allreduce_real_array4 (a)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    include 'mpif.h'
    real, dimension (size(a,1),size(a,2),size(a,3),size(a,4)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array4

! ********************* barrier **********************

  subroutine barrier
    implicit none
    include 'mpif.h'
    integer :: ierror
    call mpi_barrier (MPI_COMM_WORLD, ierror)
  end subroutine barrier

! ********************* sends **********************

  subroutine send_integer (i, dest, tag)
    implicit none
    include 'mpif.h'
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer

  subroutine send_integer_array (i, dest, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array

  subroutine send_integer_array2 (i, dest, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:,:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array2

  subroutine send_integer_array3 (i, dest, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:,:,:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array3

  subroutine send_integer_array4 (i, dest, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:,:,:,:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array4

  subroutine send_real (a, dest, tag)
    implicit none
    include 'mpif.h'
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real

  subroutine send_real_array (a, dest, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array

  subroutine send_real_array2 (a, dest, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:,:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array2

  subroutine send_real_array3 (a, dest, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:,:,:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array3

  subroutine send_real_array4 (a, dest, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:,:,:,:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array4

  subroutine send_complex (z, dest, tag)
    implicit none
    include 'mpif.h'
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex

  subroutine send_complex_array (z, dest, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array

  subroutine send_complex_array2 (z, dest, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:,:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array2

  subroutine send_complex_array3 (z, dest, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:,:,:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array3

  subroutine send_complex_array4 (z, dest, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:,:,:,:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array4

  subroutine send_logical (f, dest, tag)
    implicit none
    include 'mpif.h'
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical

  subroutine send_logical_array (f, dest, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array

  subroutine send_logical_array2 (f, dest, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:,:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array2

  subroutine send_logical_array3 (f, dest, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:,:,:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array3

  subroutine send_logical_array4 (f, dest, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:,:,:,:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array4

  subroutine send_character (s, dest, tag)
    implicit none
    include 'mpif.h'
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send &
         (s, len(s), MPI_CHARACTER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_character

! ********************* receives  **********************

  subroutine receive_integer (i, src, tag)
    implicit none
    include 'mpif.h'
    integer, intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer

  subroutine receive_integer_array (i, src, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array

  subroutine receive_integer_array2 (i, src, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:,:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array2

  subroutine receive_integer_array3 (i, src, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:,:,:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array3

  subroutine receive_integer_array4 (i, src, tag)
    implicit none
    include 'mpif.h'
    integer, dimension (:,:,:,:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array4

  subroutine receive_real (a, src, tag)
    implicit none
    include 'mpif.h'
    real, intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real

  subroutine receive_real_array (a, src, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array

  subroutine receive_real_array2 (a, src, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:,:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array2

  subroutine receive_real_array3 (a, src, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:,:,:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array3

  subroutine receive_real_array4 (a, src, tag)
    implicit none
    include 'mpif.h'
    real, dimension (:,:,:,:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array4

  subroutine receive_complex (z, src, tag)
    implicit none
    include 'mpif.h'
    complex, intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex

  subroutine receive_complex_array (z, src, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array

  subroutine receive_complex_array2 (z, src, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:,:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array2

  subroutine receive_complex_array3 (z, src, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:,:,:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array3

  subroutine receive_complex_array4 (z, src, tag)
    implicit none
    include 'mpif.h'
    complex, dimension (:,:,:,:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array4

  subroutine receive_logical (f, src, tag)
    implicit none
    include 'mpif.h'
    logical, intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical

  subroutine receive_logical_array (f, src, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array

  subroutine receive_logical_array2 (f, src, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:,:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array2

  subroutine receive_logical_array3 (f, src, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:,:,:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array3

  subroutine receive_logical_array4 (f, src, tag)
    implicit none
    include 'mpif.h'
    logical, dimension (:,:,:,:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array4

end module mp
