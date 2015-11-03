# include "define.inc"

module mp
!
! <doc> Easier Fortran90 interface to the MPI Message Passing Library. </doc>
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
# ifndef MPIINC
  use mpi
#endif
  ! TT: I experienced a problem on Ranger with mvapich-devel.
  ! TT: Compiler complained an inconsistent variable type for reorder below.
  ! TT: In this case, the problem was solved by commenting out the above line
  ! TT: and use mpif.h below.  (4/16/08)
# endif
  implicit none
  private

  public :: init_mp, finish_mp
  public :: broadcast, sum_reduce, sum_allreduce, sum_allreduce_sub
  public :: broadcast_sub, sum_reduce_sub
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0, job
  public :: send, ssend, receive
  public :: send_init, recv_init, start_comm, free_request !For persistent comms
  public :: barrier
  public :: comm_type
  public :: waitany
  public :: test_driver_flag
! JH> new abort method
  public :: mp_abort
! <JH
! MAB> needed by Trinity
  public :: scope, allprocs, subprocs
  public :: all_to_group, group_to_all
  public :: trin_flag
! <MAB
!<DD> added
  public :: split
  public :: free_comm
  public :: waitall, testall
  public :: nbsend
  public :: nbrecv
  public :: nproc_comm
  public :: rank_comm
  public :: allgatherv, allgather
!</DD>
  public :: init_jobs
  public :: mp_comm
  public :: mp_info


  ! <EGH needed for functional_tests
  public :: grp0
  !EGH>

# ifdef MPI
# ifdef MPIINC
! CMR: defined MPIINC for machines where need to include mpif.h
  include 'mpif.h'
#endif
  integer, pointer :: nproc
  integer, target :: ntot_proc, ngroup_proc

  integer, pointer :: iproc
  integer, target :: aproc, gproc

  logical, pointer :: proc0
  logical, target :: aproc0, gproc0

  integer, pointer :: mp_comm
  integer, target :: comm_all, comm_group

  integer, parameter :: mp_info = MPI_INFO_NULL

  integer :: job = 0
  integer (kind(MPI_REAL)) :: mpireal, mpicmplx
# else
  integer, parameter :: nproc = 1, iproc = 0
  logical, parameter :: proc0 = .true.

  integer, parameter :: mp_info = -1
  integer, parameter :: job = 0, mp_comm = -1
# endif
  integer, parameter :: allprocs = 0, subprocs = 1

! needed for Trinity -- MAB
  integer, dimension (:), allocatable :: grp0
  logical :: trin_flag = .false.

  logical :: test_driver_flag = .false.

  !>A simple object for storing details of a communicator
  type comm_type
     integer :: id=-1 !The communicator id, used in calls to MPI routines
     integer :: iproc=-1 !The procs local rank
     integer :: nproc=-1 !The total number of processors in the communicator
     logical :: proc0=.false. !Is iproc equal to 0?
  end type comm_type

!<DD>
  interface waitall
     module procedure waitall_stat
     module procedure waitall_nostat
  end interface

  interface testall
     module procedure testall_stat
     module procedure testall_nostat
     module procedure testall_nostat_noflag
  end interface
     
!<DD>
  interface nbsend
     module procedure nbsend_real_array
     module procedure nbsend_real_array_count
     module procedure nbsend_complex_array
     module procedure nbsend_complex_array_sub
     module procedure nbsend_complex_array_count
  end interface

!<DD>
  interface nbrecv
     module procedure nbrecv_real_array
     module procedure nbrecv_real_array_count
     module procedure nbrecv_complex_array
     module procedure nbrecv_complex_array_sub
     module procedure nbrecv_complex_array_count
  end interface

!<DD>
  interface send_init
     module procedure send_init_complex_array
  end interface

!<DD>
  interface recv_init
     module procedure recv_init_complex_array
  end interface

!<DD>
  interface start_comm
     module procedure start_persist
     module procedure startall_persist
  end interface

!<DD>
  interface free_request
     module procedure free_handle_persist
     module procedure free_handles_persist
  end interface

!<DD>
  interface split
     module procedure split_nokey
     module procedure split_key
     module procedure split_nokey_to_commtype
     module procedure split_key_to_commtype
     module procedure split_nokey_to_commtype_sub
     module procedure split_key_to_commtype_sub
  end interface

  interface free_comm
     module procedure free_comm_id
     module procedure free_comm_type
  end interface

!<LA> gather could be faster that gatherv, it's worth using it when possible
interface allgather
    module procedure allgather_integer_array_1to1
end interface
!<DD>
  interface allgatherv
     module procedure allgatherv_complex_array_1to1
     module procedure allgatherv_complex_array_1to3
     module procedure allgatherv_complex_array_1to1_sub
     module procedure allgatherv_complex_array_1to3_sub
  end interface

  interface broadcast
     module procedure broadcast_integer 
     module procedure broadcast_integer_array 
     module procedure broadcast_integer_2array 

     module procedure broadcast_real    
     module procedure broadcast_real_array   
     module procedure broadcast_real_2array 
     module procedure broadcast_real_3array 
     module procedure broadcast_real_4array 

     module procedure broadcast_complex 
     module procedure broadcast_complex_array
     module procedure broadcast_complex_2array
     module procedure broadcast_complex_3array
     module procedure broadcast_complex_4array

     module procedure broadcast_logical 
     module procedure broadcast_logical_array 
     module procedure broadcast_logical_2array 

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

  interface broadcast_sub
     module procedure bcastfrom_complex_array_sub
     module procedure bcastfrom_complex_2array_sub  
     module procedure bcastfrom_complex_3array_sub  
     module procedure bcastfrom_complex_4array_sub  
  end interface

  interface sum_reduce_sub
     module procedure sum_reduce_complex_array_sub
     module procedure sum_reduce_complex_4array_sub
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer
     module procedure sum_reduce_integer_array

     module procedure sum_reduce_real
     module procedure sum_reduce_real_array
     module procedure sum_reduce_real_2array
     module procedure sum_reduce_real_3array
     module procedure sum_reduce_real_4array
     module procedure sum_reduce_real_5array

     module procedure sum_reduce_complex
     module procedure sum_reduce_complex_array
     module procedure sum_reduce_complex_2array
     module procedure sum_reduce_complex_3array
     module procedure sum_reduce_complex_4array
     module procedure sum_reduce_complex_5array
  end interface
  
  !KDN 100526: Allows summing into alternate variable
  !rather than overwriting local data
  interface sum_reduce_alt
     module procedure sum_reduce_alt_complex_3array
  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer
     module procedure sum_allreduce_integer_array
     module procedure sum_allreduce_integer_2array

     module procedure sum_allreduce_real
     module procedure sum_allreduce_real_array
     module procedure sum_allreduce_real_2array
     module procedure sum_allreduce_real_3array
     module procedure sum_allreduce_real_4array
     module procedure sum_allreduce_real_5array

     module procedure sum_allreduce_complex
     module procedure sum_allreduce_complex_array
     module procedure sum_allreduce_complex_2array
     module procedure sum_allreduce_complex_3array
     module procedure sum_allreduce_complex_4array
     module procedure sum_allreduce_complex_5array
  end interface

  interface sum_allreduce_sub
     module procedure sum_allreduce_sub_integer
     module procedure sum_allreduce_sub_integer_array

     module procedure sum_allreduce_sub_real
     module procedure sum_allreduce_sub_real_array
     module procedure sum_allreduce_sub_real_2array
     module procedure sum_allreduce_sub_real_3array
     module procedure sum_allreduce_sub_real_4array
     module procedure sum_allreduce_sub_real_5array

     module procedure sum_allreduce_sub_complex
     module procedure sum_allreduce_sub_complex_array
     module procedure sum_allreduce_sub_complex_2array
     module procedure sum_allreduce_sub_complex_3array
     module procedure sum_allreduce_sub_complex_4array
     module procedure sum_allreduce_sub_complex_5array
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

     module procedure nonblocking_send_complex_array

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
     module procedure receive_complex_2array
     module procedure nonblocking_receive_complex_array

     module procedure receive_logical
     module procedure receive_logical_array

     module procedure receive_character
  end interface

! MAB> needed for Trinity
! synchronous sends
  interface ssend
     module procedure ssend_integer
     module procedure ssend_integer_array

     module procedure ssend_real
     module procedure ssend_real_array

     module procedure ssend_complex
     module procedure ssend_complex_array
     module procedure ssend_complex_2array

     module procedure ssend_logical
     module procedure ssend_logical_array
  end interface

! send stuff from global proc0 to group proc0s
  interface all_to_group
     module procedure all_to_group_real
     module procedure all_to_group_real_array
  end interface

! send stuff from group proc0s to global proc0
  interface group_to_all
     module procedure group_to_all_real
     module procedure group_to_all_real_array
  end interface
! <MAB
  
  interface barrier
     module procedure barrier_nocomm
     module procedure barrier_comm
  end interface
contains

!MAB> changed to allow for calling within trinity (coupled flux tubes)
!  subroutine init_mp
  subroutine init_mp (comm_in)
# ifdef MPI
    use constants, only: pi, kind_rs, kind_rd
    use shm_mpi3
    use file_utils, only: error_unit
    implicit none
# endif
    integer, intent (in), optional :: comm_in ! MAB
# ifdef MPI
    integer :: ierror!, rank
    logical :: init ! MAB

    call mpi_initialized (init, ierror)  ! MAB
    if (.not. init) call mpi_init (ierror)  ! MAB
!    call mpi_init (ierror)
!MAB>
    if (present(comm_in)) then
       comm_all = comm_in
       if (.not. test_driver_flag) trin_flag = .true. 
       !EGH  we can pass in a communicator if we want to run tests
    else
       comm_all = mpi_comm_world
    end if
    call mpi_comm_size (comm_all, ntot_proc, ierror)
    call mpi_comm_rank (comm_all, aproc, ierror)
!    call mpi_comm_size (mpi_comm_world, ntot_proc, ierror)
!    call mpi_comm_rank (mpi_comm_world, aproc, ierror)
!    comm_all = mpi_comm_world
!<MAB
    aproc0 = aproc == 0

    call scope (allprocs)
    call shm_init(comm_all)

    if ( (kind(pi)==kind_rs) .and. (kind_rs/=kind_rd) ) then
       mpireal = MPI_REAL
       mpicmplx = MPI_COMPLEX
    else if (kind(pi)==kind_rd) then
       mpireal = MPI_DOUBLE_PRECISION
       mpicmplx = MPI_DOUBLE_COMPLEX
    else
       write (error_unit(),*) 'ERROR: precision mismatch in mpi'
    end if
# endif

  end subroutine init_mp

!<DD>How many procs in passed communicator
  subroutine nproc_comm(comm,nproc)
    integer, intent(in) :: comm
    integer, intent(out) :: nproc
#ifdef MPI
    integer :: ierror
    call mpi_comm_size(comm,nproc,ierror)
#else
    nproc=1
#endif
  end subroutine nproc_comm

!<DD>What's rank of proc in passed communicator
  subroutine rank_comm(comm,rank)
    integer, intent(in) :: comm
    integer, intent(out) :: rank
#ifdef MPI
    integer :: ierror
    call mpi_comm_rank(comm,rank,ierror)
#else
    rank=0
#endif
  end subroutine rank_comm

  subroutine scope (focus)

    integer, intent (in) :: focus

# ifdef MPI
    if (focus == allprocs) then
       mp_comm => comm_all
       nproc => ntot_proc
       iproc => aproc
       proc0 => aproc0
    else
       mp_comm => comm_group
       nproc => ngroup_proc
       iproc => gproc
       proc0 => gproc0
    end if
# endif

  end subroutine scope

  subroutine finish_mp
# ifdef MPI
    use shm_mpi3, only : shm_clean
    implicit none
    integer :: ierror
    
    call shm_clean
    call mpi_finalize (ierror)
# endif
  end subroutine finish_mp

! ************** allgathers *****************************
  subroutine allgather_integer_array_1to1(arr,count,out,recvcnts)
    !<DD>A subroutine to do a allgatherv operation, sending recvcnts(iproc)
    !data from the iproc'th processor to all others starting at arr(start).
    implicit none
    integer, dimension(:), intent(in) :: arr  !The data to gather
    integer, intent(in) :: count !How much data to gather, <=SIZE(arr)
    integer, dimension(:), intent(out) :: out !The gathered data
    integer, intent(in) :: recvcnts !Array detailing how much data to expect from each proc
# ifdef MPI
    integer :: ierror
    !Do the gather
    call mpi_allgather(arr,count,MPI_INTEGER,out,recvcnts,&
         MPI_INTEGER,mp_comm,ierror)
# else
    out=RESHAPE(arr,SHAPE(out))    
#endif
  end subroutine allgather_integer_array_1to1


  subroutine allgatherv_complex_array_1to1(arr,count,out,recvcnts,displs)
    !<DD>A subroutine to do a allgatherv operation, sending recvcnts(iproc)
    !data from the iproc'th processor to all others starting at arr(start).
    implicit none
    complex, dimension(:), intent(in) :: arr  !The data to gather
    integer, intent(in) :: count !How much data to gather, <=SIZE(arr)
    complex, dimension(:), intent(out) :: out !The gathered data
    integer, dimension(:), intent(in) :: recvcnts !Array detailing how much data to expect from each proc
    integer, dimension(:), intent(in) :: displs !Array detailing offset in array where gathered data is to be stored
# ifdef MPI
    integer :: ierror
    !Do the gather
    call mpi_allgatherv(arr,count,mpicmplx,out,recvcnts,displs,&
         mpicmplx,mp_comm,ierror)
# else
    out=RESHAPE(arr,SHAPE(out))    
#endif
  end subroutine allgatherv_complex_array_1to1

  subroutine allgatherv_complex_array_1to3(arr,count,out,recvcnts,displs)
    !<DD>A subroutine to do a allgatherv operation, sending recvcnts(iproc)
    !data from the iproc'th processor to all others starting at arr(start).
    implicit none
    complex, dimension(:), intent(in) :: arr  !The data to gather
    integer, intent(in) :: count !How much data to gather, <=SIZE(arr)
    complex, dimension(:,:,:), intent(out) :: out !The gathered data
    integer, dimension(:), intent(in) :: recvcnts !Array detailing how much data to expect from each proc
    integer, dimension(:), intent(in) :: displs !Array detailing offset in array where gathered data is to be stored
# ifdef MPI
    integer :: ierror
    !Do the gather
    call mpi_allgatherv(arr,count,mpicmplx,out,recvcnts,displs,&
         mpicmplx,mp_comm,ierror)
# else
    out=RESHAPE(arr,SHAPE(out))
#endif
  end subroutine allgatherv_complex_array_1to3

  subroutine allgatherv_complex_array_1to1_sub(arr,count,out,recvcnts,displs,sub_comm)
    !<DD>A subroutine to do a allgatherv operation, sending recvcnts(iproc)
    !data from the iproc'th processor to all others starting at arr(start).
    implicit none
    complex, dimension(:), intent(in) :: arr  !The data to gather
    integer, intent(in) :: count !How much data to gather, <=SIZE(arr)
    complex, dimension(:), intent(out) :: out !The gathered data
    integer, dimension(:), intent(in) :: recvcnts !Array detailing how much data to expect from each proc
    integer, dimension(:), intent(in) :: displs !Array detailing offset in array where gathered data is to be stored
    integer, intent(in) :: sub_comm !Sub-communicator handle
# ifdef MPI
    integer :: ierror
    !Do the gather
    call mpi_allgatherv(arr,count,mpicmplx,out,recvcnts,displs,&
         mpicmplx,sub_comm,ierror)
# else
    out=RESHAPE(arr,SHAPE(out))
#endif
  end subroutine allgatherv_complex_array_1to1_sub

  subroutine allgatherv_complex_array_1to3_sub(arr,count,out,recvcnts,displs,sub_comm)
    !<DD>A subroutine to do a allgatherv operation, sending recvcnts(iproc)
    !data from the iproc'th processor to all others starting at arr(start).
    implicit none
    complex, dimension(:), intent(in) :: arr  !The data to gather
    integer, intent(in) :: count !How much data to gather, <=SIZE(arr)
    complex, dimension(:,:,:), intent(out) :: out !The gathered data
    integer, dimension(:), intent(in) :: recvcnts !Array detailing how much data to expect from each proc
    integer, dimension(:), intent(in) :: displs !Array detailing offset in array where gathered data is to be stored
    integer, intent(in) :: sub_comm !Sub-communicator handle
# ifdef MPI
    integer :: ierror
    !Do the gather
    call mpi_allgatherv(arr,count,mpicmplx,out,recvcnts,displs,&
         mpicmplx,sub_comm,ierror)
#else
    out=RESHAPE(arr,SHAPE(out))
#endif
  end subroutine allgatherv_complex_array_1to3_sub

! ************** comm utils *****************************
  subroutine free_comm_id (comm)
    !<DD>A routine to free the communicator with id comm
    implicit none
    integer, intent(inout) :: comm !Communicator id
# ifdef MPI
    integer :: ierror
    call mpi_comm_free(comm,ierror)
# endif
  end subroutine free_comm_id

  subroutine free_comm_type (comm)
    !<DD>A routine to free the communicator represented by comm
    implicit none
    type(comm_type), intent(inout) :: comm !Communicator object
# ifdef MPI
    integer :: ierror
    call mpi_comm_free(comm%id,ierror)
    !Reset variables
    comm%id=-1
    comm%iproc=-1
    comm%nproc=-1
    comm%proc0=.false.
# endif
  end subroutine free_comm_type

!In future we may wish to make split an interface to allow for
!user specific keys (to reorder processor ranks) and to specify
!a different communicator to split
  subroutine split_nokey (col,new_comm)
    !<DD>A routine to split the mp_comm communicator into sub-groups
    !based on each procs specific colour "col". The sub communicator's
    !handle is passed back in new_comm
    implicit none
    integer, intent(inout) :: col !Processors colour
    integer, intent(out) :: new_comm !The new sub communicator's handle
# ifdef MPI
    integer :: ierror
    !Split the comm group, note we use a constant key of 0 across all procs
    !meaning that the rank order is the same in the old and new communicators
    call mpi_comm_split(mp_comm,col,0,new_comm,ierror)
# endif
  end subroutine split_nokey

  subroutine split_key (col,key,new_comm)
    !<DD>A routine to split the mp_comm communicator into sub-groups
    !based on each procs specific colour "col" and ranked by key. The sub communicator's
    !handle is passed back in new_comm
    implicit none
    integer, intent(inout) :: col !Processors colour
    integer, intent(inout) :: key !Processors key, used to determine rank
    integer, intent(out) :: new_comm !The new sub communicator's handle
# ifdef MPI
    integer :: ierror
    !Split the comm group
    call mpi_comm_split(mp_comm,col,key,new_comm,ierror)
# endif
  end subroutine split_key

  subroutine split_nokey_to_commtype (col,new_comm)
    !<DD>A routine to split the mp_comm communicator into sub-groups
    !based on each procs specific colour "col". The sub communicator's
    !handle is passed back in new_comm
    implicit none
    integer, intent(inout) :: col !Processors colour
    type(comm_type), intent(out) :: new_comm !The new sub communicator's handle
# ifdef MPI
    integer :: ierror
    integer :: comm_id, nproc, iproc
    !Split the comm group, note we use a constant key of 0 across all procs
    !meaning that the rank order is the same in the old and new communicators
    call mpi_comm_split(mp_comm,col,0,comm_id,ierror)
    new_comm%id=comm_id
    call nproc_comm(comm_id,nproc)
    new_comm%nproc=nproc
    call rank_comm(comm_id,iproc)
    new_comm%iproc=iproc
    new_comm%proc0=iproc.eq.0
# endif
  end subroutine split_nokey_to_commtype

  subroutine split_key_to_commtype (col,key,new_comm)
    !<DD>A routine to split the mp_comm communicator into sub-groups
    !based on each procs specific colour "col" and ranked by key. The sub communicator's
    !handle is passed back in new_comm
    implicit none
    integer, intent(inout) :: col !Processors colour
    integer, intent(inout) :: key !Processors key, used to determine rank
    type(comm_type), intent(out) :: new_comm !The new sub communicator's handle
# ifdef MPI
    integer :: ierror
    integer :: comm_id, nproc
    !Split the comm group
    call mpi_comm_split(mp_comm,col,key,comm_id,ierror)
    new_comm%id=comm_id
    call nproc_comm(comm_id,nproc)
    new_comm%nproc=nproc
    call rank_comm(comm_id,nproc)
    new_comm%iproc=nproc
    new_comm%proc0=new_comm%iproc.eq.0
# endif
  end subroutine split_key_to_commtype

  subroutine split_nokey_to_commtype_sub (col,new_comm,sub)
    !<DD>A routine to split a subcommunicator into sub-groups
    !based on each procs specific colour "col". The sub communicator's
    !handle is passed back in new_comm
    implicit none
    integer, intent(inout) :: col !Processors colour
    type(comm_type), intent(out) :: new_comm !The new sub communicator's handle
    integer, intent(inout) :: sub
# ifdef MPI
    integer :: ierror
    integer :: comm_id, nproc, iproc
    !Split the comm group, note we use a constant key of 0 across all procs
    !meaning that the rank order is the same in the old and new communicators
    call mpi_comm_split(sub,col,0,comm_id,ierror)
    new_comm%id=comm_id
    call nproc_comm(comm_id,nproc)
    new_comm%nproc=nproc
    call rank_comm(comm_id,iproc)
    new_comm%iproc=iproc
    new_comm%proc0=iproc.eq.0
# endif
  end subroutine split_nokey_to_commtype_sub

  subroutine split_key_to_commtype_sub (col,key,new_comm,sub)
    !<DD>A routine to split a subcommunicator into sub-groups
    !based on each procs specific colour "col" and ranked by key. The sub communicator's
    !handle is passed back in new_comm
    implicit none
    integer, intent(inout) :: col !Processors colour
    integer, intent(inout) :: key !Processors key, used to determine rank
    type(comm_type), intent(out) :: new_comm !The new sub communicator's handle
    integer, intent(inout) :: sub !Subcommunicator to split
# ifdef MPI
    integer :: ierror
    integer :: comm_id, nproc
    !Split the comm group
    call mpi_comm_split(sub,col,key,comm_id,ierror)
    new_comm%id=comm_id
    call nproc_comm(comm_id,nproc)
    new_comm%nproc=nproc
    call rank_comm(comm_id,nproc)
    new_comm%iproc=nproc
    new_comm%proc0=new_comm%iproc.eq.0
# endif
  end subroutine split_key_to_commtype_sub

! ************** broadcasts *****************************

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
# ifdef MPI
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_character

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_integer_array

  subroutine broadcast_integer_2array (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_integer_2array

  subroutine broadcast_real (x)
    implicit none
    real, intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, 1, mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real, dimension (:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_array
  
  subroutine broadcast_real_2array(x)
  implicit none
  real, dimension(:,:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_2array

  subroutine broadcast_real_3array(x)
  implicit none
  real, dimension(:,:,:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_3array

  subroutine broadcast_real_4array(x)
  implicit none
  real, dimension(:,:,:,:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_4array

  subroutine broadcast_complex (z)
    implicit none
    complex, intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, 1, mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_array
  
  subroutine broadcast_complex_2array (z)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_2array

  subroutine broadcast_complex_3array (z)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_3array

  subroutine broadcast_complex_4array (z)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_4array

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_logical_array

  subroutine broadcast_logical_2array (f)
    implicit none
    logical, dimension (:,:), intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_logical_2array

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, mp_comm, ierror)
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
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, mp_comm, ierror)
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
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, mp_comm, ierror)
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
    call mpi_bcast (i, 1, MPI_INTEGER, src, mp_comm, ierror)
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
    call mpi_bcast (i, size(i), MPI_INTEGER, src, mp_comm, ierror)
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
    call mpi_bcast (x, 1, mpireal, src, mp_comm, ierror)
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
    call mpi_bcast (x, size(x), mpireal, src, mp_comm, ierror)
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
    call mpi_bcast (z, 1, mpicmplx, src, mp_comm, ierror)
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
    call mpi_bcast (z, size(z), mpicmplx, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_complex_array

!******************BROADCAST SUB*************************
  subroutine bcastfrom_complex_array_sub (z, src, sub)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    integer, intent (in) :: sub
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, src, sub, ierror)
# else
    if (src /= 0) call error ("broadcast from sub")
# endif
  end subroutine bcastfrom_complex_array_sub

  subroutine bcastfrom_complex_2array_sub (z, src, sub)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    integer, intent (in) :: src
    integer, intent (in) :: sub
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, src, sub, ierror)
# else
    if (src /= 0) call error ("broadcast from sub")
# endif
  end subroutine bcastfrom_complex_2array_sub

  subroutine bcastfrom_complex_3array_sub (z, src, sub)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    integer, intent (in) :: src
    integer, intent (in) :: sub
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, src, sub, ierror)
# else
    if (src /= 0) call error ("broadcast from sub")
# endif
  end subroutine bcastfrom_complex_3array_sub

  subroutine bcastfrom_complex_4array_sub (z, src, sub)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent (in) :: src
    integer, intent (in) :: sub
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, src, sub, ierror)
# else
    if (src /= 0) call error ("broadcast from sub")
# endif
  end subroutine bcastfrom_complex_4array_sub

! ************** reductions ***********************

  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, 1, MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, size(i), MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_integer_array

  subroutine sum_reduce_logical (a, dest)
    implicit none
    logical, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, MPI_LOGICAL, MPI_LOR, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, MPI_LOGICAL, MPI_LOR, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_logical

  subroutine sum_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_real_2array (a, dest)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_2array

  subroutine sum_reduce_real_3array (a, dest)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_3array

  subroutine sum_reduce_real_4array (a, dest)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_4array

  subroutine sum_reduce_real_5array (a, dest)
    implicit none
    real, dimension (:,:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_5array

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, 1, mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, 1, mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_array

  subroutine sum_reduce_complex_2array (z, dest)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_2array

  subroutine sum_reduce_complex_3array (z, dest)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_3array

  subroutine sum_reduce_complex_4array (z, dest)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_4array

  subroutine sum_reduce_complex_5array (z, dest)
    implicit none
    complex, dimension (:,:,:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if (iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_5array

! Sum all z1's into z2 at dest
  subroutine sum_reduce_alt_complex_3array (z1,z2,dest)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z1
    complex, dimension (:,:,:), intent (in out) :: z2
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    call mpi_reduce &
         (z1, z2, size(z1), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_alt_complex_3array

!******************SUM REDUCE SUB*************************
  subroutine sum_reduce_complex_array_sub (z, dest, sub)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    type(comm_type), intent(in) :: sub
# ifdef MPI
    integer :: ierror
    if(sub%iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, sub%id, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, sub%id, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_array_sub

  subroutine sum_reduce_complex_4array_sub (z, dest, sub)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
    type(comm_type), intent(in) :: sub
# ifdef MPI
    integer :: ierror
    if(sub%iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, sub%id, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, sub%id, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_4array_sub

!******************SUM ALLREDUCE*************************  
  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_integer_2array (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_integer_2array

  subroutine sum_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_real_2array (a)
    implicit none
    real, dimension (:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_2array

  subroutine sum_allreduce_real_3array (a)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_3array

  subroutine sum_allreduce_real_4array (a)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_4array

  subroutine sum_allreduce_real_5array (a)
    implicit none
    real, dimension (:,:,:,:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_5array

  subroutine sum_allreduce_complex (z)
    implicit none
    complex, intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, 1, mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_array

  subroutine sum_allreduce_complex_2array (z)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_2array

  subroutine sum_allreduce_complex_3array (z)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_3array

  subroutine sum_allreduce_complex_4array (z)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_4array

  subroutine sum_allreduce_complex_5array (z)
    implicit none
    complex, dimension (:,:,:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_5array
!/Sub-communicator allreduce
  subroutine sum_allreduce_sub_integer (i,sub_comm)
    implicit none
    integer, intent (in out) :: i
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_integer

  subroutine sum_allreduce_sub_integer_array (i,sub_comm)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_integer_array

  subroutine sum_allreduce_sub_real (a,sub_comm)
    implicit none
    real, intent (in out) :: a
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_real

  subroutine sum_allreduce_sub_real_array (a,sub_comm)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_real_array

  subroutine sum_allreduce_sub_real_2array (a,sub_comm)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_real_2array

  subroutine sum_allreduce_sub_real_3array (a,sub_comm)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_real_3array

  subroutine sum_allreduce_sub_real_4array (a,sub_comm)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_real_4array

  subroutine sum_allreduce_sub_real_5array (a,sub_comm)
    implicit none
    real, dimension (:,:,:,:,:), intent (in out) :: a
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_real_5array

  subroutine sum_allreduce_sub_complex (z,sub_comm)
    implicit none
    complex, intent (in out) :: z
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, 1, mpicmplx, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_complex

  subroutine sum_allreduce_sub_complex_array (z,sub_comm)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_complex_array

  subroutine sum_allreduce_sub_complex_2array (z,sub_comm)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_complex_2array

  subroutine sum_allreduce_sub_complex_3array (z,sub_comm)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_complex_3array

  subroutine sum_allreduce_sub_complex_4array (z,sub_comm)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_complex_4array

  subroutine sum_allreduce_sub_complex_5array (z,sub_comm)
    implicit none
    complex, dimension (:,:,:,:,:), intent (in out) :: z
    integer, intent(in) :: sub_comm
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, sub_comm, ierror)
# endif
  end subroutine sum_allreduce_sub_complex_5array

  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, 1, MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, size(i), MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_integer_array

  subroutine max_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, mpireal, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, mpireal, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, size(a), mpireal, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, size(a), mpireal, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_real_array

  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_real_array

  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, 1, MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_integer

  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, size(i), MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_integer_array

  subroutine min_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, mpireal, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, mpireal, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, size(a), mpireal, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, size(a), mpireal, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_real_array

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_real_array

! ********************* barrier **********************

  subroutine barrier_nocomm
# ifdef MPI
    implicit none
    integer :: ierror
    call mpi_barrier (mp_comm, ierror)
# endif
  end subroutine barrier_nocomm

  subroutine barrier_comm(comm)
    implicit none
    integer,intent(in)::comm
# ifdef MPI
    integer :: ierror
    call mpi_barrier (comm, ierror)
# endif
  end subroutine barrier_comm

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
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, mp_comm, ierror)
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
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, mp_comm, ierror)
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
    call mpi_send (a, 1, mpireal, dest, tagp, mp_comm, ierror)
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
    call mpi_send (a, size(a), mpireal, dest, tagp, mp_comm, ierror)
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
    call mpi_send (z, 1, mpicmplx, dest, tagp, mp_comm, ierror)
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
    call mpi_send (z, size(z), mpicmplx, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_complex_array

  subroutine nonblocking_send_complex_array (z, dest, tag, request)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer, intent (out) :: request
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, size(z), mpicmplx, dest, tagp, mp_comm, request, ierror)
# else
    call error ("send")
# endif
  end subroutine nonblocking_send_complex_array

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
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, mp_comm, ierror)
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
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, mp_comm, ierror)
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
         (s, len(s), MPI_CHARACTER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_character

! MAB> needed for Trinity
! ********************* synchronous sends **********************

  subroutine ssend_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (i, 1, MPI_INTEGER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_integer

  subroutine ssend_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (i, size(i), MPI_INTEGER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_integer_array

  subroutine ssend_real (a, dest, tag)
    implicit none
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_real

  subroutine ssend_real_array (a, dest, tag)
    implicit none
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_real_array

  subroutine ssend_complex (z, dest, tag)
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_complex

  subroutine ssend_complex_array (z, dest, tag)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_complex_array
  
  subroutine ssend_complex_2array (z, dest, tag)
    implicit none
    complex, dimension (:,:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_complex_2array
  
  subroutine ssend_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (f, 1, MPI_LOGICAL, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_logical

  subroutine ssend_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (f, size(f), MPI_LOGICAL, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_logical_array

  subroutine ssend_character (s, dest, tag)
    implicit none
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend &
         (s, len(s), MPI_CHARACTER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_character
! <MAB

! ********************* non-blocking sends **********************
  subroutine nbsend_real_array(z,dest,tag,handle)
    !<DD>Routine for nonblocking send of z to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    real, dimension(:), intent(in) :: z
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_isend(z,size(z),mpireal,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbsend_real_array

  subroutine nbsend_real_array_count(z,count,dest,tag,handle)
    !<DD>Routine for nonblocking send of z (size=count) to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    real, dimension(:), intent(in) :: z
    integer, intent(in) :: count
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_isend(z,count,mpireal,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbsend_real_array_count

  subroutine nbsend_complex_array(z,dest,tag,handle)
    !<DD>Routine for nonblocking send of z to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    complex, dimension(:), intent(in) :: z
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_isend(z,size(z),mpicmplx,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbsend_complex_array

  subroutine nbsend_complex_array_sub(z,dest,tag,sub,handle)
    !<DD>Routine for nonblocking send of z to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    complex, dimension(:), intent(in) :: z
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    type(comm_type), intent(in) :: sub
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_isend(z,size(z),mpicmplx,dest,tag,sub%id,handle,ierror)
# endif
  end subroutine nbsend_complex_array_sub

  subroutine nbsend_complex_array_count(z,count,dest,tag,handle)
    !<DD>Routine for nonblocking send of z (size=count) to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    complex, dimension(:), intent(in) :: z
    integer, intent(in) :: count
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_isend(z,count,mpicmplx,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbsend_complex_array_count

! ********************* persistent sends **********************
  
  subroutine send_init_complex_array(z,dest,tag,handle)
    !Routine to initialise a persistent send operations
    implicit none
    complex, dimension(:), intent(in) :: z
    integer, intent(in) :: dest, tag
    integer, intent(inout) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_send_init(z,size(z),mpicmplx,dest,tag,mp_comm,handle,ierror)
#endif
  end subroutine send_init_complex_array

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
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, mp_comm, &
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
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, mp_comm, &
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
    call mpi_recv (a, 1, mpireal, src, tagp, mp_comm, &
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
    call mpi_recv (a, size(a), mpireal, src, tagp, mp_comm, &
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
    call mpi_recv (z, 1, mpicmplx, src, tagp, mp_comm, &
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
    call mpi_recv (z, size(z), mpicmplx, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex_array

  subroutine receive_complex_2array (z, src, tag)
    implicit none
# ifdef MPI
    complex, dimension (:,:), intent (out) :: z
# else
    complex, dimension (:,:) :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), mpicmplx, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex_2array
  
  subroutine nonblocking_receive_complex_array (z, src, tag, request)
    implicit none
# ifdef MPI
    complex, dimension (:), intent (inout) :: z
# else
    complex, dimension (:) :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer, intent (out) :: request
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, size(z), mpicmplx, src, tagp, mp_comm, &
        request, ierror)
# else
    call error ("receive")
# endif
  end subroutine nonblocking_receive_complex_array

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
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, mp_comm, &
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
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, mp_comm, &
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
    call mpi_recv (s, len(s), MPI_CHARACTER, src, tagp, mp_comm, &
         status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_character

! ********************* non-blocking receives **********************
  subroutine nbrecv_real_array(z,dest,tag,handle)
    !<DD>Routine for nonblocking recv of z to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    real, dimension(:), intent(out) :: z
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_irecv(z,size(z),mpireal,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbrecv_real_array

  subroutine nbrecv_real_array_count(z,count,dest,tag,handle)
    !<DD>Routine for nonblocking recv of z (size=count) to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    real, dimension(:), intent(out) :: z
    integer, intent(in) :: count
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_irecv(z,count,mpireal,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbrecv_real_array_count

  subroutine nbrecv_complex_array(z,dest,tag,handle)
    !<DD>Routine for nonblocking recv of z to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    complex, dimension(:), intent(out) :: z
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_irecv(z,size(z),mpicmplx,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbrecv_complex_array

  subroutine nbrecv_complex_array_sub(z,dest,tag,sub,handle)
    !<DD>Routine for nonblocking recv of z to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    complex, dimension(:), intent(out) :: z
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    type(comm_type), intent(in) :: sub
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_irecv(z,size(z),mpicmplx,dest,tag,sub%id,handle,ierror)
# endif
  end subroutine nbrecv_complex_array_sub

  subroutine nbrecv_complex_array_count(z,count,dest,tag,handle)
    !<DD>Routine for nonblocking recv of z (size=count) to dest. Use
    !tag to label message and return handle for later checking.
    implicit none
    complex, dimension(:), intent(out) :: z
    integer, intent(in) :: count
    integer, intent(in) :: dest
    integer, intent(in) :: tag
    integer,intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_irecv(z,count,mpicmplx,dest,tag,mp_comm,handle,ierror)
# endif
  end subroutine nbrecv_complex_array_count

! ********************* persistent recvs **********************
  
  subroutine recv_init_complex_array(z,dest,tag,handle)
    !Routine to initialise a persistent recv operations
    implicit none
    complex, dimension(:), intent(in) :: z
    integer, intent(in) :: dest, tag
    integer, intent(out) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_recv_init(z,size(z),mpicmplx,dest,tag,mp_comm,handle,ierror)
#endif
  end subroutine recv_init_complex_array

! ********************* non-blocking checks **********************
  subroutine waitall_nostat (count, requests)
    !<DD>A routine to wait for all count communications, given by the
    !message handles in requests, to complete. We ignore the status
    !information
    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent (inout) :: requests
    !Note mpi_wait will set the request handle to MPI_NULL (or similar)
    !when message complete so have to set requests as inout
# ifdef MPI
    integer :: ierror
    call mpi_waitall(count,requests,MPI_STATUSES_IGNORE,ierror)
# endif
  end subroutine waitall_nostat

  subroutine waitall_stat (count, requests, status)
    !<DD>A routine to wait for all count communications, given by the
    !message handles in requests, to complete. We return the message
    !statuses
    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent (inout) :: requests
# ifdef MPI
    integer, dimension(MPI_STATUS_SIZE,count), intent(out) :: status
# else
    integer,intent(out) :: status
# endif
# ifdef MPI
    !Note mpi_wait will set the request handle to MPI_NULL (or similar)
    !when message complete so have to set requests as inout
    integer :: ierror
    call mpi_waitall(count,requests,status,ierror)
# endif
  end subroutine waitall_stat

  subroutine waitany (count, requests, requestindex, status)

    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent(inout) :: requests
    integer, intent(out) :: requestindex
# ifdef MPI
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
# else
    integer, dimension(1), intent(out) :: status
# endif

# ifdef MPI
    integer :: ierror

    call mpi_waitany(count, requests, requestindex, status, ierror)
# else
    call error ("waitany")
# endif

  end subroutine waitany

  subroutine testall_nostat_noflag (count, requests)
    !<DD>A routine to test for all count communications, given by the
    !message handles in requests, to complete. We ignore the status
    !and flag information
    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent (inout) :: requests
    !Note mpi_test will set the request handle to MPI_NULL (or similar)
    !when message complete so have to set requests as inout
# ifdef MPI
    logical :: flag !Should this be integer?
    integer :: ierror
    call mpi_testall(count,requests,flag,MPI_STATUSES_IGNORE,ierror)
# endif
  end subroutine testall_nostat_noflag

  subroutine testall_nostat (count, requests, flag)
    !<DD>A routine to test for all count communications, given by the
    !message handles in requests, to complete. We ignore the status
    !information
    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent (inout) :: requests
    logical, intent(out) :: flag
    !Note mpi_test will set the request handle to MPI_NULL (or similar)
    !when message complete so have to set requests as inout
# ifdef MPI
    integer :: ierror
    call mpi_testall(count,requests,flag,MPI_STATUSES_IGNORE,ierror)
# endif
  end subroutine testall_nostat

  subroutine testall_stat (count, requests, status, flag)
    !<DD>A routine to test for all count communications, given by the
    !message handles in requests, to complete. We return the message
    !statuses and completion flag
    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent (inout) :: requests
# ifdef MPI
    integer, dimension(MPI_STATUS_SIZE,count), intent(out) :: status
# else
    integer,intent(out) :: status
# endif
    logical, intent(out) :: flag
# ifdef MPI
    !Note mpi_test will set the request handle to MPI_NULL (or similar)
    !when message complete so have to set requests as inout
    integer :: ierror
    call mpi_testall(count,requests,flag,status,ierror)
# endif
  end subroutine testall_stat

! ********************* persistent utilities **********************
  subroutine start_persist(handle)
    !Starts a single persistent communication represented by handle
    integer, intent(inout) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_start(handle,ierror)
# endif
  end subroutine start_persist

  subroutine startall_persist(handles)
    !Starts multiple persistent communications represented by handles
    integer, dimension(:), intent(inout) :: handles
# ifdef MPI
    integer :: ierror
    call mpi_startall(size(handles),handles,ierror)
# endif
  end subroutine startall_persist

  subroutine free_handle_persist(handle)
    !Frees a single persistent request handle
    integer, intent(inout) :: handle
# ifdef MPI
    integer :: ierror
    call mpi_request_free(handle,ierror)
# endif
  end subroutine free_handle_persist

  subroutine free_handles_persist(handles)
    !Frees multiple persistent request handles
    integer, dimension(:), intent(inout) :: handles
# ifdef MPI
    integer :: i
    do i=1,size(handles)
       call free_handle_persist(handles(i))
    enddo
# endif
  end subroutine free_handles_persist

  
! ********************* other routines **********************

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
    integer, dimension(0:), intent (out) :: group0
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

    if (.not. allocated(grp0)) allocate (grp0(0:size(group0)-1))
    
! calculate dimensions  mpi processor grid will have and check that 
! ncolumns*nrows = number of processes

! nrows is # of processors per job (or group)    
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

!    call mpi_cart_create(mpi_comm_world, ndim, dims, period, reorder, comm2d, ierr)
    call mpi_cart_create(comm_all, ndim, dims, period, reorder, comm2d, ierr)
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
     
! MAB> following two lines were incorrect
!    j=1
!    group0(0) = 0
! replace with
    j = 0
    if (proc0 .and. gproc0) then
       group0(0) = 0
       j = 1
    end if
! <MAB
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

    grp0 = group0    

! TT> brought down here from init_job_name in file_utils.fpp
    call scope (subprocs)
! <TT

# endif
  end subroutine init_jobs

  subroutine all_to_group_real (all, group, njobs)
    
    implicit none

    real, dimension (:), intent (in) :: all
    real, intent (out) :: group
    integer, intent (in) :: njobs

    integer :: ik, tag

    ! TESTING -- MAB
    integer :: idx

    tag = 1000

# ifndef MPI
    call error("all_to_group")
# else
!    do ik = 0, njobs-1
!       if (proc0) then
!          if (iproc == grp0(ik)) then
!             group = all(ik+1)
!          else
!             call ssend (all(ik+1), grp0(ik), tag)
!          end if
!       else if (iproc == grp0(ik)) then
!          call receive (group, 0, tag)
!       end if
!!       call barrier
!    end do
    do ik = 0, njobs-1
       if (proc0) then
          idx = mod(ik,size(all))
          if (iproc == grp0(ik)) then
             group = all(idx+1)
          else
             call ssend (all(idx+1), grp0(ik), tag)
          end if
       else if (iproc == grp0(ik)) then
          call receive (group, 0, tag)
       end if
!       call barrier
    end do

# endif
  end subroutine all_to_group_real

  subroutine all_to_group_real_array (all, group, njobs)
    
    implicit none

    real, dimension (:,:), intent (in) :: all
    real, dimension (:), intent (out) :: group
    integer, intent (in) :: njobs

    integer :: ik, tag

    ! TESTING -- MAB
    integer :: idx

# ifndef MPI
    call error ("all_to_group")
# else
    tag = 1001

!    do ik = 0, njobs-1
!       if (proc0) then
!          if (iproc == grp0(ik)) then
!             group = all(ik+1,:)
!          else
!             call ssend (all(ik+1,:), grp0(ik), tag)
!          end if
!       else if (iproc == grp0(ik)) then
!          call receive (group, 0, tag)
!       end if
!!       call barrier
!    end do
    do ik = 0, njobs-1
       if (proc0) then
          idx = mod(ik,size(all,dim=1))
          if (iproc == grp0(ik)) then
             group = all(idx+1,:)
          else
             call ssend (all(idx+1,:), grp0(ik), tag)
          end if
       else if (iproc == grp0(ik)) then
          call receive (group, 0, tag)
       end if
!       call barrier
    end do

# endif
  end subroutine all_to_group_real_array

  subroutine group_to_all_real (group, all, njobs)
    
    implicit none

    real, intent (in) :: group
    real, dimension (:), intent (out) :: all
    integer, intent (in) :: njobs

    integer :: ik, tag

    integer :: idx

    tag = 1002

# ifndef MPI
    call error("group_to_all")
# else
    do ik = 0, njobs-1
       if (iproc == grp0(ik)) then
          if (.not. proc0) then
             call ssend (group, 0, tag)
          else
             idx = mod(ik,size(all))
             all(idx+1) = group
          end if
       else if (proc0) then
          idx = mod(ik,size(all))
          call receive (all(idx+1), grp0(ik), tag)
       end if
!       call barrier
    end do

# endif    
  end subroutine group_to_all_real

  subroutine group_to_all_real_array (group, all, njobs)
    
    implicit none

    real, dimension (:), intent (in) :: group
    real, dimension (:,:), intent (out) :: all
    integer, intent (in) :: njobs

    integer :: ik, tag

    integer :: idx

    tag = 1003

# ifndef MPI
    call error("group_to_all")
# else
    do ik = 0, njobs-1
       if (iproc == grp0(ik)) then
          if (.not. proc0) then
             call ssend (group, 0, tag)
          else
             idx = mod(ik,size(all))
             all(idx+1,:) = group
          end if
       else if (proc0) then
          idx = mod(ik,size(all))
          call receive (all(idx+1,:), grp0(ik), tag)
       end if
!       call barrier
    end do
    
# endif
  end subroutine group_to_all_real_array

  subroutine mp_abort (msg,to_screen)
    use file_utils, only: error_unit, flush_output_file
    implicit none
    character(len=*), intent (in) :: msg
    logical, intent(in), optional :: to_screen
# ifdef MPI
    integer :: ierror
    integer, parameter :: error_code=MPI_ERR_UNKNOWN
# endif

    if (proc0) then
       if(present(to_screen))then
          if(to_screen) write(6,*) "Error:"//msg
       endif
       write (error_unit(),*) "Error: "//msg
       call flush_output_file (error_unit())
    end if

# ifndef MPI
# ifndef NO_ABORT
    call abort
# endif
# else
    call mpi_abort(comm_all, error_code, ierror)
# endif
  end subroutine mp_abort


# ifndef MPI
  subroutine error (msg)
    use file_utils, only: error_unit
    implicit none
    character(len=*), intent (in) :: msg

    write (error_unit(),*) "mp error: "//msg
# ifndef NO_ABORT
    call abort
# endif
    stop
  end subroutine error
# endif




end module mp
