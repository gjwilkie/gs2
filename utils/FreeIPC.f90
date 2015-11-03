!> \mainpage
!! FreeIPC - A module to allow Fortran programs that use MPI
!!           to easily access shared memory within multicore nodes
!! \date 28th Nov 2009
!! \version 0.0 
!! \author  &copy; Ian Bush, Nitya Hariharan
!!
!! \section Legal Legal
!!
!!    This library is free software; you can redistribute it and/or
!!    modify it under the terms of the GNU Lesser General Public
!!    License as published by the Free Software Foundation; either
!!    version 3.0 of the License, or (at your option) any later version.
!!
!!    This library is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!    Lesser General Public License for more details.
!!
!!    You should have received a copy of the GNU Lesser General Public
!!    License along with this library; if not, write to the Free Software
!!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!!
!!
!! This library is released under LGPL. The full text of LGPL may be found
!! at http://www.gnu.org/licenses/lgpl-3.0.txt
!!
!! In particular, the terms of the GNU LGPL imply that you CAN use this library with
!! proprietary code with some restrictions (you can link, but not use the code directly). 
!! Please consult the above URL for exact details.
!!
!! \section Introduction Introduction
!!
!! FreeIPC is a Fortran module that attempts to allow MPI programs easy access to
!! stard system V IPC facilities and related concepts. At present the following facilities
!! are supported: 
!! - Creation and deletion of shared memory segments
!! - Critical regions
!! - intra-node synchronisation
!!
!! See
!!
!! http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_07.html#tag_02_07
!!
!! for more details.
!!
!! \section Compilation Compilation
!!
!! To module is written in standard fortran 95, except for using the the following fortran 2003
!! features
!! - Interoperability with C
!! - Intent for Pointer dummy arguments
!! - Allocatable components of derived types ( i.e. TR15581 )
!!
!! It has been shown to compile succesfully with the following compilers:
!! - gfortran ( requires version 4.3 or later )
!! - g95
!! - Intel
!! - Portland Group
!! - Pathscale
!! - Cray
!! - Sun
!! - IBM
!!
!! MPI must also be available
!!
!! \section Use Use
!!
!! As the code is a Fortran module to access the provided facilities simply
!!
!! \c Use \c FIPC_module
!!
!! is all that is required. FreeIPC reserves all symbols that start with the
!! characters \c FIPC_ for its own use, and programs that use FreeIPC should avoid
!! using such symbols to avoid name clashes.
!!
!! Before any FreeIPC facility is used MPI must be initialised. FreeIPC is itself
!! initialized by FIPC_init, and finalized with FIPC_finalize.
!!
!! The interface is deliberately very similar to MPI. In FreeIPC the place of the MPI
!! communicator is taken by a \e context.  
!! 
!! \subsection Contexts Contexts
!!
!! A context can be viewed as an extension to an MPI communicator
!! as it contains not only a set of processes but also information about how these map
!! onto any shared memory hardware. In fact it contains 3 communicators
!! - The world communicator: This is simply all processes spanned by the context
!! - The intra-communicator: Processes in the same intra-communicator may communicate with
!!                           each other through shared memory facilities. Therefore all proceeses
!!                           in an intra-communicator may all access, for example, the same shared
!!                           shared memory segment
!! - The extra-communicator: This contains all processes which have rank 0 in an intra communicator.
!!                           It is useful for communicating shared memory segments between nodes
!!
!! Programs may create their own contexts, and any shared memory features of the architecture
!! are automatically detected. Note that an intra-communicator need not span the whole of a shared
!! memory node, and it is possible to have many different contexts on a shared memory node each
!! spanning differnet subsets of the processors on the node.
!!
!! Unlike MPI a context is an opaque entity of derived type \c FIPC_ctxt. Once FreeIPC is initialised
!! the context \c FIPC_ctxt_world is available for use. This is the equivalent of MPI_comm_world; it spans
!! all processors over which FreeIPC is initialised and the intra-communicators spann all processes on the
!! ( autmoaticaly detected ) shared memory node.
!!
!! \subsection Seg_creat Shared Memory Segment Creation And Deletion
!!
!! Shared memory segments are created through the interface FIP_seg_create, and freed by
!! FIPC_seg_free. Arrays of 3 data types are supported
!! - Integer of kind c_int i.e. default integer
!! - Real of kind c_double i.e. "Double Precision"
!! - Complex of kind c_double_complex i.e. "Double Complex"
!! The arrays may be 1 to 7 dimensional. Segment creation occurs within a context,
!! and one segment is created per intra-communicator.
!! So, for example, if the context supplied is FIPC_ctxt_world, one segment is created
!! on each of the ( shared memory ) nodes of the machine.
!!
!! FIPC_seg_create returns a fortran pointer with the appropriate attributes for the
!! data type and kind required. Thus the program now views the segment as a normal
!! Fortran datatype, and the full power of the language is available to process data in it!
!!
!! \subsection Synchronization Synchronization
!!
!! As with MPI and OpenMP it is the duty of the application programmer to ensure that
!! the processes in the job are properly synchronized so as to avoid race conditions.
!! FreeIPC provides two synchronization mechanisms over and above those provided by MPI:
!! - Critical regions
!! - Intra-node barriers
!!
!! \subsection Example An Example Program
!!
!! \verbatim
!! Program FreeIPC_example
!! 
!!   Use MPI
!!   Use FIPC_module
!! 
!!   Implicit None
!! 
!!   Integer, Dimension( : ), Pointer :: seg
!! 
!!   Integer :: rank
!!   Integer :: error
!!   Integer :: i
!! 
!!   Call MPI_init( error )
!!   Call MPI_comm_rank( MPI_comm_world, rank, error )
!! 
!!   ! Start up FreeIPC
!!   Call FIPC_init( MPI_comm_world, error )
!! 
!!   ! Create a segment that consists of one integer
!!   Call FIPC_seg_create( FIPC_ctxt_world, (/ 1 /), seg, error )
!! 
!!   seg = 0
!!   ! Make sure everybody has initialized seg
!!   Call FIPC_node_barrier( FIPC_ctxt_world, error )
!! 
!!   ! A deliberate race condition !
!!   Do i = 1, 1000
!!      seg = seg + 1
!!   End Do
!!   ! Make sure everybody has finished the loop
!!   Call FIPC_node_barrier( FIPC_ctxt_world, error )
!!   If( rank == 0 ) Then
!!      Write( *, * ) seg
!!   End If
!!   ! Make sure the result has been written out
!!   Call FIPC_node_barrier( FIPC_ctxt_world, error )
!! 
!!   ! Now do it properly
!!   seg = 0
!!   ! Make sure everybody has initialized seg
!!   Call FIPC_node_barrier( FIPC_ctxt_world, error )
!!   Do i = 1, 1000
!!      Call FIPC_critical_start( FIPC_ctxt_world, error )
!!      seg = seg + 1
!!      Call FIPC_critical_end( FIPC_ctxt_world, error )
!!   End Do
!!   ! Make sure everybody has finished the loop
!!   Call FIPC_node_barrier( FIPC_ctxt_world, error )
!!   If( rank == 0 ) Then
!!      Write( *, * ) seg
!!   End If
!!   ! Make sure the result has been written out
!!   Call FIPC_node_barrier( FIPC_ctxt_world, error )
!! 
!!   ! Free our segment
!!   Call FIPC_seg_free( seg, error )
!!
!!   ! Close down FreeIPC
!!   Call FIPC_finalize( error )
!! 
!!   Call MPI_finalize( error )
!! 
!! End Program FreeIPC_example
!! \endverbatim
!!
!! Example output:
!! \verbatim
!! Wot now ? mpirun -np 2  ./a.out
!!        1073
!!        2000
!! \endverbatim


Module FIPC_module

  ! A library to allow Fortran programs that use MPI to easily access 
  ! shared memory within multicore nodes

  ! It uses standard Fortran 2003 and a few thin C wrappers. Fortran
  ! 2003 C interoperability is used extensively, however the rest of
  ! the code is standard Fortran 95, expcept for the use of allocatble
  ! components of derived types - again a standard Fortran 2003 language 
  ! feature.

  ! The shared memory facilities are provided by use of the standard System V
  ! IPC facilities. See 
  !
  ! http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_07.html#tag_02_07
  !
  ! for details

  Use, Intrinsic :: iso_c_binding, Only : c_int, c_long, c_double, c_ptr, &
       c_f_pointer, c_null_ptr, c_associated, c_loc

!!  Use mpi

  Implicit None

  Include 'mpif.h'

  !
  ! Little type to handle communicators and associated data
  !
  !!!> \cond
  Type, Private :: communicator
     Private
     !> \private
     Logical :: initialized = .False.
     !> \private
     Integer :: handle
     !> \private
     Integer :: size
     !> \private
     Integer :: rank
  End Type communicator
!!!> \endcond

  !>
  !! \brief An opaque type that represents a FIPC context.
  !!
  !! This type represents a FIPC context. It is an opaque data type with no public components
  !! All operations occur within a context. It is very similar to
  !! to a mpi communicator with a bit of extra stuff held to look
  !! after the shared memory parts of the architecture. 
  !!
  !! \if for_fipc_implementors
  !! It consists of a:
  !! INITIALIZED: Is this context set up ?
  !!              Note we only use the allocation status of the pointer to check this
  !! WORLD_COMM : The communicator spanning all processes in this context
  !! INTRA_COMM : A communicator spanning all process  in the context on this shared memory node
  !! EXTRA_COMM : A communicator spanning all the process zeros in all the intra_comms
  !!              on all the nodes which hold processes in WORLD_COMM
  !! SEMID      : The handle of a semaphore shared by members of the INTRA_COMM
  !> \endif
  !
  Type, Public :: FIPC_ctxt
     Private
     !> \private
     Logical             , Pointer :: initialized => Null() 
     !> \private
     Type( communicator ), Pointer :: world_comm
     !> \private
     Type( communicator ), Pointer :: intra_comm
     !> \private
     Type( communicator ), Pointer :: extra_comm
     !> \private
     Integer (c_int)               :: semid
  End Type FIPC_ctxt

  !> The "default" context. Compare MPI_COMM_WORLD.
  Type( FIPC_ctxt ), Save, Public :: FIPC_ctxt_world

  ! Error flags. Note we avoid clashes with the mpi error flags.

  !> Return value to indicate succesful completion
  Integer, Parameter, Public :: FIPC_success               =  0 

  !> Return value indicating an attempt was made to initialise FreeIPC when it was already initialised
  Integer, Parameter, Public :: FIPC_already_initialized   =  1 + MPI_ERR_LASTCODE

  !> Return value indicating a memory allocation failed within FreeIPC
  Integer, Parameter, Public :: FIPC_allocation_failed     =  2 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to get a shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_get_failed        =  3 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to attach to a shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_attach_failed     =  4 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to inquire the properties shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_inquire_failed    =  5 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to detach from shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_detach_failed     =  6 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to remove a shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_remove_failed     =  7 + MPI_ERR_LASTCODE

  !> Return value indicating that a FreeIPC internal consistency check failed
  Integer, Parameter, Public :: FIPC_sanity_failed         =  8 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt was made to free FIPC_ctxt_world
  Integer, Parameter, Public :: FIPC_freeing_ctxt_world    =  9 + MPI_ERR_LASTCODE

  !> Return value indicating that FreeIPC could not identify which shared memory segment is to be freed
  Integer, Parameter, Public :: FIPC_seg_not_found         = 10 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt was made to free a NULL pointer
  Integer, Parameter, Public :: FIPC_free_null_pointer     = 11 + MPI_ERR_LASTCODE

  !> Return value indicating that an insufficent dimensions were supplied when allocating a shared memory segment
  Integer, Parameter, Public :: FIPC_insufficient_dims     = 12 + MPI_ERR_LASTCODE

  !> Return value indicating that FreeIPC was not initalized
  Integer, Parameter, Public :: FIPC_not_initialized       = 13 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt was made to free a context that was still in use, for instance
  !> a shared memory segment still exists within that context
  Integer, Parameter, Public :: FIPC_ctxt_in_use           = 14 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to get a semaphore failed
  Integer, Parameter, Public :: FIPC_sem_get_failed        = 15 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to remove a semaphore failed
  Integer, Parameter, Public :: FIPC_sem_remove_failed     = 16 + MPI_ERR_LASTCODE

  !> Return value indicating that FreeIPC could not identify which semaphore is to be freed
  Integer, Parameter, Public :: FIPC_sem_not_found         = 17 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to start a critical region failed
  Integer, Parameter, Public :: FIPC_critical_start_failed = 18 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to start a critical region failed
  Integer, Parameter, Public :: FIPC_critical_end_failed   = 19 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to initialise posix shared memory data structure failed
  Integer, Parameter, Public :: FIPC_node_init_failed    = 20 + MPI_ERR_LASTCODE

  !> Return value indicating that the node posix data were not initialised
  Integer, Parameter, Public :: FIPC_node_uninitialized  = 21 + MPI_ERR_LASTCODE

  !> Value to indicate that this process should not be in the context that results
  !> from a FIPC_ctxt_split
  Integer, Parameter, Public :: FIPC_undefined = MPI_undefined

  !> Value to indicate that an attempt has been made to extract a non-existant communicator
  !> from a context. The most common case is trying to extract the extra communicator on a process
  !> that is not rank zero in an intra communicator
  Integer, Parameter, Public :: FIPC_comm_null = MPI_comm_null

  ! Reductions available
  !> Handle to indicate a global sum is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_sum  = mpi_sum
  !> Handle to indicate a global product is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_prod = mpi_prod
  !> Handle to indicate a global minimum is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_max  = mpi_max
  !> Handle to indicate a global maximum is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_min  = mpi_min

!Added by Nitya
  !>File to write out the segment and semaphores id to.
  Character(len=50), Parameter, Public :: id_file = 'sms_id.txt'
  !>Type for id, set to 1 for Shared Memory Segment, 2 for Semaphore
  Integer :: id_type  

  ! Public interfaces
  Public :: FIPC_node_init
  Public :: FIPC_ctxt_init
  Public :: FIPC_init
  Public :: FIPC_initialized
  Public :: FIPC_finalize
  Public :: FIPC_ctxt_finalize
  Public :: FIPC_node_finalize
  Public :: FIPC_ctxt_dup
  Public :: FIPC_ctxt_split
  Public :: FIPC_ctxt_free
  Public :: FIPC_ctxt_intra_comm
  Public :: FIPC_ctxt_world_comm
  Public :: FIPC_ctxt_extra_comm
  Public :: FIPC_seg_create
  Public :: FIPC_seg_free
  Public :: FIPC_node_barrier
  Public :: FIPC_critical_start
  Public :: FIPC_critical_end

!Added by Nitya, allreduce was uncommented
  Public :: FIPC_allreduce
  Private :: FIPC_write_sms

  Private

  !> \brief Initialise data structure needed for Posix flavor
  Interface FIPC_node_init
      Module Procedure FIPC_node_init_world
      Module Procedure FIPC_node_init_comm
  End Interface

  !> \brief Create a shared memory segment.
  !!
  !! Sets up a shared memory segment on each of the shared memory nodes spanned by the context ctxt. The interface is\n
  !! Subroutine FIPC_seg_create( 
  !!                             Type( FIPC_ctxt ),intent(in)    ctxt,\n 
  !!                             Integer, Dimension( : ),intent(in) n,\n 
  !!                             <choice>,Dimension( <1-7d> ),intent(out) a,\n
  !!                             Integer,intent(out) error )
  !! \param ctxt  The context within which to create the segment    
  !! \param n     An array containg the dimensions of the arry which will be stored in the segment
  !! \param a     A pointer to a 1-7 dimensionsal array of type Integer( c_int ), Real( c_double ) or
  !!              Complex( complex ). On exit this points to the shared memory segment
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Interface FIPC_seg_create
     Module Procedure seg_create_integer_1d_size_in_int
     Module Procedure seg_create_integer_2d_size_in_int
     Module Procedure seg_create_integer_3d_size_in_int
     Module Procedure seg_create_integer_4d_size_in_int
     Module Procedure seg_create_integer_5d_size_in_int
     Module Procedure seg_create_integer_6d_size_in_int
     Module Procedure seg_create_integer_7d_size_in_int
     Module Procedure seg_create_double_1d_size_in_int
     Module Procedure seg_create_double_2d_size_in_int
     Module Procedure seg_create_double_3d_size_in_int
     Module Procedure seg_create_double_4d_size_in_int
     Module Procedure seg_create_double_5d_size_in_int
     Module Procedure seg_create_double_6d_size_in_int
     Module Procedure seg_create_double_7d_size_in_int
     Module Procedure seg_create_complex_1d_size_in_int
     Module Procedure seg_create_complex_2d_size_in_int
     Module Procedure seg_create_complex_3d_size_in_int
     Module Procedure seg_create_complex_4d_size_in_int
     Module Procedure seg_create_complex_5d_size_in_int
     Module Procedure seg_create_complex_6d_size_in_int
     Module Procedure seg_create_complex_7d_size_in_int
  End Interface

  !> 
  !!
  !! \brief Frees a shared memory segment. 
  !!
  !! This routine frees a shared memory segment on each of the shared memory nodes spanned by the context ctxt
  !>
  Interface FIPC_seg_free
     Module Procedure seg_free_integer_1d
     Module Procedure seg_free_integer_2d
     Module Procedure seg_free_integer_3d
     Module Procedure seg_free_integer_4d
     Module Procedure seg_free_integer_5d
     Module Procedure seg_free_integer_6d
     Module Procedure seg_free_integer_7d
     Module Procedure seg_free_double_1d
     Module Procedure seg_free_double_2d
     Module Procedure seg_free_double_3d
     Module Procedure seg_free_double_4d
     Module Procedure seg_free_double_5d
     Module Procedure seg_free_double_6d
     Module Procedure seg_free_double_7d
     Module Procedure seg_free_complex_1d
     Module Procedure seg_free_complex_2d
     Module Procedure seg_free_complex_3d
     Module Procedure seg_free_complex_4d
     Module Procedure seg_free_complex_5d
     Module Procedure seg_free_complex_6d
     Module Procedure seg_free_complex_7d
  End Interface

  !!!!> \cond for_fipc_implementors

  ! Error flags that can be returned by the System V routines
  ! Need to read these from the C header files
  Integer, Private :: EACCES
  Integer, Private :: EEXIST
  Integer, Private :: EINVAL
  Integer, Private :: ENFILE
  Integer, Private :: ENOENT
  Integer, Private :: ENOMEM
  Integer, Private :: ENOSPC

  ! Flags for control of the creation of shared beasties
  Integer( c_int ), Private :: SEG_NOCREATE   = 0 
  Integer( c_int ), Private :: SEG_CREATE     = 1
  Integer( c_int ), Private :: SEG_NOEXCLUDE  = 0
  Integer( c_int ), Private :: SEG_EXCLUDE    = 1
  Integer( c_int ), Private :: SEG_UREAD      = 4 * 8 * 8
  Integer( c_int ), Private :: SEG_UWRITE     = 2 * 8 * 8
  Integer( c_int ), Private :: SEG_GREAD      = 4 * 8
  Integer( c_int ), Private :: SEG_GWRITE     = 2 * 8
  Integer( c_int ), Private :: SEG_WREAD      = 4
  Integer( c_int ), Private :: SEG_WWRITE     = 2
  Integer( c_int ), Private :: SEG_NOREADONLY = 0
  Integer( c_int ), Private :: SEG_READONLY   = 1

  ! Elements of seg inquire array
  Integer( c_long ), Parameter, Private :: SEG_SEGSZ  = 1 ! size of segment in bytes
  Integer( c_long ), Parameter, Private :: SEG_LPID   = 2 ! process ID of last shared memory operation
  Integer( c_long ), Parameter, Private :: SEG_CPID   = 3 ! process ID of creator
  Integer( c_long ), Parameter, Private :: SEG_NATTCH = 4 ! number of current attaches
  Integer( c_long ), Parameter, Private :: SEG_ATIME  = 5 ! time of last shmat()
  Integer( c_long ), Parameter, Private :: SEG_DTIME  = 6 ! time of last shmdt()
  Integer( c_long ), Parameter, Private :: SEG_CTIME  = 7 ! time of last change by shmctl()

  ! Parameterized kind for use in functions taking an array "shape"
  Integer, Parameter, Private :: shape_t = kind(1)
  ! For differentiating between data types
  Integer, Parameter, Private :: integer_vals = 1
  Integer, Parameter, Private :: double_vals  = 2
  Integer, Parameter, Private :: complex_vals = 3

  ! Debugging Flag. Set to false for production
!  Logical, Parameter, Private :: debug = .True.
  Logical, Parameter, Private :: debug = .False.

  ! Derived type for storing data about a segment
  Type, Private :: segment
     Integer( c_int  )                              :: shmid
     Type( c_ptr     )                              :: shmaddr
     Type( FIPC_ctxt )                              :: ctxt
     Integer                                        :: type
     Integer( shape_t ), Dimension( : ), Allocatable :: sizes
  End Type segment

  ! Derived type for linked list for saving data about created segments
  Type, Private :: segment_list_type
     Type( segment           )          :: data
     Type( segment_list_type ), Pointer :: next => Null()
  End Type segment_list_type

  ! Linked list of data about created segments
  Type( segment_list_type ), Pointer, Private :: seg_list => Null()

  ! Derived type for storing data about a semaphore
  Type, Private :: semaphore
     Integer( c_int  ) :: semid
     Type( FIPC_ctxt ) :: ctxt
  End Type semaphore

  ! Derived type for linked list for saving data about created semaphores
  Type, Private :: semaphore_list_type
     Integer( c_int )                     :: semid
     Integer                              :: intra_handle
     Type( semaphore_list_type ), Pointer :: next => Null()
  End Type semaphore_list_type

  ! Linked list of data about created segements
  Type( semaphore_list_type ), Pointer, Private :: sem_list => Null()

  ! Largest allowed size for temporary buffers
  Integer, Parameter, Private :: max_buff_size = 1024 * 1024 / 8 ! 1 Mbyte of reals

  ! Interfaces to C wrappers.
  Interface

     Subroutine FIPC_get_errval( EACCES, EEXIST, EINVAL, ENFILE, ENOENT, ENOMEM, ENOSPC ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
!Changed intent to Out instead of In. The value gets set in a C routine

        Integer( c_int ), Intent( Out ) :: EACCES
        Integer( c_int ), Intent( Out ) :: EEXIST
        Integer( c_int ), Intent( Out ) :: EINVAL
        Integer( c_int ), Intent( Out ) :: ENFILE
        Integer( c_int ), Intent( Out ) :: ENOENT
        Integer( c_int ), Intent( Out ) :: ENOMEM
        Integer( c_int ), Intent( Out ) :: ENOSPC

     End Subroutine FIPC_get_errval

     Function FIPC_node_initialized() Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int ) :: FIPC_node_initialized
     End Function FIPC_node_initialized
     
     Function FIPC_get_seg( fcomm, size, disp, shmaddr, shmid ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int, c_long, c_ptr
       Implicit None
       Integer( c_int  )                      :: FIPC_get_seg
       Integer( c_int  ), Value, Intent( In ) :: fcomm
       Integer( c_long ), Value, Intent( In ) :: size
       Integer( c_long ), Value, Intent( In ) :: disp
       Type   ( c_ptr  )                      :: shmaddr
       Integer( c_int  ), Intent (Out)        :: shmid
     End Function FIPC_get_seg



     Function FIPC_crit_start( semid, rank ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int )                      :: FIPC_crit_start
       Integer( c_int ), Value, Intent( In ) :: semid, rank
     End Function FIPC_crit_start

     Function FIPC_crit_end( semid, rank ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int )                      :: FIPC_crit_end
       Integer( c_int ), Value, Intent( In ) :: semid, rank
     End Function FIPC_crit_end


     Function FIPC_sizeof_c_int() Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int ) :: FIPC_sizeof_c_int
     End Function FIPC_sizeof_c_int

     Function FIPC_sizeof_c_long() Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int ) :: FIPC_sizeof_c_long
     End Function FIPC_sizeof_c_long

     Function FIPC_sizeof_c_double() Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int ) :: FIPC_sizeof_c_double
     End Function FIPC_sizeof_c_double

     Function FIPC_sizeof_c_complex() Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int ) :: FIPC_sizeof_c_complex
     End Function FIPC_sizeof_c_complex

  End Interface

  ! Interfaces for context comparison function
  Interface operator( == )
     Module Procedure test_eq_ctxt
  End Interface
  Interface operator( /= )
     Module Procedure test_ne_ctxt
  End Interface

  ! Interfaces for communicator comparison function
  Interface operator( == ) 
     Module Procedure test_eq_comm
  End Interface
  Interface operator( /= )
     Module Procedure test_ne_comm
  End Interface

  !!!!> \endcond

  ! Overloaded interfaces

!Nitya, uncommented FIPC_allreduce
  Interface FIPC_allreduce
     Module Procedure allreduce_double
  End Interface

  ! alias name in the new initalisation scheme
  Interface FIPC_ctxt_init
     Module Procedure FIPC_init
  End Interface 

  Interface FIPC_ctxt_finalize
     Module Procedure FIPC_finalize
  End Interface

Contains


  Subroutine FIPC_node_init_world(error)
    Implicit none

    Integer, Intent(   Out ) :: error
     
     error = FIPC_success

   End Subroutine FIPC_node_init_world


   Subroutine FIPC_node_init_comm(comm, error)
     Implicit none
     
    Integer, Intent(   In  ) :: comm
    Integer, Intent(   Out ) :: error
     
    ! create shared memory structure only on comm
    ! the programmer must ensure that comm covers 
    ! all nodes that plan to use shared memory
    error = FIPC_success
    
  End Subroutine FIPC_node_init_comm
  

  !> \brief Initialise FreeIPC. 
  !!
  !! This routine initializes FreeIPC. It sets up a context spanning all process in the MPI communicator 
  !! universe_comm, i.e. FIPC_ctxt_world, and determines which processes are on the smae shared memory node.
  !! \param universe_comm The communicator that spans all the processes within which FreeIPC will work 
  !! \param error         On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>

  Subroutine FIPC_init( universe_comm, error )

    Integer, Intent( In    ) :: universe_comm
    Integer, Intent(   Out ) :: error

    Integer :: world_comm, intra_comm, extra_comm

    ! Check if the node data structure are initialised
    ! relevant for posix version

    !If ( fipc_node_initialized() == 0) then
    !   error = FIPC_node_uninitialized
    !   Return
    !End If

    ! Can only initialise once
    If( Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_already_initialized
       Return
    End If


    ! From the universe communicator generate the 3 communicators required
    ! to specify FIPC_ctxt_world
    Call generate_base_comms( universe_comm, world_comm, intra_comm, extra_comm, error )
    If( error /= FIPC_success ) Then
       Return
    End If

    ! From those communicators generate FIPC_ctxt_world
    Call ctxt_create( world_comm, intra_comm, extra_comm, FIPC_ctxt_world, error )

    If( error /= FIPC_success ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_init

  !> \brief Test if FreeIPC is initialised. 
  !!
  !! Test if FreeIPC is initialised. 
  !! \param flag  On succesful exit this is set to .true. if FIPC is initialised, , .false. otherwise
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_initialized( flag, error )

    ! Test whether FreeIPC has been initialized. 
    ! On success return .True. in FLAG if is has, otherwise .FALSE., and 
    ! set ERROR to FIPC_SUCCESS
    ! On error FLAG is undefined and ERROR is any value but FIPC_SUCCESS

    Logical, Intent( Out ) :: flag
    Integer, Intent( Out ) :: error

    flag = Associated( FIPC_ctxt_world%initialized )

    error = FIPC_success

  End Subroutine FIPC_initialized

  !> \brief Finalise FreeIPC.
  !!
  !! This routine closes FreeIPC down.  It frees all shared memory segments and semaphores
  !! that the user has created through FreeIPC, and frees FIPC_ctxt_world.
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_finalize( error )

    ! Finalize FIPC.
    !
    ! On success ERROR is set to FIPC_SUCCESS. Any other value
    ! indicates error. These can be compared to the symbolic constants
    ! defined above for better diagnosis

    Integer, Intent( Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Free all the segments we know about
    Do While( Associated( seg_list ) )
       Call segment_free( seg_list%data%shmid, error )
       If( error /= FIPC_success ) Then
          Return
       End If
    End Do

    ! And free the context
    Call ctxt_free( FIPC_ctxt_world, error )
    If( error /= 0 ) Then
       Return
    End If

!!$    ! Finally free any outstanding semaphores we know about
!!$    Do While( Associated( sem_list ) )
!!$       Call semaphore_free( sem_list%semid, error )
!!$       If( error /= FIPC_success ) Then
!!$          Return
!!$       End If
!!$    End Do

!!$    Call node_key_and_sem_free(error)
!!$    If( error /= FIPC_success ) Then
!!$       Return
!!$    End If

    error = FIPC_success

  End Subroutine FIPC_finalize


  Subroutine FIPC_node_finalize(error)
    Implicit none

!!$    Interface
!!$       Function FIPC_node_keys_and_sems_free() Bind( C )
!!$         Use, Intrinsic :: iso_c_binding, Only : c_int
!!$         Implicit None
!!$         Integer( c_int ) :: FIPC_node_keys_and_sems_free
!!$       End Function FIPC_node_keys_and_sems_free
!!$    End Interface
!!$
    Integer, Intent( Out ) :: error
!!$
!!$    if ( fipc_node_keys_and_sems_free() /= 0) then
!!$       error = FIPC_undefined
!!$       Return
!!$    End If

    error = FIPC_success

  End Subroutine FIPC_node_finalize

  !> \brief Duplicate a context
  !!
  !! This routine duplicates the context ctxt_1, returning a new context in ctxt_2. This is very similar to
  !! mpi_comm_dup - see http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Comm_dup.html
  !! \param ctxt_1 The input context
  !! \param ctxt_2 The result context
  !! \param error  On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_ctxt_dup( ctxt_1, ctxt_2, error )

    ! Duplicate the context CTXT_1, returning th new context in CTXT_2
    !
    ! On success ERROR is set to FIPC_SUCCESS. Any other value
    ! indicates error. These can be compared to the symbolic constants
    ! defined above for better diagnosis

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt_1
    Type( FIPC_ctxt ), Intent(   Out ) :: ctxt_2
    Integer          , Intent(   Out ) :: error

    Integer :: world_comm_2
    Integer :: intra_comm_2
    Integer :: extra_comm_2

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Duplicate the communicators
    Call mpi_comm_dup( ctxt_1%world_comm%handle, world_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_dup( ctxt_1%intra_comm%handle, intra_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    ! Extra comm only defined on rank 0 of the intra comm
    If( ctxt_1%intra_comm%rank == 0 ) Then
       Call mpi_comm_dup( ctxt_1%extra_comm%handle, extra_comm_2, error )
       If( error /= 0 ) Then
          Return
       End If
    End If

    ! Now create the new context

    Call ctxt_create( world_comm_2, intra_comm_2, extra_comm_2, ctxt_2, error )
    If( error /= 0 ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_dup

  !> \brief Split a context
  !!
  !! This routine splits the context CTXT_1 according to the colours and keys given
  !! by COLOUR and KEY, with the. Processes with the same value of COLOUR end up
  !! in the same new context, the rank order being controlled by KEY. The special value FIPC_undefined
  !! may be used to indicate that a process should not be in any of the new contexts that are
  !! created by this routine. This is very
  !! similar to MPI_COMM_SPLIT, see http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Comm_split.html
  !! \param ctxt_1 The input context
  !! \param colour Control of subset assignment (nonnegative integer or FIPC_undefined). Processes with the same color are 
  !! in the same new context
  !! \param key    Control of rank assigment
  !! \param ctxt_2 The result context
  !! \param error  On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_ctxt_split( ctxt_1, colour, key, ctxt_2, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt_1
    Integer          , Intent( In    ) :: colour
    Integer          , Intent( In    ) :: key
    Type( FIPC_ctxt ), Intent(   Out ) :: ctxt_2
    Integer          , Intent(   Out ) :: error

    Integer :: world_comm_2
    Integer :: intra_comm_2
    Integer :: extra_comm_2

    Integer :: intra_rank, is_node_0

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Split the communicators
    Call mpi_comm_split( ctxt_1%world_comm%handle, colour, key, world_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_split( ctxt_1%intra_comm%handle, colour, key, intra_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    ! Extra comm only defined on rank 0 of the intra comm
    ! BUT WHAT IF SPLIT WITHIN THE INTRA COMM. Need to think here ....
    Call mpi_comm_rank( intra_comm_2, intra_rank, error )
    If( error /= 0 ) Then
       Return
    End If
    is_node_0 = merge( 1, MPI_UNDEFINED, intra_rank == 0 )
!    Call mpi_comm_split( world_comm_2, colour, key, extra_comm_2, error )

!TESTING, THE COMM SPLIT FOR EXTRA COMM SHOULD BE OVER INTRA_RANK==0 AND NOT COLOR
    Call  mpi_comm_split( world_comm_2, is_node_0, key, extra_comm_2, error )
!END TESTING 

    If( error /= 0 ) Then
       Return
    End If

    ! Now create the new context if required
    If( colour /= FIPC_undefined ) Then
       Call ctxt_create( world_comm_2, intra_comm_2, extra_comm_2, ctxt_2, error )
    End If
    If( error /= 0 ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_split

  !> \brief Free a context
  !!
  !! This routine frees the context CTXT. Compare MPI_COMM_FREE - 
  !! http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Comm_free.html
  !! \param ctxt  The context to be freed
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_ctxt_free( ctxt, error )

    Type( FIPC_ctxt ), Intent( InOut ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Type( segment_list_type ), Pointer :: p

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Check not trying to free the base context
    If( ctxt == FIPC_ctxt_world ) Then
       error = FIPC_freeing_ctxt_world
       Return
    End If

    ! Check if anybody is still using this context
    p => seg_list
    Do While( Associated( p ) )
       If( p%data%ctxt == ctxt ) Then
          error = FIPC_ctxt_in_use
          Return
       End If
       p => p%next
    End Do

    Call ctxt_free( ctxt, error )
    If( error /= FIPC_SUCCESS ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_free

  !> \brief Extract the intra-node communicator
  !!
  !! This routine extracts from the context CTXT a communicator that spans the processes
  !! within the context on the same shared memory node. The general syntax is similar to MPI_COMM_GROUP
  !! http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Comm_group.html
  !! \param ctxt  The context
  !! \param comm  The intra-node communicator
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_ctxt_intra_comm( ctxt, comm, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: comm
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    comm = ctxt%intra_comm%handle

    error = FIPC_success

  End Subroutine FIPC_ctxt_intra_comm

  !> \brief Extract the world communicator
  !!
  !! This routine extracts from the context CTXT a communicator that spans all the processes
  !! within the context. The general syntax is similar to MPI_COMM_GROUP
  !! http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Comm_group.html
  !! \param ctxt  The context
  !! \param comm  The world communicator
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_ctxt_world_comm( ctxt, comm, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: comm
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    comm = ctxt%world_comm%handle

    error = FIPC_success

  End Subroutine FIPC_ctxt_world_comm

  !> \brief Extract the extra-node communicator
  !!
  !! This routine extracts from the context CTXT a communicator that spans the processes
  !! within the context that have rank 0 in a intra-communicator. If the process does not
  !! have rank 0 in the intra-communicator FIPC_comm_null is returned. The general syntax is similar to MPI_COMM_GROUP
  !! http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Comm_group.html
  !! \param ctxt  The context
  !! \param comm  The extra-node communicator
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_ctxt_extra_comm( ctxt, comm, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: comm
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    If( ctxt%intra_comm%rank == 0 ) Then
       comm = ctxt%extra_comm%handle
    Else
       comm = FIPC_comm_null
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_extra_comm

  !> \brief Start a crtitical region
  !!
  !! This routine starts a critical region for the processes spanned by the context CTXT.
  !! The nett result is that within a shared memory node only one process can be executing
  !! code within a critical region. Note there is nothing wrong with nested critical regions.
  !! Also see FIPC_critical_end
  !! \param ctxt  The context
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_critical_start( ctxt, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! is ctxt%semid a fortran integer?
    call MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, ctxt%semid, error)
    !error = fipc_crit_start( ctxt%semid, ctxt%intra_comm%rank )
    
    If( error /= 0 ) Then
       error = FIPC_critical_start_failed
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_critical_start

  !> \brief End a crtitical region
  !!
  !! This routine ends a critical region for the processes spanned by the context CTXT.
  !! Also see FIPC_critical_start
  !! \param ctxt  The context
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_critical_end( ctxt, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Is ctxt%semid a fortran integer ?
    call MPI_Win_unlock(0, ctxt%semid, error)
    !error = fipc_crit_end( ctxt%semid, ctxt%intra_comm%rank )
    
    If( error /= 0 ) Then
       error = FIPC_critical_end_failed
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_critical_end

  !> \brief Synchronize within a shared memory node
  !!
  !! This routine synchronises the processes on the smae shared memory node that are
  !! in the context CTXT
  !! \param ctxt  The context
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error. 
  !>
  Subroutine FIPC_node_barrier( ctxt, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Call mpi_barrier( ctxt%intra_comm%handle, error )

  End Subroutine FIPC_node_barrier

  Subroutine generate_base_comms( universe_comm, world_comm, intra_comm, extra_comm, error )

    ! Try to work out the mapping of the processors onto the
    ! physical nodes, and from that generate the communicators
    ! needed to set up FIPC_ctxt_world

    Integer, Intent( In    ) :: universe_comm
    Integer, Intent(   Out ) :: world_comm
    Integer, Intent(   Out ) :: intra_comm
    Integer, Intent(   Out ) :: extra_comm
    Integer, Intent(   Out ) :: error

    ! First set the world comm. 
    Call generate_base_world_comm( universe_comm, world_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Now the intra node communicator
    Call generate_base_intra_comm( world_comm, intra_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    ! We now have an intra comm. From that it's easy to generate an extra comm.
    Call generate_base_extra_comm( world_comm, intra_comm, extra_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    If( debug ) Then
       Call write_base_info( intra_comm, extra_comm, error )
       If( error /= 0 ) Then
          Return
       End If
    End If

    error = FIPC_success

  Contains

    Subroutine generate_base_world_comm( universe_comm, world_comm, error )

      ! From the universe communicator generate the base world communicator.
      ! Trivial !

      Integer, Intent( In    ) :: universe_comm
      Integer, Intent(   Out ) :: world_comm
      Integer, Intent(   Out ) :: error

      Call mpi_comm_dup( universe_comm, world_comm, error )
      If( error /= 0 ) Then
         Return
      End If

      error = 0

    End Subroutine generate_base_world_comm

    Subroutine generate_base_intra_comm( world_comm, intra_comm, error )

      ! Generate the intra node communicator. This is the difficult one.
      ! The basic idea is to find which other processes can attach to
      ! a shared memory segment created by a reference process.

      Integer, Intent( In    ) :: world_comm
      Integer, Intent(   Out ) :: intra_comm
      Integer, Intent(   Out ) :: error
      
      Integer( c_int ) :: colour

      Integer :: world_rank, world_size

!      error = fipc_generate_intra_colour( Int( world_comm, c_int ), colour )
!      If( error /= 0 ) Then
!         Return
!      End If
      
      Call mpi_comm_rank( world_comm, world_rank, error )
!      Call mpi_comm_size( world_comm, world_size, error )

      Call mpi_comm_split_type( world_comm, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, intra_comm, error )

      if (error /= MPI_SUCCESS) then
         error = -1 ! needed to improve this
      else
         error = 0
      endif
     
    End Subroutine generate_base_intra_comm

    Subroutine generate_base_extra_comm( world_comm, intra_comm, extra_comm, error )

      ! From the world and intra communicators generate the extra communicator
      ! Note that this is defined only on processor zero of each multicore node.

      Integer, Intent( In    ) :: world_comm
      Integer, Intent( In    ) :: intra_comm
      Integer, Intent(   Out ) :: extra_comm
      Integer, Intent(   Out ) :: error

      Integer :: world_rank
      Integer :: intra_rank
      Integer :: split_val

      Call mpi_comm_rank( world_comm, world_rank, error )
      If( error /= MPI_SUCCESS ) Then
         Return
      End If
      Call mpi_comm_rank( intra_comm, intra_rank, error )
      If( error /= MPI_SUCCESS ) Then
         Return
      End If
      split_val = Merge( 1, MPI_UNDEFINED, intra_rank == 0 )
      Call mpi_comm_split( world_comm, split_val, world_rank, extra_comm, error )
      If( error /= MPI_SUCCESS ) Then
         Return
      End If
      
      error = 0

    End Subroutine generate_base_extra_comm

    Subroutine write_base_info( intra_comm, extra_comm, error )

      Integer, Intent( In    ) :: intra_comm
      Integer, Intent( In    ) :: extra_comm
      Integer, Intent(   Out ) :: error

      Integer :: intra_size, intra_rank
      Integer :: extra_size, extra_rank
      Integer :: i

      ! Write out a bit more data
      Call mpi_comm_size( intra_comm, intra_size, error )
      If( error /= 0 ) Then
         Return
      End If
      Call mpi_comm_rank( intra_comm, intra_rank, error )
      If( error /= 0 ) Then
         Return
      End If
      If( intra_rank == 0 ) Then
         Call mpi_comm_size( extra_comm, extra_size, error )
         If( error /= 0 ) Then
            Return
         End If
         Call mpi_comm_rank( extra_comm, extra_rank, error )
         If( error /= 0 ) Then
            Return
         End If
         If( extra_rank == 0 ) Then
            Write( *, * ) 'Number of nodes found: ', extra_size
         End If
         Do i = 0 , extra_size
            Call mpi_barrier( extra_comm, error )
            If( i == extra_rank ) Then
               Write( *, * ) 'Number of processes on node ', extra_rank, ' is ', intra_size
            End If
         End Do
      End If

      error = 0

    End Subroutine write_base_info

  End Subroutine generate_base_comms

!!$  Subroutine get_new_seg(comm, n, shmaddr, shmid )
!!$
!!$    ! Get a segment of size N bytes in exclusive mode
!!$    !
!!$    ! On success the segment id is returned in SHMID
!!$    ! On error SHMID is negative. Possible values are
!!$    ! the negative of the values that errno can be set to
!!$    ! by shmget ( see 
!!$    !
!!$    ! http://www.opengroup.org/onlinepubs/009695399/functions/shmget.html
!!$    !
!!$    ! ) if the shmget failed, or - Huge( 1 ) if despite trying
!!$    ! multiple times we failed to create a seg
!!$
!!$    Integer, Intent( In              ) :: comm
!!$    Integer( c_long ), Intent( In    ) :: n
!!$    Type( c_ptr     )                  :: shmaddr  ! No intent for Cray compiler
!!$    Integer( c_int  ), Intent(   Out ) :: shmid
!!$
!!$    Integer, Parameter :: max_tries = 200
!!$
!!$    Integer            :: tries
!!$
!!$    tries = 1
!!$
!!$    Do
!!$       ! Try to create a seg in exlusive mode. Note that fipc_seg_create
!!$       ! tries to use a different value of shmid each time it is called.
!!$       shmid = fipc_get_seg( comm, n, SEG_CREATE, SEG_EXCLUDE, SEG_UREAD + SEG_UWRITE )
!!$       If( shmid /= EEXIST ) Then
!!$          ! Either fipc_get_seg was succesful, or it failed for a reason other
!!$          ! than the seg already exisiting
!!$          Exit
!!$       End If
!!$       ! The seg we tried to create already exists. Try again ...
!!$       tries = tries + 1
!!$       If( tries > max_tries ) Then
!!$          ! ... unless we've got bored. Make sure we can't go on for ever;
!!$          ! As in the old geek joke we're sort of putting a known elephant in Cairo.
!!$          ! http://en.wikipedia.org/wiki/Elephant_in_Cairo
!!$          shmid = - Huge( tries )
!!$          Exit
!!$       End If
!!$    End Do
!!$
!!$  End Subroutine get_new_seg


  Subroutine ctxt_create( world_comm, intra_comm, extra_comm, ctxt, error )

    ! From the basic data create a context.
    !
    ! On success ERROR is set to FIPC_success and ctxt holds the new context
    ! On error error is any other value than FIPC_success

    Integer          , Intent( In    ) :: world_comm
    Integer          , Intent( In    ) :: intra_comm
    Integer          , Intent( In    ) :: extra_comm
    Type( FIPC_ctxt ), Intent(   Out ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Integer( c_int ) :: semid

    Allocate( ctxt%world_comm, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    Call set_comm( world_comm, ctxt%world_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    Allocate( ctxt%intra_comm, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    Call set_comm( intra_comm, ctxt%intra_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    ! The extra communicator is only defined on processor zero 
    ! of the intra processors. However it is useful to know the 
    ! extra rank and size of the zero proc in the intra comm on all procs 
    ! in the inter comm - this is the "node" number for this multicore 
    ! node in the universe, and the number of nodes. Hence allocate the extra comm on all procs
    Allocate( ctxt%extra_comm, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    If( ctxt%intra_comm%rank == 0 ) Then
       Call set_comm( extra_comm, ctxt%extra_comm, error )
       If( error /= 0 ) Then
          Return
       End If
    End If
    Call mpi_bcast( ctxt%extra_comm%rank, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_bcast( ctxt%extra_comm%size, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Now set up the semaphore for this context
!!$    Call semaphore_create( ctxt, semid, error )
!!$    If( error /= FIPC_success ) Then
!!$       Return
!!$    End If
    ctxt%semid = semid

!TESTING 
!   print *, ' wr ', ctxt%world_comm%rank, ' ws ',  ctxt%world_comm%size, ' ir ',  ctxt%intra_comm%rank, ' is ',  ctxt%intra_comm%size, &
!         ' er ', ctxt%extra_comm%rank, ' es ',  ctxt%extra_comm%size
!END TESTING

!Added code Nitya
!We mark the segment to be deleted, this will be done
!after all the processes have detached from the segment.

!   If( ctxt%intra_comm%rank == 0 ) Then

   !Write the semaphore id out to a file, incase the semaphore is not deleted after abnormal execution,
   !They will have to removed forcefully by ipcrm
   !The file has to be written out by the MPI process with rank 0 in the root communicator.

!!    call FIPC_write_sms( ctxt%extra_comm, semid, 1, ctxt=ctxt)
!   End If

   Allocate( ctxt%initialized, Stat = error )
   If( error /= 0 ) Then
      error = FIPC_allocation_failed
      Return
   End If

   error = FIPC_success

  Contains

    Subroutine set_comm( handle, comm, error )

      Integer             , Intent( In    ) :: handle
      Type( communicator ), Intent(   Out ) :: comm
      Integer             , Intent(   Out ) :: error

      comm%handle = handle

      Call mpi_comm_size( comm%handle, comm%size, error )
      If( error /= 0 ) Then
         Return
      End If

      Call mpi_comm_rank( comm%handle, comm%rank, error )
      If( error /= 0 ) Then
         Return
      End If

      comm%initialized = .True.

    End Subroutine set_comm

  End Subroutine ctxt_create

  Subroutine ctxt_free( ctxt, error )

    ! Free a context.
    !
    ! On success ERROR is set to FIPC_success 
    ! On error error is any other value than FIPC_success

    Type( FIPC_ctxt ), Intent( InOut ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Integer :: comm

!!$    ! First free the semaphore
!!$    Call semaphore_free( ctxt%semid, error )
!!$    If( error /= 0 ) Then
!!$       Return
!!$    End If

    ! And now the communicators
    If( ctxt%intra_comm%rank == 0 ) Then
       Call mpi_comm_free( ctxt%extra_comm%handle, error )
       If( error /= 0 ) Then
          Return
       End If
       ctxt%extra_comm%initialized = .False.
       Deallocate( ctxt%extra_comm )
    End If

    ! Take a copy of the intra_node communicator so can sync at end
    Call mpi_comm_dup( ctxt%intra_comm%handle, comm, error ) 
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_free( ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If
    ctxt%intra_comm%initialized = .False.
    Deallocate( ctxt%intra_comm )

    Call mpi_comm_free( ctxt%world_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If
    ctxt%world_comm%initialized = .False.
    Deallocate( ctxt%world_comm )

    ! Make sure all in sync on the node, carefully using the
    ! copy of the intra node context
    Call mpi_barrier( comm, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_free( comm, error )
    If( error /= 0 ) Then
       Return
    End If

    Deallocate( ctxt%initialized )

    error = FIPC_success

  End Subroutine ctxt_free


  Subroutine segment_create( what, rank, sizes, ctxt, shmaddr, error )

    ! Create a segment of size N bytes, attach it to all procs
    ! on this node in the context CTXT, and return a c pointer
    ! to it in SHMADDR. On success ERROR is equal to FIPC_SUCCESS,
    ! any other value indicates error. On failure SHMADDR is also
    ! set to C_NULL_PTR

    Integer                          , Intent( In    ) :: what
    Integer                          , Intent( In    ) :: rank
    Integer( shape_t ), Dimension( : ), Intent( In    ) :: sizes
    Type( FIPC_ctxt )                , Intent( InOut  ) :: ctxt
    Type( c_ptr     )                                  :: shmaddr ! No intent for Cray compiler
    Integer                          , Intent(   Out ) :: error

    Type( segment_list_type ), Pointer :: p, this

    Integer( c_int ) :: shmid
    Integer :: rm_seg

     error=FIPC_success

    ! Create the seg on proc 0 in this context on this node
    !If( ctxt%intra_comm%rank == 0 ) Then
    If ( fipc_get_seg( ctxt%intra_comm%handle, &
         sizeof_what( what ) * Product( sizes ), sizeof_what(what), &
         shmaddr, shmid) /= 0) then
       error = -1
    End If

    ! Tell the other procs on the node
    !Call mpi_bcast( shmid, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    !If( error /= 0 ) Then
    !   shmaddr = c_null_ptr
    !   Return
    !End If
    
    If( error /= 0 ) Then
       error = FIPC_seg_get_failed
       shmaddr = c_null_ptr
       Return
    End If
    
    ! use MPI3 win for locks instead of semaphores
    ctxt%semid = shmid

    ! Attach to the seg
    !shmaddr = fipc_attach_seg( shmid, SEG_NOREADONLY )
    If( debug ) Then
       If (c_associated( shmaddr )) then
          Write( *, * ) 'assoc status of seg in segment_create: true'
       Else
          Write( *, * ) 'assoc status of seg in segment_create: false'
       End If
    End If
    If( .Not.  c_associated( shmaddr ) ) Then 
       error = FIPC_seg_attach_failed
       shmaddr = c_null_ptr
       Return
    End If

    ! Add data on the seg to the end of linked list about created segs
    p => seg_list
    If( .Not. Associated( p ) ) Then
       ! First element
       Allocate( seg_list, Stat = error )
       If( error /= 0 ) Then
          error = FIPC_allocation_failed
          shmaddr = c_null_ptr
          Return
       End If
       this => seg_list
    Else
       Do While( Associated( p%next ) )
          p => p%next
       End Do       
       Allocate( p%next, Stat = error )
       If( error /= 0 ) Then
          error = FIPC_allocation_failed
          shmaddr = c_null_ptr
          Return
       End If
       this => p%next
    End If
       
    this%data%shmid   = shmid
    this%data%shmaddr = shmaddr
    this%data%ctxt    = ctxt
    this%data%type    = what
    Allocate( this%data%sizes( 1:rank ), Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    this%data%sizes = sizes( 1:rank )
    this%next         => Null()

    ! Make sure all up to date
    Call mpi_barrier( ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       shmaddr = c_null_ptr
       Return
    End if

    If( debug ) Then
       Call print_seg_list( ctxt )
    End If

!Added code Nitya
!We mark the segment to be deleted, this will be done
!after all the processes have detached from the segment.

!!$    If( ctxt%intra_comm%rank == 0 ) Then
!!$     rm_seg = fipc_remove_seg(shmid)
!!$
!!$     If(rm_seg /= 0) Then
!!$       error = FIPC_sem_remove_failed
!!$       Return
!!$     End If
!!$
!!$   !Write the segment id out to a file, incase the segments are not deleted after abnormal execution,
!!$   !They will have to removed forcefully by ipcrm
!!$   !The file has to be written out by the MPI process with rank 0 in the root communicator.
!!$ 
!!$!!    call FIPC_write_sms( ctxt%extra_comm, shmid , ctxt=ctxt)
!!$   End If
!!$
!!$    ! Make sure all up to date
!!$    Call mpi_barrier( ctxt%intra_comm%handle, error )
!!$    If( error /= 0 ) Then
!!$       shmaddr = c_null_ptr
!!$       Return
!!$    End if

    error = FIPC_success

  End Subroutine segment_create

  Subroutine segment_free( shmid, error )

    ! Free a segment

    Integer( c_int ), Intent( In   ) :: shmid
    Integer         , Intent(  Out ) :: error

    Type( segment_list_type ), Pointer :: p, prev

    Type( c_ptr ) :: shmaddr

    Integer :: retval

    ! First find the seg
    prev => Null()
    p    => seg_list
    Do While( Associated( p ) )
       If( p%data%shmid == shmid ) Then
          Exit
       End If
       prev => p
       p    => p%next
    End Do
    If( .Not. Associated( p ) ) Then
       error = FIPC_seg_not_found
       Return
    End If

    If( debug ) Then
       ! Not really neccessary, but can make debug printing a bit clearer if
       ! things get too out of step
       Call mpi_barrier( p%data%ctxt%intra_comm%handle, error )
    End If

    ! inconsistent, shmid is generated in C !!!
    call MPI_Win_free(shmid, error)
    If ( error /= MPI_SUCCESS ) Then
       return
    End If

    ! And remove the segment from the linked list
    ! Link up list around the one to die
    If( Associated( prev ) ) Then
       ! Not first in list
       prev%next => p%next
    Else
       ! First in list
       seg_list => p%next
    End If
    Deallocate( p%data%sizes )
    Deallocate( p )
    p => Null()

    error = FIPC_success

  End Subroutine segment_free

!!$  Subroutine semaphore_create( ctxt, semid, error )
!!$
!!$    ! Create a semaphore of size N bytes in the context CTXT, and return 
!!$    ! it's id in SEMID. On success ERROR is equal to FIPC_SUCCESS,
!!$    ! any other value indicates error. 
!!$
!!$    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
!!$    Integer( c_int ) , Intent(   Out ) :: semid
!!$    Integer          , Intent(   Out ) :: error
!!$
!!$    Type( semaphore_list_type ), Pointer :: p, this
!!$
!!$    ! Create the sem on proc 0 in this context on this node
!!$    If( ctxt%intra_comm%rank == 0 ) Then
!!$       Call get_new_sem( semid )
!!$       If( debug ) Then
!!$          Write( *, * ) 'semid ', semid
!!$       End If
!!$       If( semid < 0 ) Then
!!$          error = FIPC_sem_get_failed
!!$          Return
!!$       End If
!!$    End If
!!$
!!$    ! Tell the other procs on the node
!!$    Call mpi_bcast( semid, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
!!$    If( error /= 0 ) Then
!!$       Return
!!$    End If
!!$
!!$    ! Add data on the sem to the end of linked list about created sems
!!$    p => sem_list
!!$    If( .Not. Associated( p ) ) Then
!!$       ! First element
!!$       Allocate( sem_list, Stat = error )
!!$       If( error /= 0 ) Then
!!$          error = FIPC_allocation_failed
!!$          Return
!!$       End If
!!$       this => sem_list
!!$    Else
!!$       Do While( Associated( p%next ) )
!!$          p => p%next
!!$       End Do       
!!$       Allocate( p%next, Stat = error )
!!$       If( error /= 0 ) Then
!!$          error = FIPC_allocation_failed
!!$          Return
!!$       End If
!!$       this => p%next
!!$    End If
!!$       
!!$    this%semid        = semid
!!$    this%intra_handle = ctxt%intra_comm%handle
!!$
!!$    ! Make sure all up to date
!!$    Call mpi_barrier( ctxt%intra_comm%handle, error )
!!$    If( error /= 0 ) Then
!!$       Return
!!$    End if
!!$
!!$    error = FIPC_success
!!$
!!$  End Subroutine semaphore_create
!!$
!!$  Subroutine semaphore_free( semid, error )
!!$
!!$    ! Free a semaphore
!!$
!!$    Integer( c_int ), Intent( In   ) :: semid
!!$    Integer         , Intent(  Out ) :: error
!!$
!!$    Type( semaphore_list_type ), Pointer :: p, prev
!!$
!!$    Integer :: retval
!!$    Integer :: rank
!!$
!!$    ! First find the sem
!!$    prev => Null()
!!$    p    => sem_list
!!$    Do While( Associated( p ) )
!!$       If( p%semid == semid ) Then
!!$          Exit
!!$       End If
!!$       prev => p
!!$       p    => p%next
!!$    End Do
!!$    If( .Not. Associated( p ) ) Then
!!$       error = FIPC_sem_not_found
!!$       Return
!!$    End If
!!$
!!$    If( debug ) Then
!!$       ! Not really neccessary, but can make debug printing a bit clearer if
!!$       ! things get too out of step
!!$       Call mpi_barrier( p%intra_handle, error )
!!$    End If
!!$
!!$    ! And delete the semaphore
!!$    Call mpi_comm_rank( p%intra_handle, rank, error )
!!$    If( error /= 0 ) Then
!!$       Return
!!$    End if
!!$    If( rank == 0 ) Then
!!$       retval = fipc_remove_sem( semid )
!!$       If( retval /= 0 ) Then
!!$          error = FIPC_sem_remove_failed
!!$          Return
!!$       End If
!!$    End If
!!$
!!$    ! Make sure all up to date
!!$    Call mpi_barrier( p%intra_handle, error )
!!$    If( error /= 0 ) Then
!!$       Return
!!$    End if
!!$
!!$    ! And remove the semaphore from the linked list
!!$    ! Link up list around the one to die
!!$    If( Associated( prev ) ) Then
!!$       ! Not first in list
!!$       prev%next => p%next
!!$    Else
!!$       ! First in list
!!$       sem_list => p%next
!!$    End If
!!$    Deallocate( p )
!!$    p => Null()
!!$
!!$    error = FIPC_success
!!$
!!$  End Subroutine semaphore_free

!!$  subroutine node_key_and_sem_free(error)
!!$    implicit none
!!$    Integer :: error
!!$    
!!$    Interface
!!$       Function FIPC_node_keys_and_sems_free() Bind( C )
!!$         Use, Intrinsic :: iso_c_binding, Only : c_int
!!$         Implicit None
!!$         Integer( c_int ) :: FIPC_node_keys_and_sems_free
!!$       End Function FIPC_node_keys_and_sems_free
!!$    End Interface
!!$
!!$    error = fipc_node_keys_and_sems_free()
!!$    
!!$  end subroutine node_key_and_sem_free

  Subroutine print_seg_list( ctxt )

    Type( FIPC_ctxt ), Intent( In ) :: ctxt

    Type( segment_list_type ), Pointer :: p

    Integer( c_long ), Dimension( 1:7 ) :: buf

    Integer( c_int ) :: shmid

    Integer :: seg_n_at, rank
    Integer :: retval
    Integer :: mod_place

    Character( Len = 10 ) :: type
    Character( Len = 60 ) :: format

    If( ctxt%intra_comm%rank == 0 ) Then
       Write( *, '( "------ Shared Memory Segments --------" )' )
       Write( *, '( "shmid", t10, "type", t20, "Dimensions", t40, "Attached"  )' )
       format = '( i0, t10, a, t20, ?( i0, 1x ), t40, 1x, i0 )'
       mod_place = scan( format, '?' )
       p => seg_list
       Do While( Associated( p ) )
          shmid  = p%data%shmid
          !retval = fipc_inquire_seg( shmid, Size( buf ), buf )
          seg_n_at = buf( SEG_NATTCH )
          Select Case( p%data%type )
          Case( integer_vals )
             type = 'Integer'
          Case( double_vals )
             type = 'Double'
          Case( complex_vals )
             type = 'Complex'
          Case Default
             type = 'UNKNOWN'
          End Select
          rank = Size( p%data%sizes )
          Write( format( mod_place:mod_place ), '( i1 )' ) rank
          Write( *, format ) shmid, type, p%data%sizes, seg_n_at
          p => p%next
       End Do
    End If

  End Subroutine print_seg_list

  ! Now follow the segment creation routines

  Subroutine seg_create_integer_1d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )               , Intent( InOut ) :: ctxt
    Integer, Dimension( : )         , Intent( In    ) :: n
    Integer, Dimension( : ), Pointer, Intent(   Out ) :: a
    Integer                         , Intent(   Out ) :: error

    Integer, Parameter :: rank = 1
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_1d_size_in_int

  Subroutine seg_create_integer_2d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                  , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )            , Intent( In    ) :: n
    Integer, Dimension( :, : ), Pointer, Intent(   Out ) :: a
    Integer                            , Intent(   Out ) :: error

    Integer, Parameter :: rank = 2
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_2d_size_in_int

  Subroutine seg_create_integer_3d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                     , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )               , Intent( In    ) :: n
    Integer, Dimension( :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                               , Intent(   Out ) :: error

    Integer, Parameter :: rank = 3
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_3d_size_in_int

  Subroutine seg_create_integer_4d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                        , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                  , Intent( In    ) :: n
    Integer, Dimension( :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                  , Intent(   Out ) :: error

    Integer, Parameter :: rank = 4
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_4d_size_in_int

  Subroutine seg_create_integer_5d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                           , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                     , Intent( In    ) :: n
    Integer, Dimension( :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: rank = 5
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_5d_size_in_int

  Subroutine seg_create_integer_6d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                           , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                     , Intent( In    ) :: n
    Integer, Dimension( :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: rank = 6
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_6d_size_in_int

  Subroutine seg_create_integer_7d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                 , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                           , Intent( In    ) :: n
    Integer, Dimension( :, :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: rank = 7
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_integer_7d_size_in_int

  Subroutine seg_create_double_1d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                        , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                  , Intent( In    ) :: n
    Real( c_double ), Dimension( : ), Pointer, Intent(   Out ) :: a
    Integer                                  , Intent(   Out ) :: error

    Integer, Parameter :: rank = 1
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_1d_size_in_int

  Subroutine seg_create_double_2d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                           , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                     , Intent( In    ) :: n
    Real( c_double ), Dimension( :, : ), Pointer, Intent(   Out ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: rank = 2
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_2d_size_in_int

  Subroutine seg_create_double_3d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                              , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                        , Intent( In    ) :: n
    Real( c_double ), Dimension( :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                        , Intent(   Out ) :: error

    Integer, Parameter :: rank = 3
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_3d_size_in_int

  Subroutine seg_create_double_4d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                 , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                           , Intent( In    ) :: n
    Real( c_double ), Dimension( :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: rank = 4
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_4d_size_in_int

  Subroutine seg_create_double_5d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                    , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                              , Intent( In    ) :: n
    Real( c_double ), Dimension( :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                              , Intent(   Out ) :: error

    Integer, Parameter :: rank = 5
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_5d_size_in_int

  Subroutine seg_create_double_6d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                       , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                                 , Intent( In    ) :: n
    Real( c_double ), Dimension( :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                                 , Intent(   Out ) :: error

    Integer, Parameter :: rank = 6
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_6d_size_in_int

  Subroutine seg_create_double_7d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                          , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                                    , Intent( In    ) :: n
    Real( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                                    , Intent(   Out ) :: error

    Integer, Parameter :: rank = 7
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_double_7d_size_in_int

  Subroutine seg_create_complex_1d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                           , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                     , Intent( In    ) :: n
    Complex( c_double ), Dimension( : ), Pointer, Intent(   Out ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: rank = 1
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_1d_size_in_int

  Subroutine seg_create_complex_2d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                              , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                        , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, : ), Pointer, Intent(   Out ) :: a
    Integer                                        , Intent(   Out ) :: error

    Integer, Parameter :: rank = 2
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_2d_size_in_int

  Subroutine seg_create_complex_3d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                 , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                           , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: rank = 3
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_3d_size_in_int

  Subroutine seg_create_complex_4d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                    , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                              , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                              , Intent(   Out ) :: error

    Integer, Parameter :: rank = 4
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_4d_size_in_int

  Subroutine seg_create_complex_5d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                       , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                                 , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                                 , Intent(   Out ) :: error

    Integer, Parameter :: rank = 5
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_5d_size_in_int

  Subroutine seg_create_complex_6d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                          , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                                    , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                                    , Intent(   Out ) :: error

    Integer, Parameter :: rank = 6
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_6d_size_in_int

  Subroutine seg_create_complex_7d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                             , Intent( InOut    ) :: ctxt
    Integer, Dimension( : )                                       , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                                       , Intent(   Out ) :: error

    Integer, Parameter :: rank = 7
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.tmpl'

  End Subroutine seg_create_complex_7d_size_in_int

  ! Now follow the segment freeing

  Subroutine seg_free_integer_1d( a, error )

    Integer, Dimension( : ), Pointer, Intent( InOut ) :: a
    Integer                         , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_1d

  Subroutine seg_free_integer_2d( a, error )

    Integer, Dimension( :, : ), Pointer, Intent( InOut ) :: a
    Integer                            , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_2d

  Subroutine seg_free_integer_3d( a, error )

    Integer, Dimension( :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                               , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_3d

  Subroutine seg_free_integer_4d( a, error )

    Integer, Dimension( :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                  , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_4d

  Subroutine seg_free_integer_5d( a, error )

    Integer, Dimension( :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_5d

  Subroutine seg_free_integer_6d( a, error )

    Integer, Dimension( :, :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                        , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_6d

  Subroutine seg_free_integer_7d( a, error )

    Integer, Dimension( :, :, :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, :, :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_integer_7d

  Subroutine seg_free_double_1d( a, error )

    Real( c_double ), Dimension( : ), Pointer, Intent( InOut ) :: a
    Integer                                  , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_1d

  Subroutine seg_free_double_2d( a, error )

    Real( c_double ), Dimension( :, : ), Pointer, Intent( InOut ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_2d

  Subroutine seg_free_double_3d( a, error )

    Real( c_double ), Dimension( :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                        , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_3d

  Subroutine seg_free_double_4d( a, error )

    Real( c_double ), Dimension( :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_4d

  Subroutine seg_free_double_5d( a, error )

    Real( c_double ), Dimension( :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                              , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_5d

  Subroutine seg_free_double_6d( a, error )

    Real( c_double ), Dimension( :, :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                                 , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_6d

  Subroutine seg_free_double_7d( a, error )

    Real( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                                    , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_double_7d

  Subroutine seg_free_complex_1d( a, error )

    Complex( c_double ), Dimension( : ), Pointer, Intent( InOut ) :: a
    Integer                                     , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_1d

  Subroutine seg_free_complex_2d( a, error )

    Complex( c_double ), Dimension( :, : ), Pointer, Intent( InOut ) :: a
    Integer                                        , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_2d

  Subroutine seg_free_complex_3d( a, error )

    Complex( c_double ), Dimension( :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_3d

  Subroutine seg_free_complex_4d( a, error )

    Complex( c_double ), Dimension( :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                              , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_4d

  Subroutine seg_free_complex_5d( a, error )

    Complex( c_double ), Dimension( :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                                 , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_5d

  Subroutine seg_free_complex_6d( a, error )

    Complex( c_double ), Dimension( :, :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                                    , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_6d

  Subroutine seg_free_complex_7d( a, error )

    Complex( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                                       , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.tmpl'

  End Subroutine seg_free_complex_7d

!Nitya - uncommented the routine
!cannot guarantee it is bug free!

  Subroutine allreduce_double( a, n, op, ctxt, error )

    Real( c_double ), Dimension( * ), Intent( InOut ) :: a
    Integer                         , Intent( In    ) :: n
    Integer                         , Intent( In    ) :: op
    Type( FIPC_ctxt )               , Intent( In    ) :: ctxt
    Integer                         , Intent(   Out ) :: error

    Real( c_double ), Dimension( : ), Allocatable :: buffer

    Integer :: buff_size
    Integer :: begin, finish

    Allocate( buffer( 1:Min( max_buff_size, n ) ), Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If

    buff_size = Size( buffer )

    begin = 1
    Do While( begin <= n )

       finish = Min( begin + buff_size - 1, n )

!Nitya ANOTHER STATEMENT TO MAKE THIS WORK
!the mpi_allreduce has to be done only on those procs that have rank 0
!in the intra_comm

      if(ctxt%intra_comm%rank == 0) then 

        Call mpi_allreduce( a( begin:finish ), buffer, finish - begin + 1, &
             MPI_DOUBLE_PRECISION, op, ctxt%extra_comm%handle, error )
        If( error /= 0 ) Then
           Return
        End If

      end if 

!copy contents of buffer to a, only rank 0 has to do it
      if(ctxt%world_comm%rank == 0) then 
        a( begin:finish ) = buffer
      end if 

!END OF CHANGE

       begin = finish + 1
    End Do

    Call FIPC_node_barrier( ctxt, error )
    If( error /= 0 ) Then
       Return
    End If

    Deallocate( buffer )

    error = FIPC_success
   
  End Subroutine allreduce_double

  Function sizeof_what( what ) Result( r )

    Integer( c_long ) :: r

    Integer, Intent( In ) :: what

    Select Case( what )
    Case( integer_vals )
       r = FIPC_sizeof_c_int()
    Case( double_vals )
       r = FIPC_sizeof_c_double()
    Case( complex_vals )
       r = FIPC_sizeof_c_complex()
    Case Default
       r = -1
    End Select

  End Function sizeof_what

  Function test_eq_ctxt( a, b )

    ! Test two contexts for equality

    Logical                         :: test_eq_ctxt
    Type( FIPC_ctxt ), Intent( In ) :: a
    Type( FIPC_ctxt ), Intent( In ) ::b 

    If( Associated( a%initialized ) .And. Associated( b%initialized ) ) Then
       test_eq_ctxt = a%world_comm == b%world_comm .And. &
                      a%intra_comm == b%intra_comm
       If( a%intra_comm%rank == 0 ) Then
          test_eq_ctxt = test_eq_ctxt .And. a%extra_comm == b%extra_comm
       End If
    Else
       test_eq_ctxt = .False.
    End If

  End Function test_eq_ctxt

  Function test_ne_ctxt( a, b )

    ! Test two contexts for inequality

    Logical                    :: test_ne_ctxt
    Type( FIPC_ctxt ), Intent( In ) :: a
    Type( FIPC_ctxt ), Intent( In ) :: b

    test_ne_ctxt = .Not. a == b

  End Function test_ne_ctxt

  Function test_eq_comm( a, b )

    ! Test two communicators for equality

    Logical                            :: test_eq_comm
    Type( communicator ), Intent( In ) :: a
    Type( communicator ), Intent( In ) :: b

    If( a%initialized .And. b%initialized ) Then
       test_eq_comm = a%handle == b%handle .And. &
                      a%size   == b%size   .And. &
                      a%rank   == b%rank
    Else
       test_eq_comm = .False.
    End If

  End Function test_eq_comm

  Function test_ne_comm( a, b )

    ! Test two communicators for inequality

    Logical                            :: test_ne_comm
    Type( communicator ), Intent( In ) :: a
    Type( communicator ), Intent( In ) :: b

    test_ne_comm = .Not. a == b

  End Function test_ne_comm

!Nitya - Subroutine to write out the SMS ids to a file 
!In case of abnormal termination, the ids can be read
!from the file and deleted.
  Subroutine FIPC_write_sms( extra_comm, id, type, ctxt)

    implicit none

    Include 'mpif.h'

    Type(Communicator), Pointer, Intent( In    ) :: extra_comm
    Integer( c_int ) , Intent( In    ) :: id
    Integer, Optional, Intent( In    ) :: type
    Type( FIPC_ctxt ), Intent( In ) :: ctxt

    Integer( c_int ) :: ipc_id
    Integer, Dimension(MPI_STATUS_SIZE) :: status1
    Integer :: rank, ierr
    Logical :: exists

       If( extra_comm%rank == 0 ) then

          inquire(file=id_file, exist=exists)

          If(exists) then
             open(unit=6, file=trim(id_file), status="old", action = "write", &
                  position = "append", iostat=ierr)
          Else
             open(unit=6, file=trim(id_file), status="new", action = "write", &
                  position = "rewind", iostat=ierr)
          End If

          !If type is present, then its a semaphore
          !Since the type is passed only while writing the semaphore id, we assume (correctly?)
          !that the other processes would have also called the subroutine in the same context.
          If(present(type)) then 
             write(unit=6, fmt="(a9,i15)") 'ipcrm -s ' , id
          Else 
             write(unit=6, fmt="(a9,i15)") 'ipcrm -m ' , id
          End If

          !Receive ids from all ranks in root_comm.
          Do rank=1, extra_comm%size - 1
             call mpi_recv(ipc_id, 1, MPI_INTEGER, rank, 0, extra_comm%handle, status1, ierr)

             If(present(type)) then
                write(unit=6, fmt="(a9,i15)") 'ipcrm -s ' , ipc_id
             Else
                write(unit=6, fmt="(a9,i15)") 'ipcrm -m ', ipc_id
             End If

          End Do

          close(unit=6, status="keep", iostat=ierr)

       Else

          call mpi_ssend(id, 1, MPI_INTEGER, 0, 0, extra_comm%handle, ierr)   !send ids to rank 0 in root comm

       End If

  End Subroutine FIPC_write_sms

End Module FIPC_module
