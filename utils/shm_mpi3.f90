! This module handles share memory using MPI3 and records the allocation in a linked list
! such that shared memory can be accessed later in arbitrary location
! of the code (hopefully)
!
!
! To do:
! 1) node barrier with mpi_win_fence (?)
! 2) extend the redistribution to David's non-blocking mpi

! Questions:
! Can the reference from application array to shm array done safe in fortran (allowing for 
! copy in). Does target attribute helps? What happens if the size on a node is 0?
! Is base_ptr set to c_nul_ptr?

! More on the above questions
! A clearer way to identify the segments is to return a tag whem shm_alloc is called
! and to use that tag to get the asociated arrays inside the node.
! the oply drawback of this is that the application must keep a record of the tags, but 
! this seems similar to other MPI-ish kind of tags.

! Lucian Anton May 2014


module shm_mpi3
  use, intrinsic :: iso_c_binding, only : c_ptr
  use mpi
implicit none
private

!<LA shared memory public procedures and data types
  public :: shm_init, shm_alloc, shm_free, &
       shm_onnode, shm_node_id, shm_get_node_pointer, &
       shm_node_barrier, shm_clean, shm_fence
       
!LA>
  
  interface shm_alloc
     module procedure shm_alloc_c1
     module procedure shm_alloc_c2
     module procedure shm_alloc_c3
     module procedure shm_alloc_r1
     module procedure shm_alloc_r2
     module procedure shm_alloc_r3
  end interface shm_alloc
  
  interface shm_free
     module procedure shm_free_c1
     module procedure shm_free_c2
     module procedure shm_free_c3
     module procedure shm_free_r1
     module procedure shm_free_r2
     module procedure shm_free_r3
  end interface shm_free
  
  interface shm_get_node_pointer
     module procedure shm_get_node_pointer_c1
     module procedure shm_get_node_pointer_c2
     module procedure shm_get_node_pointer_c3
     module procedure shm_get_node_pointer_r1
     module procedure shm_get_node_pointer_r2
     module procedure shm_get_node_pointer_r3
  end interface shm_get_node_pointer

  interface remap_bounds
     module procedure remap_bounds_1c
     module procedure remap_bounds_2c
     module procedure remap_bounds_3c
     module procedure remap_bounds_1r
     module procedure remap_bounds_2r
     module procedure remap_bounds_3r
  end interface remap_bounds

  interface shm_fence
     module procedure shm_fence_c
     module procedure shm_fence_r
  end interface shm_fence

  integer, parameter :: maxlen=127

  type shm_info_t
     integer comm, wcomm, size, id ! comm, number of mpi ranks and id in the node
     integer, allocatable :: wranks(:) ! maps the node ranks to world ranks 
  end type shm_info_t

  type shm_node_pointers_t
     integer id, win, ndim
     type(c_ptr), allocatable :: nd(:)
     type(shm_node_pointers_t), pointer :: next => null()
     integer tag ! counter that can be used to identify a given shared segment
     integer, allocatable :: se(:,:) ! lower upper bounds of the pointer associated in the node
     character(len=maxlen) label ! can be used to store info about application
                                 ! (e.g. in which subroutine the segment was created) 
  end type shm_node_pointers_t
  
  type(shm_node_pointers_t), pointer :: shm_pointers => null(), shm_ptr_head => null()

  type (shm_info_t), save, public :: shm_info
! LA> 
  
  integer, save :: counter = 0 ! segment counter
  integer, save :: info_noncontig = MPI_INFO_NULL
  logical, parameter :: debug=.false.

contains

  subroutine shm_init(wcomm)
    implicit none
    
    integer, intent(in) :: wcomm

    integer comm, id_world, id_node, n, ierr
    integer(kind=MPI_ADDRESS_KIND) ta

    ! test for MPI version 

    call mpi_comm_split_type(wcomm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm, ierr)
    call  mpi_comm_size(comm, n, ierr)
    call mpi_comm_rank(comm, id_node, ierr)
    shm_info%id = id_node
    shm_info%size = n
    shm_info%comm = comm
    shm_info%wcomm = wcomm
    
    allocate(shm_info%wranks(0:n-1))        
    call mpi_comm_rank(wcomm, id_world, ierr)        
    call mpi_allgather(id_world, 1, mpi_integer, &
         shm_info%wranks, 1, mpi_integer, comm, ierr)

    ! use contigous block for accelerated ffts
    !call mpi_info_create(info_noncontig, ierr)
    !call mpi_info_set(info_noncontig, "alloc_shared_noncontig", "true", ierr)

    ! check the size of MPI_ADRESS_KIND
    if (id_world == 0) then 
       write(*,*) "shm_mpi3 init: test MPI_ADDRESS_KIND vs integer ", kind(ta), kind(n)
    endif
       

  end subroutine shm_init

  subroutine shm_alloc_c1(a, lubd, tag, label, ierror)
    integer, intent(in) :: lubd(:) ! upper * lower bound in the following format (s1, e1, s2, e2, ...)
    complex, pointer, intent(inout) :: a(:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    integer ndim
    
    ndim = size(lubd) /2
    
    call shm_alloc_c(ndim, lubd, a1=a, tag=tag, label=label, ierror=ierror)

  end subroutine shm_alloc_c1

  subroutine shm_alloc_c2(a, lubd, tag, label, ierror)
    
    integer, intent(in) :: lubd(:) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    complex, pointer,intent(inout) :: a(:,:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    integer ndim
    
    ndim = size(lubd) /2
    
    call shm_alloc_c(ndim, lubd, a2=a, tag=tag, label=label, ierror=ierror)

  end subroutine shm_alloc_c2

  subroutine shm_alloc_c3(a, lubd, tag, label, ierror)

    integer, intent(in) :: lubd(:) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    complex, pointer,intent(inout) :: a(:,:,:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    integer ndim
    
    ndim = size(lubd) /2
    
    call shm_alloc_c(ndim, lubd, a3=a, tag=tag, label=label, ierror=ierror)
    
  end subroutine shm_alloc_c3

  subroutine shm_alloc_c(ndim, lubd, a1, a2, a3, tag, label, ierror)
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_loc
    implicit none
    ! allocated shared array and stores the information in shm_info
    integer, intent(in) :: ndim ! array rank
    integer, intent(in) :: lubd(2*ndim) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    complex, pointer, optional, intent(inout) :: a1(:), a2(:,:), a3(:,:,:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 

    type(c_ptr) :: base_ptr
    integer(KIND=MPI_ADDRESS_KIND) asize
    integer disp, id, i, win, ierr
    integer, allocatable :: ashape(:)
    complex x

    !sanity checks
    if ( ndim == 1 .and. .not. present(a1) .or. &
         ndim == 2 .and. .not. present(a2) .or. &
         ndim == 3 .and. .not. present(a3) ) then
       call error_abort("inconsistent ndim and optional array argumnent in shm_alloc call")
    end if

    if (present(ierror)) ierror = 0

    !call MPI_type_size(MPI_DOUBLE_COMPLEX,disp,ierr)

    include "shm_mpi3_alloc_tmpl.inc"

  end subroutine shm_alloc_c

  subroutine shm_alloc_r1(a, lubd, tag, label, ierror)
    integer, intent(in) :: lubd(:) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    real, pointer,intent(inout) :: a(:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    integer ndim
    
    ndim = size(lubd) /2
    
    call shm_alloc_r(ndim, lubd, a1=a, tag=tag, label=label, ierror=ierror)

  end subroutine shm_alloc_r1

  subroutine shm_alloc_r2(a, lubd, tag, label, ierror)
    
    integer, intent(in) :: lubd(:) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    real, pointer,intent(inout) :: a(:,:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    integer ndim
    
    ndim = size(lubd) /2
    
    call shm_alloc_r(ndim, lubd, a2=a, tag=tag, label=label, ierror=ierror)

  end subroutine shm_alloc_r2

  subroutine shm_alloc_r3(a, lubd, tag, label, ierror)

    integer, intent(in) :: lubd(:) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    real, pointer,intent(inout) :: a(:,:,:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    integer ndim
    
    ndim = size(lubd) /2
    
    call shm_alloc_r(ndim, lubd, a3=a, tag=tag, label=label, ierror=ierror)
    
  end subroutine shm_alloc_r3

  subroutine shm_alloc_r(ndim, lubd, a1, a2, a3, tag, label, ierror)
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_loc
    implicit none
    ! allocated shared array and stores the information in shm_info
    integer, intent(in) :: ndim ! array rank
    integer, intent(in) :: lubd(2*ndim) ! upper * lowee bound in the following format (s1, e1, s2, e2, ...)
    real, pointer, optional, intent(inout) :: a1(:), a2(:,:), a3(:,:,:)
    integer, optional, intent(out) :: tag
    character(len=maxlen), optional, intent(in) :: label
    integer, intent(out), optional :: ierror 
    
    type(c_ptr) :: base_ptr
    integer(KIND=MPI_ADDRESS_KIND) asize
    integer disp, id, i, win, ierr
    integer, allocatable :: ashape(:)
    real x ! for sizeof 

    !sanity checks
    if ( ndim == 1 .and. .not. present(a1) .or. &
         ndim == 2 .and. .not. present(a2) .or. &
         ndim == 3 .and. .not. present(a3) ) then
       call error_abort("inconsistent ndim and optional array argumnent in shm_alloc call")
    end if

    if (present(ierror)) ierror = 0

    !call MPI_type_size(MPI_DOUBLE_PRECISION,disp,ierr)
   

    include "shm_mpi3_alloc_tmpl.inc"

  end subroutine shm_alloc_r

  subroutine shm_free_c1(a)
    implicit none
    complex, intent(inout) :: a(:) 
    call shm_free_c(a)
  end subroutine shm_free_c1

  subroutine shm_free_c2(a)
    implicit none
    complex, intent(inout) :: a(:,:) 
    call shm_free_c(a)
  end subroutine shm_free_c2

  subroutine shm_free_c3(a)
    implicit none
    complex, intent(inout) :: a(:,:,:) 
    call shm_free_c(a)
  end subroutine shm_free_c3


  subroutine shm_free_c(a)
    use, intrinsic :: iso_c_binding, only : c_loc, c_associated
    implicit none
    complex, target, intent(inout) :: a(*)

    include "shm_mpi3_free_tmpl.inc"
       
  end subroutine shm_free_c


  subroutine shm_free_r1(a)
    implicit none
    real, intent(inout) :: a(:) 
    call shm_free_r(a)
  end subroutine shm_free_r1

  subroutine shm_free_r2(a)
    implicit none
    real, intent(inout) :: a(:,:) 
    call shm_free_r(a)
  end subroutine shm_free_r2

  subroutine shm_free_r3(a)
    implicit none
    real, intent(inout) :: a(:,:,:) 
    call shm_free_r(a)
  end subroutine shm_free_r3

subroutine shm_free_r(a)
    use, intrinsic :: iso_c_binding, only : c_loc, c_associated
    implicit none
    real, target, intent(inout) :: a(*)

    include "shm_mpi3_free_tmpl.inc"
       
  end subroutine shm_free_r


  subroutine shm_node_barrier(tag)
!# ifdef MPI
    implicit none
    integer, optional, intent(in) :: tag ! to be used for MPI_Win_fence
    integer ierr
    call mpi_barrier(shm_info%comm, ierr)
!#endif
  end subroutine shm_node_barrier


 subroutine shm_fence_r(a)
    use, intrinsic :: iso_c_binding, only : c_loc, c_associated
    implicit none
    real, target, intent(in) :: a

    include "shm_mpi3_fence_tmpl.inc"

 end subroutine shm_fence_r


 subroutine shm_fence_c(a)
    use, intrinsic :: iso_c_binding, only : c_loc, c_associated
    implicit none
    complex, target, intent(in) :: a

    include "shm_mpi3_fence_tmpl.inc"

 end subroutine shm_fence_c


!********************* shm critical regions ********
!!$  subroutine shm_node_critical_start
!!$# ifdef MPI
!!$   use FIPC_module
!!$     implicit none
!!$
!!$    integer ierr
!!$    Call FIPC_critical_start(FIPC_ctxt_world, ierr)
!!$#endif
!!$  end subroutine shm_node_critical_start
!!$
!!$
!!$  subroutine shm_node_critical_end
!!$# ifdef MPI
!!$    use FIPC_module
!!$    implicit none
!!$
!!$    integer ierr
!!$    Call FIPC_critical_end(FIPC_ctxt_world, ierr)
!!$#endif
!!$  end subroutine shm_node_critical_end
!!$

!********************* shm auxiliaries **************************

  function shm_onnode(ip)
! cheks if a node rank belong to the curent node
    implicit none
    integer, intent(in) :: ip
    logical shm_onnode

    integer i

    shm_onnode = .false.
    
    do i = 0, shm_info%size - 1
       if (ip == shm_info%wranks(i)) then
          shm_onnode = .true.
          exit
       end if
    enddo
    
  end function shm_onnode

  function shm_node_id(ip)
    ! returns node id corresponding to ip
    ! -1 in case of error
    integer, intent(in) :: ip
    integer shm_node_id
    
    integer i
    
    !write(0,*) 'node_id ', ip, shm_info%size,shm_info%wranks(0:)
    
    shm_node_id = -1
    do i = 0, shm_info%size-1
       if ( ip == shm_info%wranks(i)) then
          shm_node_id = i
          exit
       endif
    end do
  end function shm_node_id

  function shm_get_node_pointer_c1(pin, id, tag) result(ptr)
    implicit none
    complex, target, intent(in) :: pin(:)
    integer, intent(in) :: id ! rank id in node comm
    integer, optional, intent(in) :: tag
    complex, pointer :: ptr(:)

    call shm_get_node_pointer_c(pin, id, a1=ptr, tag=tag)
    
  end function shm_get_node_pointer_c1

  function shm_get_node_pointer_c2(pin, id, tag) result(ptr)
    implicit none
    complex, target, intent(in) :: pin(:,:)
    integer, intent(in) :: id ! rank id in node comm
    integer, optional, intent(in) :: tag
    complex, pointer :: ptr(:,:)
    
    call shm_get_node_pointer_c(pin, id, a2=ptr, tag=tag)
    
  end function shm_get_node_pointer_c2
  
  function shm_get_node_pointer_c3(pin, id, tag) result(ptr)
    implicit none
    complex, target, intent(in) :: pin(:,:,:)
    integer, intent(in) :: id ! rank id in node comm
    integer, optional, intent(in) :: tag
    complex, pointer :: ptr(:,:,:)

    call shm_get_node_pointer_c(pin, id, a3=ptr, tag=tag)
    
  end function shm_get_node_pointer_c3

  
  function shm_get_node_pointer_r1(pin, id, tag) result(ptr)
    implicit none
    real, target, intent(in) :: pin(:)
    integer, intent(in) :: id ! rank id in node comm
    integer, optional, intent(in) :: tag
    real, pointer :: ptr(:)

    call shm_get_node_pointer_r(pin, id, a1=ptr, tag=tag)
    
  end function shm_get_node_pointer_r1

  function shm_get_node_pointer_r2(pin, id, tag) result(ptr)
    implicit none
    real, target, intent(in) :: pin(:,:)
    integer, intent(in) :: id ! rank id in node comm
    integer, optional, intent(in) :: tag
    real, pointer :: ptr(:,:)
    
    call shm_get_node_pointer_r(pin, id, a2=ptr, tag=tag)
    
  end function shm_get_node_pointer_r2
  
  function shm_get_node_pointer_r3(pin, id, tag) result(ptr)
    implicit none
    real, target, intent(in) :: pin(:,:,:)
    integer, intent(in) :: id ! rank id in node comm
    integer, optional, intent(in) :: tag
    real, pointer :: ptr(:,:,:)
    
    call shm_get_node_pointer_r(pin, id, a3=ptr, tag=tag)
    
  end function shm_get_node_pointer_r3


  subroutine shm_get_node_pointer_c(pin, id, a1, a2, a3, tag)
    use, intrinsic :: iso_c_binding, only : c_loc, c_associated, c_f_pointer
    implicit none
    complex, target, intent(in) :: pin(*)
    integer, intent(in) :: id ! rank id in node comm
    complex, pointer, optional, intent(out) :: a1(:)
    complex, pointer, optional, intent(out) :: a2(:,:)
    complex, pointer, optional, intent(out) :: a3(:,:,:)
    integer, optional, intent(in) :: tag
    
    include "shm_mpi3_get_node_pointer_tmpl.inc"

  end subroutine shm_get_node_pointer_c


  subroutine shm_get_node_pointer_r(pin, id, a1, a2, a3, tag)
    use, intrinsic :: iso_c_binding, only : c_loc, c_associated, c_f_pointer
    implicit none
    real, target, intent(in) :: pin(*)
    integer, intent(in) :: id ! rank id in node comm
    real, pointer, optional, intent(out) :: a1(:)
    real, pointer, optional, intent(out) :: a2(:,:)
    real, pointer, optional, intent(out) :: a3(:,:,:)
    integer, optional, intent(in) :: tag
    
    include "shm_mpi3_get_node_pointer_tmpl.inc"

  end subroutine shm_get_node_pointer_r

  subroutine shm_clean()
    ! the linked list used o store shm data
    ! must be called before mpi_finalize or after
    ! all shm business is finished
    ! Probabale mpi finalise does a similar think
    implicit none

    integer ierr
    type (shm_node_pointers_t), pointer :: current => null(), aux => null()

    call shm_node_barrier

    current => shm_ptr_head
    do while (associated(current))
       if ( current%win /= MPI_WIN_NULL) then 
          call mpi_win_free(current%win,ierr)
       endif
       deallocate(current%nd, current%se)
       aux => current
       current => current%next
       deallocate(aux)
    end do
  end subroutine shm_clean

FUNCTION remap_bounds_1c(lb1, array) RESULT(ptr)
  INTEGER, INTENT(IN)                          :: lb1
  complex, DIMENSION(lb1:), INTENT(IN), TARGET :: array
  complex, DIMENSION(:), POINTER               :: ptr
  ptr => array
END FUNCTION remap_bounds_1c

FUNCTION remap_bounds_2c(lb1, lb2, array) RESULT(ptr)
  INTEGER, INTENT(IN)                               :: lb1,lb2
  complex, DIMENSION(lb1:, lb2:), INTENT(IN), TARGET :: array
  complex, DIMENSION(:,:), POINTER                  :: ptr
  ptr => array
END FUNCTION remap_bounds_2c

FUNCTION remap_bounds_3c(lb1, lb2, lb3, array) RESULT(ptr)
  INTEGER, INTENT(IN)                            :: lb1,lb2,lb3
  complex, DIMENSION(lb1:,lb2:,lb3:), INTENT(IN), TARGET :: array
  complex, DIMENSION(:,:,:), POINTER                  :: ptr
  ptr => array
END FUNCTION remap_bounds_3c

FUNCTION remap_bounds_1r(lb1, array) RESULT(ptr)
  INTEGER, INTENT(IN)                          :: lb1
  real, DIMENSION(lb1:), INTENT(IN), TARGET :: array
  real, DIMENSION(:), POINTER               :: ptr
  ptr => array
END FUNCTION remap_bounds_1r

FUNCTION remap_bounds_2r(lb1, lb2, array) RESULT(ptr)
  INTEGER, INTENT(IN)                               :: lb1,lb2
  real, DIMENSION(lb1:, lb2:), INTENT(IN), TARGET :: array
  real, DIMENSION(:,:), POINTER                  :: ptr
  ptr => array
END FUNCTION remap_bounds_2r

FUNCTION remap_bounds_3r(lb1, lb2, lb3, array) RESULT(ptr)
  INTEGER, INTENT(IN)                            :: lb1,lb2,lb3
  real, DIMENSION(lb1:,lb2:,lb3:), INTENT(IN), TARGET :: array
  real, DIMENSION(:,:,:), POINTER                  :: ptr
  ptr => array
END FUNCTION remap_bounds_3r

subroutine error_abort(s)
  implicit none
  character(len=*), intent(in) :: s

  integer r, ierr

  write(0,'(a)') s
  call mpi_comm_rank(shm_info%wcomm, r, ierr)
  write(0,*) 'rank ', r, 'node id ', shm_info%id
  call mpi_abort(shm_info%wcomm, -1, ierr)
end subroutine error_abort

end module shm_mpi3
