! Modifications for optimised local copy in c_redist_22 and c_redist_32 
! (and their inverse routines):
! (c) The Numerical Algorithms Group (NAG) Ltd, 2012
! on behalf of EPSRC for the HECToR project
module redistribute
!
! Redistribute distributed (integer, real, complex or logical) 
! (1, 2, 3, or 4) dimensional arrays into two dimensional arrays with 
! first index on local processor, and vice versa.
!
! The first operation is called 'gather' and the second is called 'scatter.'
!  
! One can also do a 'fill' operation.  This consists of copying 
! values from a (2, 3, or 4) dimensional array of 
! (integer, real, complex, or logical ) values into 
! another array with the same number of dimensions.
!
! One can also do a three index to four index redistribution for complex numbers.
!
  implicit none
  public :: opt_redist_nbk, using_measure_scatter, opt_redist_persist, opt_redist_persist_overlap
  private

  logical :: opt_local_copy = .false.
  logical :: opt_redist_nbk = .false. !Do we use non-blocking sends?
  logical :: opt_redist_persist = .false. !Do we use persistent sends?
  logical :: opt_redist_persist_overlap = .false. !Do we try to do a bit of overlapping with persistent sends?
  logical :: using_measure_scatter=.true. !If true then don't time in gather/scatter directly
  integer :: naky = 0
  integer :: ntgrid = 0
  integer :: nlambda = 0
  integer :: ntheta0 = 0
  integer :: negrid = 0
  integer :: nx = 0
  integer :: xxf_lo_ulim_proc = 0
  integer :: yxf_lo_ulim_proc = 0
  integer :: g_lo_blocksize = 0
  character (len=5) :: layout = ''

  public :: index_list_type, delete_list
  public :: redist_type, delete_redist
! TT>
  public :: report_map_property, measure_gather, measure_scatter
  public :: gather_count, scatter_count, time_redist
! <TT

  public :: init_redist, gather, scatter
  public :: init_fill, fill
  public :: set_redist_character_type
  public :: set_xxf_optimised_variables
  public :: set_yxf_optimised_variables

  interface gather
     module procedure c_redist_22, r_redist_22, i_redist_22, l_redist_22
     module procedure c_redist_32, r_redist_32, i_redist_32, l_redist_32
     module procedure c_redist_42, r_redist_42, i_redist_42, l_redist_42
     module procedure c_redist_23
     module procedure c_redist_34, c_redist_36
! TT>
     module procedure c_redist_33
! <TT
  end interface

  interface scatter
     module procedure c_redist_12,     r_redist_12,     i_redist_12,     l_redist_12
     module procedure c_redist_22_inv, r_redist_22_inv, i_redist_22_inv, l_redist_22_inv
     module procedure c_redist_32_inv, r_redist_32_inv, i_redist_32_inv, l_redist_32_inv
     module procedure c_redist_42_inv, r_redist_42_inv, i_redist_42_inv, l_redist_42_inv
     module procedure c_redist_36_inv
! TT>
     module procedure c_redist_33_inv
! <TT
  end interface

! TT>
  interface measure_gather
     module procedure measure_gather_32, measure_gather_33
     module procedure measure_gather_22
  end interface

  interface measure_scatter
     module procedure measure_scatter_23, measure_scatter_33
     module procedure measure_scatter_22
  end interface

   interface init_redist
      module procedure init_redist_gs2, init_redist_agk
   end interface

  integer :: gather_count=0, scatter_count=0
  real, save :: time_redist(2)=0.
! <TT

  interface fill
     module procedure c_fill_2, c_fill_3, c_fill_4
     module procedure r_fill_2, r_fill_3, r_fill_4
     module procedure i_fill_2, i_fill_3, i_fill_4
     module procedure l_fill_2, l_fill_3, l_fill_4
  end interface

  type :: index_map
     integer :: nn
     integer, dimension (:), pointer :: k => null() 
     integer, dimension (:), pointer :: l => null() 
     integer, dimension (:), pointer :: m => null() 
     integer, dimension (:), pointer :: n => null()
     integer, dimension (:), pointer :: o => null()
     integer, dimension (:), pointer :: p => null()
  end type index_map

  ! TT: want to add map name, from_layout and to_layout
  type :: redist_type
     private
     integer, dimension (6) :: to_low, from_low, to_high, from_high
     type (index_map), dimension (:), pointer :: to   => null()
     type (index_map), dimension (:), pointer :: from => null()
     complex, dimension (:), pointer :: complex_buff  => null()
     real,    dimension (:), pointer :: real_buff     => null()
     integer, dimension (:), pointer :: integer_buff  => null()
     logical, dimension (:), pointer :: logical_buff  => null()
     integer :: nsend, nrecv !Note these swap meaning when doing inverse communication
     complex, dimension(:,:), allocatable :: cbuff_send, cbuff_recv !Should probably label these with From and To instead os Send/Recv
     integer, dimension(:), allocatable :: send_hand, recv_hand, send_inv_hand, recv_inv_hand
     integer, dimension(:), allocatable :: recv_ip, recv_inv_ip
     character (len=3) :: redistname = ""
  end type redist_type

  type :: index_list_type
     integer, dimension (:), pointer :: first  => null()
     integer, dimension (:), pointer :: second => null()
     integer, dimension (:), pointer :: third  => null()
     integer, dimension (:), pointer :: fourth => null()
     integer, dimension (:), pointer :: fifth => null()
     integer, dimension (:), pointer :: sixth => null()
  end type index_list_type


contains

  subroutine set_redist_character_type(r, chartype)

    type (redist_type), intent (inout) :: r
    character (3), intent (in) :: chartype

    r%redistname = chartype

  end subroutine set_redist_character_type

  subroutine set_yxf_optimised_variables(lyxf_lo_ulim_proc)
   
    integer, intent (in) :: lyxf_lo_ulim_proc

    yxf_lo_ulim_proc = lyxf_lo_ulim_proc

  end subroutine set_yxf_optimised_variables

  subroutine set_xxf_optimised_variables(lopt_local_copy, lnaky, lntgrid, lntheta0, lnlambda, lnegrid, lnx, &
       lxxf_lo_ulim_proc, lg_lo_blocksize, llayout)
   
    logical, intent (in) :: lopt_local_copy
    integer, intent (in) :: lnaky, lntgrid, lxxf_lo_ulim_proc, lg_lo_blocksize, lntheta0, lnlambda, lnegrid, lnx
    character (len=5), intent (in) :: llayout

    opt_local_copy = lopt_local_copy
    naky = lnaky
    ntgrid = lntgrid
    xxf_lo_ulim_proc = lxxf_lo_ulim_proc
    ntheta0 = lntheta0
    nlambda = lnlambda
    negrid = lnegrid
    nx = lnx
    g_lo_blocksize = lg_lo_blocksize
    layout = llayout

  end subroutine set_xxf_optimised_variables

  subroutine setup_persistent(r,char,buff_size)
    use mp, only: iproc, nproc, send_init, recv_init
    type (redist_type), intent (inout) :: r
    character(1), intent (in) :: char
    integer, intent(in) :: buff_size
    integer :: ip, count_r, count_s

    !Count messages
    r%nrecv=0
    do ip=0,nproc-1
       if (r%to(ip)%nn>0.and.ip.ne.iproc) r%nrecv=r%nrecv+1
    enddo
    r%nsend=0
    do ip=0,nproc-1
       if(r%from(ip)%nn>0.and.ip.ne.iproc) r%nsend=r%nsend+1
    enddo

    !//Persistent comms -- Note could use these buffers in nonblocking case as well
    if(opt_redist_persist)then
       !Now allocate some buffers
       !Only implemented complex so far
       if(char.eq.'c')then
          if(r%nrecv>0)then
             allocate(r%cbuff_recv(buff_size,r%nrecv))
             allocate(r%recv_hand(r%nrecv))
             r%recv_hand=0
             allocate(r%send_inv_hand(r%nrecv))
             allocate(r%recv_ip(r%nrecv))
          endif
          if(r%nsend>0)then
             allocate(r%cbuff_send(buff_size,r%nsend))
             allocate(r%send_hand(r%nsend))
             r%send_hand=0
             allocate(r%recv_inv_hand(r%nsend))
             allocate(r%recv_inv_ip(r%nsend))
          endif
       endif

       !Setup communications
       count_s=0
       count_r=0
       do ip=0,nproc-1
          if(r%from(ip)%nn>0.and.ip.ne.iproc)then
             count_s=count_s+1
             call send_init(r%cbuff_send(1:r%from(ip)%nn,count_s),ip,iproc,r%send_hand(count_s))
             call recv_init(r%cbuff_send(1:r%from(ip)%nn,count_s),ip,ip,r%recv_inv_hand(count_s))
             r%recv_inv_ip(count_s)=ip
          endif
          if(r%to(ip)%nn>0.and.ip.ne.iproc)then
             count_r=count_r+1
             call recv_init(r%cbuff_recv(1:r%to(ip)%nn,count_r),ip,ip,r%recv_hand(count_r))
             call send_init(r%cbuff_recv(1:r%to(ip)%nn,count_r),ip,iproc,r%send_inv_hand(count_r))
             r%recv_ip(count_r)=ip
          endif
       enddo
    endif

  end subroutine setup_persistent

  subroutine init_redist_gs2 (r, char, to_low, to_high, to_list, &
       from_low, from_high, from_list, ierr)

    use mp, only: iproc, nproc, proc0, mp_abort
    implicit none
    type (redist_type), intent (inout) :: r
    character(1), intent (in) :: char
    type (index_list_type), dimension (0:nproc-1), intent (in) :: to_list, from_list
    integer, dimension(:), intent (in) :: to_low, from_low, to_high, from_high    
    integer :: j, ip, n_to, n_from, buff_size
    integer, optional, intent (out) :: ierr

    allocate (r%to(0:nproc-1), r%from(0:nproc-1))

    if (present(ierr)) ierr = 0
    buff_size = 0

    do j = 1, size(to_high)
       r%to_low(j) = to_low(j)
    enddo

! TT> added possibility of higher dimension
!    r%to_low(1) = 1
!    r%to_low(2) = to_low ! TT: assumes two-dimensionality
    ! Actually, did above work for redist_23 or 34?
    ! It is used to determine the lower bound of the index
    ! when redist receives array
    ! The new description below hangs up if dimension of to_high is
    ! set incorrectly
!    r%to_low(:) = 1
!    r%to_low(size(to_high)) = to_low
!    print *, 'size(to_high) is', size(to_high), iproc
!    print *, 'r%to_low', r%to_low, iproc
! <TT

    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          n_to = size(to_list(ip)%first)
          r%to(ip)%nn = n_to

          allocate (r%to(ip)%k(n_to))
          allocate (r%to(ip)%l(n_to))

          r%to(ip)%k = to_list(ip)%first
          r%to(ip)%l = to_list(ip)%second
! TT>
          if (associated(to_list(ip)%third)) then
             allocate (r%to(ip)%m(n_to))
             r%to(ip)%m = to_list(ip)%third
          end if
! <TT
          if (associated(to_list(ip)%fourth)) then
             allocate (r%to(ip)%n(n_to))
             r%to(ip)%n = to_list(ip)%fourth
          end if

          if (associated(to_list(ip)%fifth)) then
             allocate (r%to(ip)%o(n_to))
             r%to(ip)%o = to_list(ip)%fifth
          end if

          if (associated(to_list(ip)%sixth)) then
             allocate (r%to(ip)%p(n_to))
             r%to(ip)%p = to_list(ip)%sixth
          end if

          if (ip /= iproc) buff_size = max(buff_size, n_to)
       else
          r%to(ip)%nn = 0
       endif
    enddo

    do j = 1, size(from_low)
       r%from_low(j) = from_low(j)
    enddo

    do j = 1, size(from_high)
       r%from_high(j) = from_high(j)
    enddo

    do j = 1, size(to_high)
       r%to_high(j) = to_high(j)
    enddo


    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then

          n_from = size(from_list(ip)%first)
          r%from(ip)%nn = n_from

          allocate (r%from(ip)%k(n_from))
          allocate (r%from(ip)%l(n_from))

          r%from(ip)%k = from_list(ip)%first
          r%from(ip)%l = from_list(ip)%second

          if (associated (from_list(ip)%third)) then
             allocate (r%from(ip)%m(n_from))
             r%from(ip)%m = from_list(ip)%third
          endif

          if (associated (from_list(ip)%fourth)) then
             allocate (r%from(ip)%n(n_from))
             r%from(ip)%n = from_list(ip)%fourth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_from)
       else
          r%from(ip)%nn = 0
       endif
    enddo


    select case (char)
       case ('c') 
          if (buff_size > 0) allocate (r%complex_buff(buff_size))
       case ('r') 
          if (buff_size > 0) allocate (r%real_buff(buff_size))
       case ('i') 
          if (buff_size > 0) allocate (r%integer_buff(buff_size))
       case ('l') 
          if (buff_size > 0) allocate (r%logical_buff(buff_size))
       case default
          if (proc0) then
             write(*,*) 'Type to be redistributed invalid.  Must stop.'
             write(*,*) char
          endif
          call mp_abort('Type to be redistributed invalid.  Must stop.')
    end select

    !Setup buffers etc.
    call setup_persistent(r,char,buff_size)
  end subroutine init_redist_gs2

  subroutine init_redist_agk (r, char, to_low, to_high, to_list, &
       from_low, from_high, from_list, ierr)

    use mp, only: iproc, nproc, proc0, mp_abort
    implicit none
    type (redist_type), intent (inout) :: r
    character(1), intent (in) :: char
    integer, intent (in) :: to_low
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension (0:) :: to_list, from_list
    type (index_list_type), dimension (0:nproc-1), intent (in) :: to_list, from_list
! <TT
    integer, dimension(:), intent (in) :: from_low, to_high, from_high    
    ! TT: to_high, from_high seem never used

    integer :: j, ip, n_to, n_from, buff_size
    integer, optional, intent (out) :: ierr

    allocate (r%to(0:nproc-1), r%from(0:nproc-1))

    if (present(ierr)) ierr = 0
    buff_size = 0


! TT> added possibility of higher dimension
!    r%to_low(1) = 1
!    r%to_low(2) = to_low ! TT: assumes two-dimensionality
    ! Actually, did above work for redist_23 or 34?
    ! It is used to determine the lower bound of the index
    ! when redist receives array
    ! The new description below hangs up if dimension of to_high is
    ! set incorrectly
    r%to_low(:) = 1
    r%to_low(size(to_high)) = to_low
!    print *, 'size(to_high) is', size(to_high), iproc
!    print *, 'r%to_low', r%to_low, iproc
! <TT

    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          n_to = size(to_list(ip)%first)
          r%to(ip)%nn = n_to

          allocate (r%to(ip)%k(n_to))
          allocate (r%to(ip)%l(n_to))

          r%to(ip)%k = to_list(ip)%first
          r%to(ip)%l = to_list(ip)%second
! TT>
          if (associated(to_list(ip)%third)) then
             allocate (r%to(ip)%m(n_to))
             r%to(ip)%m = to_list(ip)%third
          end if
! <TT

          if (ip /= iproc) buff_size = max(buff_size, n_to)
       else
          r%to(ip)%nn = 0
       endif
    enddo

    do j = 1, size(from_low)
       r%from_low(j) = from_low(j)
    enddo

    do j = 1, size(from_high)
       r%from_high(j) = from_high(j)
    enddo

    do j = 1, size(to_high)
       r%to_high(j) = to_high(j)
    enddo


    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then

          n_from = size(from_list(ip)%first)
          r%from(ip)%nn = n_from

          allocate (r%from(ip)%k(n_from))
          allocate (r%from(ip)%l(n_from))

          r%from(ip)%k = from_list(ip)%first
          r%from(ip)%l = from_list(ip)%second

          if (associated (from_list(ip)%third)) then
             allocate (r%from(ip)%m(n_from))
             r%from(ip)%m = from_list(ip)%third
          endif

          if (associated (from_list(ip)%fourth)) then
             allocate (r%from(ip)%n(n_from))
             r%from(ip)%n = from_list(ip)%fourth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_from)
       else
          r%from(ip)%nn = 0
       endif
    enddo

    select case (char)
       case ('c') 
          if (buff_size > 0) allocate (r%complex_buff(buff_size))
       case ('r') 
          if (buff_size > 0) allocate (r%real_buff(buff_size))
       case ('i') 
          if (buff_size > 0) allocate (r%integer_buff(buff_size))
       case ('l') 
          if (buff_size > 0) allocate (r%logical_buff(buff_size))
       case default
          if (proc0) then
             write(*,*) 'Type to be redistributed invalid.  Must stop.'
             write(*,*) char
          endif
          call mp_abort('Type to be redistributed invalid.  Must stop.')
    end select

    !Setup buffers etc.
    call setup_persistent(r,char,buff_size)

  end subroutine init_redist_agk

  subroutine init_fill (f, char, to_low, to_high, to_list, &
       from_low, from_high, from_list, ierr)

    use mp, only: nproc, proc0, iproc, mp_abort
    implicit none
    type (redist_type), intent (out) :: f
    character(1), intent (in) :: char
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension (0:) :: to_list, from_list
    type (index_list_type), dimension (0:nproc-1), intent (in) :: to_list, from_list
! <TT
    integer, dimension(:), intent (in) :: to_low, from_low, to_high, from_high
    integer, optional, intent (out) :: ierr
    
    integer :: j, ip, n_to, n_from, buff_size

    if (present(ierr)) ierr = 0
    do j = 1, size(to_low)
       f%to_low(j) = to_low(j)
    enddo

    do j = 1, size(from_low)
       f%from_low(j) = from_low(j)
    enddo

    allocate (f%to(0:nproc-1), f%from(0:nproc-1))

    buff_size = 0
    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          n_to = size(to_list(ip)%first)
          f%to(ip)%nn = n_to
          allocate (f%to(ip)%k(n_to))
          f%to(ip)%k = to_list(ip)%first

          if (associated (to_list(ip)%second)) then
             allocate (f%to(ip)%l(n_to))
             f%to(ip)%l = to_list(ip)%second
          endif

          if (associated (to_list(ip)%third)) then
             allocate (f%to(ip)%m(n_to))
             f%to(ip)%m = to_list(ip)%third
          endif

          if (associated (to_list(ip)%fourth)) then
             allocate (f%to(ip)%n(n_to))
             f%to(ip)%n = to_list(ip)%fourth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_to)
       else
          f%to(ip)%nn = 0
       endif
    enddo

    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then
          n_from = size(from_list(ip)%first)
          f%from(ip)%nn = n_from
          allocate (f%from(ip)%k(n_from))
          f%from(ip)%k = from_list(ip)%first

          if (associated (from_list(ip)%second)) then
             allocate (f%from(ip)%l(n_from))
             f%from(ip)%l = from_list(ip)%second
          endif

          if (associated (from_list(ip)%third)) then
             allocate (f%from(ip)%m(n_from))
             f%from(ip)%m = from_list(ip)%third
          endif

          if (associated (from_list(ip)%fourth)) then
             allocate (f%from(ip)%n(n_from))
             f%from(ip)%n = from_list(ip)%fourth
          endif

          if (ip /= iproc) buff_size = max(buff_size, n_from)
       else
          f%from(ip)%nn = 0
       endif
    enddo

    select case (char)
       case ('c') 
          if (buff_size > 0) allocate (f%complex_buff(buff_size))          
       case ('r') 
          if (buff_size > 0) allocate (f%real_buff(buff_size))
       case ('i') 
          if (buff_size > 0) allocate (f%integer_buff(buff_size))
       case ('l') 
          if (buff_size > 0) allocate (f%logical_buff(buff_size))
       case default
          if (proc0) then
             write(*,*) 'Type to be redistributed invalid.  Must stop.'
             write(*,*) char
          endif
          call mp_abort('Type to be redistributed invalid.  Must stop.')
    end select

    !Setup buffers etc.
    call setup_persistent(f,char,buff_size)
  end subroutine init_fill

  subroutine delete_redist(r)

    use mp, only: nproc, free_request
    type (redist_type), intent (in out) :: r
    
    integer :: i

    if (associated(r%to)) then
       do i = 0, nproc-1
          if (associated (r%to(i)%k)) deallocate (r%to(i)%k)
          if (associated (r%to(i)%l)) deallocate (r%to(i)%l)
          if (associated (r%to(i)%m)) deallocate (r%to(i)%m)
          if (associated (r%to(i)%n)) deallocate (r%to(i)%n)
       enddo
       deallocate (r%to)
    endif

    if (associated(r%from)) then
       do i = 0, nproc-1
          if (associated (r%from(i)%k)) deallocate (r%from(i)%k)
          if (associated (r%from(i)%l)) deallocate (r%from(i)%l)
          if (associated (r%from(i)%m)) deallocate (r%from(i)%m)
          if (associated (r%from(i)%n)) deallocate (r%from(i)%n)
       enddo
       deallocate (r%from)
    endif

    if (associated (r%complex_buff)) deallocate (r%complex_buff)
    if (associated (r%real_buff))    deallocate (r%real_buff)
    if (associated (r%integer_buff)) deallocate (r%integer_buff)
    if (associated (r%logical_buff)) deallocate (r%logical_buff)

    !Persistent related
    if(allocated(r%recv_hand))then
       call free_request(r%recv_hand)
       deallocate(r%recv_hand)
    endif
    if(allocated(r%send_hand))then
       call free_request(r%send_hand)
       deallocate(r%send_hand)
    endif
    if(allocated(r%recv_inv_hand))then
       call free_request(r%recv_inv_hand)
       deallocate(r%recv_inv_hand)
    endif
    if(allocated(r%send_inv_hand))then
       call free_request(r%send_inv_hand)
       deallocate(r%send_inv_hand)
    endif
    if(allocated(r%recv_ip)) deallocate(r%recv_ip)
    if(allocated(r%recv_inv_ip)) deallocate(r%recv_inv_ip)
    if(allocated(r%cbuff_send)) deallocate(r%cbuff_send)
    if(allocated(r%cbuff_recv)) deallocate(r%cbuff_recv)

  end subroutine delete_redist

  subroutine delete_list(list)
    use mp, only: nproc
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension(0:) :: list
    type (index_list_type), dimension(0:nproc-1), intent (inout) :: list
! <TT

    integer :: ip

    do ip = 0, nproc-1
       if(associated(list(ip)%first)) deallocate(list(ip)%first)
       if(associated(list(ip)%second)) deallocate(list(ip)%second)
       if(associated(list(ip)%third)) deallocate(list(ip)%third)
       if(associated(list(ip)%fourth)) deallocate(list(ip)%fourth)
       if(associated(list(ip)%fifth)) deallocate(list(ip)%fifth)
       if(associated(list(ip)%sixth)) deallocate(list(ip)%sixth)
    enddo

  end subroutine delete_list

  subroutine c_redist_12 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):), intent (in) :: from_here

    complex, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_12

  subroutine c_redist_22 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    
    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_22_mpi_copy_persist_start(r, from_here)

    ! redistribute from local processor to local processor
    ! The flag opt_local_copy is set by the user in the input 
    ! file to specify whether optimised local copy routines are used.
    ! These c_redist_22_new_copy routine is the new local copy 
    ! functionality where indirect addressing has largely been removed.
    ! c_redist_22_old_copy is the original local copy code.
    if(opt_local_copy .and. (r%redistname .eq. 'x2y')) then
       ! c_redist_22_new_copy is the new local copy functionality where 
       ! indirect addressing has largely been removed
       call c_redist_22_new_copy(r, from_here, to_here)
    else
       ! c_redist_22_old_copy is the original local copy functionality
       call c_redist_22_old_copy(r, from_here, to_here)
    end if
    
    ! c_redist_22_mpi_copy contains all the remote to local 
    ! copy functionality
    !<DD> Pick from non-blocking and blocking approach
    if(opt_redist_nbk)then
       if(opt_redist_persist)then
          if(.not.opt_redist_persist_overlap) call c_redist_22_mpi_copy_persist_start(r, from_here)
          call c_redist_22_mpi_copy_persist_end(r, to_here)
       else
          call c_redist_22_mpi_copy_nonblock(r, from_here, to_here)
       endif
    else
       call c_redist_22_mpi_copy(r, from_here, to_here)
    endif

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_22

  subroutine c_redist_22_old_copy (r, from_here, to_here)

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the GS2 standard FFT situation this routine maps 
!         xxf(it,ixxf) to yxf(ik,iyxf) data type 
!         where it is kx (or x) index, ik is ky (or y) index, 
!         ixxf is (y,ig,isgn,"les") and iyxf is (x,ig,isgn,"les") 
!
    do i = 1, r%from(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%from(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%to(iproc)%k(i),r%to(iproc)%l(i)) 
!    should come from (r%from(iproc)%k(i), r%from(iproc)%l(i)). 
!
! This do loop, in GS2 standard FFT situation, corresponds to:
!    to_here(ik,iyxf)=from_here(it,ixxf)
!
       to_here(r%to(iproc)%k(i), &
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

  end subroutine c_redist_22_old_copy

  subroutine c_redist_22_new_copy (r, from_here, to_here)
!=====================================================================
!AJ, June 2011: New code from DCSE project on GS2 Indirect Addressing
!=====================================================================
!
! AJ, June 2011:
!  Modified LOCAL COPY part of c_redist_22 routine, as used by GS2 to 
!  transform xxf data type to yxf data type.
!  Here we REDUCE indirect addressing and cache load by EXPLOITING 
!  understanding of xxf and yxf data types in GS2.
! 
    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i,ik,it,itmin,itmax,it_nlocal,ixxf,iyxf

    i = 1
!AJ Loop over all local copies from THIS proc (iproc) to THIS proc
    do while (i .le. r%from(iproc)%nn)
!AJ Initialise the inner loop indices from the existing indirect addresses
       itmin = r%from(iproc)%k(i)
       ixxf = r%from(iproc)%l(i)
       ik = r%to(iproc)%k(i)
       iyxf = r%to(iproc)%l(i)
!AJ it_nlocal is max #it-indices that can be accommodated on iproc
       it_nlocal = (yxf_lo_ulim_proc+1) - iyxf
!AJ itmax selects either the it_nlocal computed above or the maximum
!AJ yxf_lo%nx space available to this process.
       itmax = min((itmin-1)+it_nlocal,nx)
       do it = itmin,itmax
          to_here(ik,iyxf) = from_here(it,ixxf)
          iyxf = iyxf + 1
          i = i + 1
       end do
    end do

  end subroutine c_redist_22_new_copy

  subroutine c_redist_22_mpi_copy (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_22_mpi_copy

  !<DD>This routine does redistributes using entirely non-blocking
  !comms.
  !NOTE: This routine uses more memory than blocking approach
  !--> Can perform slower than blocking routine when redistribute
  !is particularly large, else typically gives better performance
  subroutine c_redist_22_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    use mp, only: barrier,waitall,nbsend,nbrecv

    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here
    complex, dimension(:,:),allocatable :: r_buff,s_buff
    integer,dimension(:),allocatable :: recvhand, sendhand,recvlookup
    integer :: i, ip, nrecv,count,nsend,j

    !Work out how big recv buffer has to be and allocate it
    nrecv=0 !Could really store this in the redist_object
    do i=0,nproc-1
       if (r%to(i)%nn>0.and.i.ne.iproc) nrecv=nrecv+1
    enddo
    if (nrecv>0) then
       allocate(r_buff(SIZE(r%complex_buff),nrecv))
       allocate(recvhand(nrecv))
       allocate(recvlookup(nrecv))
    endif

    !Now work out how big send buffer has to be, allocate it and fill it
    nsend=0 !Could store this in the redist_object
    do i=0,nproc-1
       if(r%from(i)%nn>0.and.i.ne.iproc) nsend=nsend+1
    enddo
    if (nsend>0) then
       allocate(s_buff(SIZE(r%complex_buff),nsend))
       allocate(sendhand(nsend))
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             !Any benefit in making this a pointer assignment instead?
             !Avoids copying data here but could lead to questionable
             !memory access during the actual send?
             do j=1,r%from(i)%nn
                s_buff(j,count)=from_here(r%from(i)%k(j), &
                                              r%from(i)%l(j))
             enddo
             call nbsend(s_buff(1:r%from(i)%nn,count),r%from(i)%nn,i,iproc,sendhand(count))
          endif
       enddo
    endif

    !<DD>Moved after sends
    !Post receives
    !Probably not much advantage to doing these ASAP as not going to do anything
    !until data is on its way --> Could post sends first so chance of data being
    !available as soon as receive posted.
    if (nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             recvlookup(count)=i
             call nbrecv(r_buff(1:r%to(i)%nn,count),r%to(i)%nn,i,i,recvhand(count))
          endif
       enddo
    endif

    if (nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(nrecv,recvhand)

       !Now unpack data
       do i=1,nrecv
          !Which processor is this data from? Would we be better looping over all procs and just checking
          !if that's one to unpack? Removes memory lookup of array but requires branching and incrementing
          !a counter to work out the equivalent of i here (i.e. which receive we're looking at)         
          ip=recvlookup(i)

          !Could this be vectorised and is it useful to do so?
          !Could we make r_buff a pointer array so that MPI put data directly into
          !to_here array? Would this be beneficial? Would need to be done before RECVs posted
          do j=1,r%to(ip)%nn
             to_here(r%to(ip)%k(j), &
                  r%to(ip)%l(j)) &
                  = r_buff(j,i)
          enddo
       enddo

       deallocate(r_bufF)
       deallocate(recvhand)
       deallocate(recvlookup)
    endif

    if (nsend>0) then
       !Shouldn't really need this as we know recvs have worked fine
       call waitall(nsend,sendhand)

       deallocate(s_buff)
       deallocate(sendhand)
    endif

  end subroutine c_redist_22_mpi_copy_nonblock

  subroutine c_redist_22_mpi_copy_persist_start(r, from_here)
    use mp, only: iproc, nproc
    use mp, only: start_comm
    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here
    integer :: i, count,j

    !Start recv comms
    if (r%nrecv>0) then
       call start_comm(r%recv_hand)
    endif

    !Pack buffer and start communication
    if (r%nsend>0) then
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             !Any benefit in making this a pointer assignment instead?
             !Avoids copying data here but could lead to questionable
             !memory access during the actual send?
             do j=1,r%from(i)%nn
                r%cbuff_send(j,count)=from_here(r%from(i)%k(j), &
                                              r%from(i)%l(j))
             enddo
             call start_comm(r%send_hand(count))
          endif
       enddo
    endif
  end subroutine c_redist_22_mpi_copy_persist_start

  subroutine c_redist_22_mpi_copy_persist_end(r, to_here)
    use mp, only: waitall
    implicit none
    type (redist_type), intent (in out) :: r
    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here
    integer :: i, ip, j

    !Wait for recv communications to complete and unpack buffer
    if (r%nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(r%nrecv,r%recv_hand)

       !Now unpack data
       do i=1,r%nrecv
          !Which processor is this data from?
          ip=r%recv_ip(i)

          do j=1,r%to(ip)%nn
             to_here(r%to(ip)%k(j), &
                  r%to(ip)%l(j)) &
                  = r%cbuff_recv(j,i)
          enddo
       enddo
    endif

    !Wait for all sends to complete
    if (r%nsend>0) then
       call waitall(r%nsend,r%send_hand)
    endif
  end subroutine c_redist_22_mpi_copy_persist_end
  !</DD>

  subroutine c_redist_22_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_22_inv_mpi_copy_persist_start(r, from_here)
    ! redistribute from local processor to local processor
    if(opt_local_copy .and. (r%redistname .eq. 'x2y')) then
       ! c_redist_22_inv_new_copy is the new local copy functionality where 
       ! indirect addressing has largely been removed
       call c_redist_22_inv_new_copy(r, from_here, to_here)
    else
       ! c_redist_22_inv_old_copy is the original local copy functionality
       call c_redist_22_inv_old_copy(r, from_here, to_here)
    end if

    ! c_redist_22_inv_mpi_copy contains all the remote to local 
    ! copy functionality
    !<DD>
    if(opt_redist_nbk)then
       if(opt_redist_persist)then
          if(.not.opt_redist_persist_overlap) call c_redist_22_inv_mpi_copy_persist_start(r, from_here)
          call c_redist_22_inv_mpi_copy_persist_end(r, to_here)
       else
          call c_redist_22_inv_mpi_copy_nonblock(r, from_here, to_here)
       endif
    else
       call c_redist_22_inv_mpi_copy(r, from_here, to_here)
    endif

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_22_inv

  subroutine c_redist_22_inv_old_copy (r, from_here, to_here)

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the GS2 standard FFT situation this routine maps 
!         yxf(ik,iyxf) to xxf(it,ixxf) data type 
!         where it is kx (or x) index, ik is ky (or y) index, 
!         ixxf is (y,ig,isgn,"les") and iyxf is (x,ig,isgn,"les") 
!
    do i = 1, r%to(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%from(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%from(iproc)%k(i),r%from(iproc)%l(i)) 
!    should come from (r%to(iproc)%k(i), r%to(iproc)%l(i)). 
!
! This do loop, in GS2 standard FFT situation, corresponds to:
!    to_here(it,ixxf)=from_here(ik,iyxf)
!
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

  end subroutine c_redist_22_inv_old_copy


  subroutine c_redist_22_inv_new_copy (r, from_here, to_here)
!=====================================================================
!AJ, June 2011: New code from DCSE project on GS2 Indirect Addressing
!=====================================================================
!
! AJ, June 2011:
!  Modified LOCAL COPY part of c_redist_22_inv routine, used by GS2 to 
!  transform yxf data type back to xxf data type.
!  Here we REDUCE indirect addressing and cache load by EXPLOITING 
!  understanding of xxf and yxf data types in GS2.
! 
    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i,ik,it,itmin,itmax,it_nlocal,ixxf,iyxf

    i = 1
!AJ Loop over all local copies from THIS proc (iproc) to THIS proc
    do while (i .le. r%to(iproc)%nn)
!AJ Initialise the inner loop indices from the existing indirect addresses
       itmin = r%from(iproc)%k(i)
       ixxf = r%from(iproc)%l(i)
       ik = r%to(iproc)%k(i)
       iyxf = r%to(iproc)%l(i)
!AJ it_nlocal is max #it-indices that can be accommodated on iproc
       it_nlocal = (yxf_lo_ulim_proc+1) - iyxf
!AJ itmax selects either the it_nlocal computed above or the maximum
!AJ yxf_lo%nx space available to this process.
       itmax = min((itmin-1)+it_nlocal,nx)
       do it = itmin,itmax
          to_here(it,ixxf) = from_here(ik,iyxf)
          iyxf = iyxf + 1
          i = i + 1
       end do
    end do

  end subroutine c_redist_22_inv_new_copy


  subroutine c_redist_22_inv_mpi_copy (r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here


    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_22_inv_mpi_copy

  subroutine c_redist_22_inv_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    use mp, only: barrier,waitall,nbsend,nbrecv
    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here
    complex, dimension(:,:),allocatable :: r_buff,s_buff
    integer,dimension(:),allocatable :: recvhand, sendhand,recvlookup
    integer :: i, ip, nrecv,count,nsend,j

    !Work out how big recv buffer has to be and allocate it
    nrecv=0 !Could really store this in the redist_object
    do i=0,nproc-1
       if (r%from(i)%nn>0.and.i.ne.iproc) nrecv=nrecv+1
    enddo
    if (nrecv>0) then
       allocate(r_buff(SIZE(r%complex_buff),nrecv))
       allocate(recvhand(nrecv))
       allocate(recvlookup(nrecv))
    endif

    !Now work out how big send buffer has to be, allocate it and fill it
    nsend=0
    do i=0,nproc-1
       if(r%to(i)%nn>0.and.i.ne.iproc) nsend=nsend+1
    enddo
    if (nsend>0) then
       allocate(s_buff(SIZE(r%complex_buff),nsend))
       allocate(sendhand(nsend))
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             do j=1,r%to(i)%nn
                s_buff(j,count)=from_here(r%to(i)%k(j), &
                                              r%to(i)%l(j))
             enddo
          endif
       enddo
    endif

    !Post receives
    if (nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             recvlookup(count)=i
             call nbrecv(r_buff(1:r%from(i)%nn,count),r%from(i)%nn,i,i,recvhand(count))
          endif
       enddo
    endif

    !Now send data where required
    count=0
    do i=0,nproc-1
       if(r%to(i)%nn>0.and.i.ne.iproc) then
          !Note we could think about using RSEND here as we should post all receives first
          !IN practice this will require syncing so probaly no good
          count=count+1
          call nbsend(s_buff(1:r%to(i)%nn,count),r%to(i)%nn,i,iproc,sendhand(count))
       endif
    enddo

    if (nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(nrecv,recvhand)

       !Now unpack data
       do i=1,nrecv
          !Which processor is this data from
          ip=recvlookup(i)
          do j=1,r%from(ip)%nn
             to_here(r%from(ip)%k(j), &
                  r%from(ip)%l(j)) &
                  = r_buff(j,i)
          enddo
       enddo

       deallocate(r_bufF)
       deallocate(recvhand)
       deallocate(recvlookup)
    endif

    if (nsend>0) then
       !Shouldn't really need this as we know recvs have worked fine
       call waitall(nsend,sendhand)

       deallocate(s_buff)
       deallocate(sendhand)
    endif

  end subroutine c_redist_22_inv_mpi_copy_nonblock

  subroutine c_redist_22_inv_mpi_copy_persist_start(r, from_here)

    use mp, only: iproc, nproc
    use mp, only: start_comm
    implicit none
    
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here
    integer :: i, count,j
    !For inverse transform SEND actually refers to received data and
    !RECV refers to data to be sent

    !Post receives
    if (r%nsend>0) then
       call start_comm(r%recv_inv_hand)
    endif

    !Pack send buffer and start sends
    if (r%nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             do j=1,r%to(i)%nn
                r%cbuff_recv(j,count)=from_here(r%to(i)%k(j), &
                     r%to(i)%l(j))
             enddo
             call start_comm(r%send_inv_hand(count))
          endif
       enddo
    endif
  end subroutine c_redist_22_inv_mpi_copy_persist_start

  subroutine c_redist_22_inv_mpi_copy_persist_end(r, to_here)
    use mp, only: waitall
    implicit none
    
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here
    integer :: i, ip, j
    !For inverse transform SEND actually refers to received data and
    !RECV refers to data to be sent

    !Wait for receives to complete and unpack data
    if (r%nsend>0) then
       call waitall(r%nsend,r%recv_inv_hand)

       !Now unpack data
       do i=1,r%nsend
          !Which processor is this data from
          ip=r%recv_inv_ip(i)
          do j=1,r%from(ip)%nn
             to_here(r%from(ip)%k(j), &
                  r%from(ip)%l(j)) &
                  = r%cbuff_send(j,i)
          enddo
       enddo
    endif

    !Wait for sends to complete
    if (r%nrecv>0) then
       call waitall(r%nrecv,r%send_inv_hand)
    endif
  end subroutine c_redist_22_inv_mpi_copy_persist_end
  !</DD>

  subroutine c_redist_32 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer, save :: optimised_choice = 0
    real :: time_optimised_loop_1(2), time_optimised_loop_2(2)
    integer :: i
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_32_mpi_copy_persist_start(r, from_here)
!AJ redistribute from local processor to local processor
!AJ The flag opt_local_copy is set by the user in the input 
!AJ file to specify whether optimised local copy routines are used.
!AJ These c_redist_32_*_copy routines are the new local copy 
!AJ functionality where indirect addressing has largely been removed
!AJ There are two different versions of the optimised local copy 
!AJ functionality which are selected at runtime (providing 
!AJ opt_local_copy is true) through a performance measurement process 
!AJ documented below.
    if(opt_local_copy .and. (r%redistname .eq.  'g2x')) then
!AJ Because there are two different optimised local copy routines 
!AJ which provide benefits for different process counts and use cases 
!AJ we use an auto-tuning method to select which routine to use.  This 
!AJ works by timing both of the optimised routines on the first run of 
!AJ this functionality and then using the quickest for the rests of the 
!AJ times this routine is called.  The optimised_choice varible is used 
!AJ to record the choice of routine (it is a SAVE variable so will maintain 
!AJ a value between calls to the routine), it is initialise to 0 when the 
!AJ code first runs.  We run each routine 4 times, only timing the last 
!AJ three runs to avoid any potential initialisation penalties.  We time 
!AJ three rather than one to deal ensure that we collect enough timing 
!AJ data.
       if(optimised_choice .eq. 0) then

          call c_redist_32_new_copy(r, from_here, to_here)
          call time_message(.false.,time_optimised_loop_1,' Optimised Loop 1')
          do i= 1,3 
             call c_redist_32_new_copy(r, from_here, to_here)
          end do
          call time_message(.false.,time_optimised_loop_1,' Optimised Loop 1')

          call c_redist_32_new_opt_copy(r, from_here, to_here)
          call time_message(.false.,time_optimised_loop_2,' Optimised Loop 2')
          do i = 1,3
             call c_redist_32_new_opt_copy(r, from_here, to_here)
          end do
          call time_message(.false.,time_optimised_loop_2,' Optimised Loop 2')
          
          if(time_optimised_loop_1(1) .gt. time_optimised_loop_2(1)) then
             optimised_choice = 2
          else
             optimised_choice = 1
          end if

!AJ This else is encountered once the optimised auto-tuning choice has 
!AJ been calculated above.
       else

          if(optimised_choice .eq. 1) then
             call c_redist_32_new_copy(r, from_here, to_here)
          else
             call c_redist_32_new_opt_copy(r, from_here, to_here)
          end if

       end if

!AJ This else is encountered when the optimised local copy functionality is 
!AJ not enabled.
    else
!AJ c_redist_32_old_copy is the original local copy functionality
       call c_redist_32_old_copy(r, from_here, to_here)
    end if
    
!AJ c_redist_32_mpi_copy contains all the remote to local 
!AJ copy functionality
    !<DD>
    if(opt_redist_nbk)then
       if(opt_redist_persist) then
          if(.not.opt_redist_persist_overlap) call c_redist_32_mpi_copy_persist_start(r, from_here) !DD
          call c_redist_32_mpi_copy_persist_end(r, to_here) !DD
       else
          call c_redist_32_mpi_copy_nonblock(r, from_here, to_here) !DD
       endif
    else
       call c_redist_32_mpi_copy(r, from_here, to_here)
    endif

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_32


  subroutine c_redist_32_old_copy(r, from_here, to_here)

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the GS2 standard FFT situation this routine maps 
!         g(ig, isgn, iglo) to xxf(it,ixxf) data type 
!         where it is kx (or x) index, ixxf is (y,ig,isgn,"les") 
!         and iglo is ("xyles") 
    do i = 1, r%from(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%from(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%to(iproc)%k(i),r%to(iproc)%l(i)) 
!    should come from (r%from(iproc)%k(i),r%from(iproc)%l(i),r%from(iproc)%m(i)). 
!
! This do loop, in GS2 standard FFT situation, corresponds to:
!    to_here(it,ixxf)=from_here(ig,isgn,iglo)
!
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

  end subroutine c_redist_32_old_copy


  subroutine c_redist_32_new_copy(r, from_here, to_here)
!=====================================================================
!AJ, June 2011: New code from DCSE project on GS2 Indirect Addressing
!=====================================================================
!
! AJ, June 2011:
!  Modified LOCAL COPY part of c_redist_32 routine, as used by GS2 to 
!  transform g data type to xxf data type.
!  Here we REDUCE indirect addressing and cache load by EXPLOITING 
!  understanding of g and xxf data types in GS2.
!
    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i,k,t2,t1,f3,f2,f1,fhigh,thigh,f2max
    real :: nakyrecip

    i = 1
!AJ We want to be able to divide by naky with floating point 
!AJ arithmetic so we need to convert naky to a float (hence
!AJ the need for nakyrecip).  Also, to remove the necessity 
!AJ to have a division in the body of the loop below we are 
!AJ calculating the reciprocal of naky so that instead of 
!AJ dividing by naky we can multiply  by 1/naky.
    nakyrecip = naky
    nakyrecip = 1/nakyrecip
    f2max = r%from_high(2)
!AJ Loop over all local copies from THIS proc (iproc) to THIS proc
    do while(i .le. r%from(iproc)%nn)
!AJ Initialise the inner loop indices.  As the from_here array has 
!AJ 3 indices (c_redist_32 goes from 3 indices to 2 indices in the
!AJ copy) we need two inner loops to be able to move through all 
!AJ the to_here and from_here indices.
!AJ t1 (which equates to the it index) and f3 (which equates to iglo)
!AJ do not change for the iterations of the inner loops.  f2 (isgn) 
!AJ can only be 1 or 2 so the first inner loop is restricted to at 
!AJ most 2 iterations.
       f2 = r%from(iproc)%l(i)
       f3 = r%from(iproc)%m(i)
       t1 = r%to(iproc)%k(i)
!AJ Ensure that isgn (f2) is either 1 or 2.
       do while (f2 .le. f2max)
!AJ Get initial value of f1 (which equates to ig) and t2 (which equates
!AJ to ixxf).
          f1 = r%from(iproc)%k(i)
          t2 = r%to(iproc)%l(i)
!AJ Work out the maximum value ixxf (t2) can have by calculating the range 
!AJ of ixxf that this process owns.  We step through t2 by naky each iteration
!AJ so the line below calculates thigh as the number of ixxf steps for these 
!AJ inner loops of the operation owned by this process.
          thigh = ceiling(((xxf_lo_ulim_proc+1) - t2)*nakyrecip)
!AJ from_here and to_here are calculated for c_redist_32 using ig, isgn, 
!AJ and iglo.  The calculation of the maximum number if ixxf steps above can, 
!AJ therefore, be used to also calculate the maximum number of ig setsp.  THis 
!AJ is what the following line does.
          thigh = thigh + (f1-1)
!AJ Finally check the actual number of inner loop steps by ensuring that the 
!AJ computed maximum bound of ig is not beyond the actual allowed maximum value 
!AJ in r%from_high(1).
          fhigh = min(thigh,r%from_high(1))
          do k = f1,fhigh
             to_here(t1,t2) = from_here(k,f2,f3)
             t2 = t2 + naky
             i = i + 1
          end do
!AJ If at the end of the inner loop we still have theoretical iterations 
!AJ left on ixxf (i.e. the calculated thigh is higher than the allowed 
!AJ maximum of ig) then move to the next isgn.  If there aren't any 
!AJ iterations left then exit the inner loops (the f2 loop).
          if(thigh .gt. r%from_high(1)) then
             f2 = f2 + 1
          else
             f2 = f2max + 1
          end if
       end do
    end do

  end subroutine c_redist_32_new_copy


  subroutine c_redist_32_new_opt_copy(r, from_here, to_here)
!=====================================================================
!AJ, June 2011: New code from DCSE project on GS2 Indirect Addressing
!=====================================================================
!
! AJ, June 2011:
!  Modified LOCAL COPY part of c_redist_32 routine, as used by GS2 to 
!  transform g data type to xxf data type.
!  Here we REDUCE indirect addressing and cache load by EXPLOITING 
!  understanding of g and xxf data types in GS2.
!  This functionality is different to c_redist_32_new_copy because 
!  that routine does not always provide optimsed performance 
!  (due to the way memory reads and writes are performed, particularly 
!  the use of a strided write which triggers a write allocate cache 
!  miss).
!  This new functionality moves through the loop indices in a different 
!  order to the other optimised local copy, so instead of following the 
!  original functionality by moving through the to_here and from_here loops
!  using the i variable from i = 1 to from(iproc)%nn we start at i = 1 and 
!  then skip through i's for a certain range, then move back to i = 2 and 
!  follow the same process until all the i range has been processed.
!
    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i,t2,t1,f3,f2,f1
    integer :: f3max,f3maxmultiple,f3incr,innermax,iincrem,iglomax,t1test
    integer :: innermaxmultiplier,outerf3limit,startf3
    real :: innermaxrealvalue,tempnaky

!AJ tempnaky is a real version (rather than integer version) of naky
!AJ It is used to enable real arthimetic rather than integer arthimetic.
    tempnaky = naky
!AJ innermaxmultiplier is the upper limit of the xxf data range owned by this process
    innermaxmultiplier = xxf_lo_ulim_proc+1

!AJ This optimised new local copy functionality is reliant on the particular 
!AJ layout used to work out how to increment through the f3 (iglo) indices.
!AJ f3maxmultiple sets the upper bound of the f3 blocks we step through.
!AJ f3incr sets the amount added to f3 at each step of the inner loops.
    select case (layout)
    case ('yxels')
       f3maxmultiple = naky*ntheta0
       f3incr = naky
    case ('yxles')
       f3maxmultiple = naky*ntheta0
       f3incr = naky
    case ('lexys')
       f3maxmultiple = nlambda*negrid*ntheta0
       f3incr = nlambda*negrid
    case ('lxyes')
       f3maxmultiple = nlambda*ntheta0
       f3incr = nlambda
    case ('lyxes')
       f3maxmultiple = nlambda*naky*ntheta0
       f3incr = nlambda*naky
    case('xyles')
       f3maxmultiple = ntheta0
       f3incr = 1
    end select
          

    i = 1
!AJ iglomax is the maximum iglo (f3) owned by this process
    iglomax = ((iproc+1)*g_lo_blocksize)
!AJ t1test is used handle the it (t1) index which has a 
!AJ range 1 -> (xxf_lo%ntheta0+1)/2) and then 
!AJ (it - xxf_lo%ntheta0 + xxf_lo%nx) -> xxf_lo%nx.  t1test is used to 
!AJ identify when this change in range happens and enable the code to 
!AJ compensate for it.
    t1test  = ((ntheta0+1)/2)+1
!AJ Loop over all local copies from THIS proc (iproc) to THIS proc
    do while(i .le. r%from(iproc)%nn)

       
!AJ Initial look up of values needed to calculate innermax
!AJ and setup the first iteration of the loop below (do while i .le. innermax)
       f1 = r%from(iproc)%k(i)
       f2 = r%from(iproc)%l(i)
       f3 = r%from(iproc)%m(i)
       t1 = r%to(iproc)%k(i)
       t2 = r%to(iproc)%l(i)

!AJ outerf3limit calculates the number of f3 iterations that can be undertaken 
!AJ before reaching the next f3 point.
       outerf3limit = f3 + f3incr
       do while(f3 .lt. outerf3limit)
          
!AJ iincrem counts how much i needs to be skipped at the end of the final inner loop
!AJ The i loop below moves through the starting i's for these inner loop iterations
!AJ but the innermost loop increments through further i's without actually moving i
!AJ so once we've finished the loop below we need to move i forward to skip all the 
!AJ elements we have just processed.
          iincrem = i
          
!AJ work out the size of the inner loops body 
!AJ innermax is the range of the t2 index (ixxf) this process 
!AJ owns.
          innermaxrealvalue = (innermaxmultiplier-t2)/tempnaky
          innermax = ceiling(innermaxrealvalue)-1
!AJ If the isgn (f2) index is at maximum (2) and the ig (f1) index 
!AJ will also be at maximum if we use the innermax value calculated above 
!AJ (ig maximum is r%from_high(1) which currently is ntgrid in the code) 
!AJ then limit innermax to the difference between ig and ntgrid (which is 
!AJ ig max).
          if(f2 .eq. 2 .and. (f1+innermax) .gt. r%from_high(1)) then
             innermax = i + (ntgrid - f1)
!AJ If the isgn (f2) index is not at maximum, but the calculated innermax would 
!AJ take the indexes beyond whaty is allowed for ig (which is the maximum of ig 
!AJ multipled by each increment of ig) then restrict innermax to the maximum 
!AJ number of allowed ig increments.
          else if((f2 .eq. 1 .and. innermax .gt. (((2*ntgrid)+1)+(ntgrid-f1)))) then
             innermax = i + ((2*ntgrid)+1) + (ntgrid - f1)
!AJ If we have reached this else then the calculated innermax does not breach 
!AJ any of the data limits it can be used for the number of iterations of the 
!AJ inner loop. 
          else
             innermax = i + innermax
          end if
          
          do while(i .le. innermax) 
             
!AJ Get the initial address variables.  There are potentially some extra data lookups 
!AJ here so there may be scope to optimise this, although this has not been checked.
             f1 = r%from(iproc)%k(i)
             f2 = r%from(iproc)%l(i)
             f3 = r%from(iproc)%m(i)
             t1 = r%to(iproc)%k(i)
             t2 = r%to(iproc)%l(i)

!AJ Record the initial f3 variable value for the inner loops for use in the bounds checking 
!AJ at the end of the inner loops.             
             startf3 = f3
!AJ Calculate the size of the innermost loop by working out the maximum that f3 (iglo) 
!AJ can go to for this process (using the f3maxmultiple as calculated previously that 
!AJ effectively defines the f3 blocksize of this process).
             f3max = ((f3/f3maxmultiple)+1)*f3maxmultiple
!AJ Ensure that the calculate f3max is not larger than the total iglo space this process
!AJ owns.
             f3max = min(f3max,iglomax)
             do while (f3 .lt. f3max)
                to_here(t1,t2) = from_here(f1,f2,f3)                   
                f3 = f3 + f3incr
                t1 = t1 + 1
!AJ Deal with the issues with t1 (it) being a split range variable as described in the 
!AJ comment for t1test.
                if(t1 .eq. t1test) then
                   t1 = t1 - ntheta0 + nx
                end if
                iincrem = iincrem + 1
             end do
             
!AJ Increment the outermost inner loop (move forwards in the i loop one step)
             i = i + 1
                          
          end do

!AJ Ensure that we have not finished the i loop altogether
          if(i .lt. r%from(iproc)%nn .and. iincrem .lt. r%from(iproc)%nn) then
!AJ The f1, f2 and t2 lookups are required to ensure that
!AJ the calculation of innermax the next time round this loop
!AJ are correctly performed as we are now moving forward in the loop
!AJ It may be that it is better to re-order the initialisation of values 
!AJ so this happens at the beginning of the loop rather than at this point.
             f1 = r%from(iproc)%k(i)
             f2 = r%from(iproc)%l(i)
             t2 = r%to(iproc)%l(i)          
!AJ Updated to ensure the f3 loop we are in is correctly setup as with the code 
!AJ above we are moving to a new i so we need to get the correct value of f3 for this i
             f3 =  startf3 + 1
          else
!AJ If we have gone beyond the end of the i loop then force the exit of the inner loops.
!AJ This should never happen as the loop functionality is designed to always perform the 
!AJ correct number of iterations, so this check could be removed from the code.
             f3 = outerf3limit
          end if

       end do

!AJ Once all the inner loops have been iterated move the i variable on to the next point 
!AJ point to be considered.
       i = iincrem
       
    end do
    
  end subroutine c_redist_32_new_opt_copy

  subroutine c_redist_32_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_32_mpi_copy

!<DDHACK>
  subroutine c_redist_32_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    use mp, only: barrier
    use mp, only: waitall, nbsend,nbrecv
    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here
    complex, dimension(:,:),allocatable :: r_buff,s_buff
    integer,dimension(:),allocatable :: recvhand, sendhand,recvlookup
    integer :: i, ip, nrecv,count,nsend,j

    !Work out how big recv buffer has to be and allocate it
    nrecv=0 !Could really store this in the redist_object
    do i=0,nproc-1
       if (r%to(i)%nn>0.and.i.ne.iproc) nrecv=nrecv+1
    enddo
    if (nrecv>0) then
       allocate(r_buff(SIZE(r%complex_buff),nrecv))
       allocate(recvhand(nrecv))
       allocate(recvlookup(nrecv))
    endif

    !Now work out how big send buffer has to be, allocate it and fill it
    nsend=0 !Could store this in the redist_object
    do i=0,nproc-1
       if(r%from(i)%nn>0.and.i.ne.iproc) nsend=nsend+1
    enddo
    if (nsend>0) then
       allocate(s_buff(SIZE(r%complex_buff),nsend))
       allocate(sendhand(nsend))
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             !Any benefit in making this a pointer assignment instead?
             !Avoids copying data here but could lead to questionable
             !memory access during the actual send?
             do j=1,r%from(i)%nn
                s_buff(j,count)=from_here(r%from(i)%k(j), &
                                              r%from(i)%l(j), &
                                              r%from(i)%m(j))
             enddo
             call nbsend(s_buff(1:r%from(i)%nn,count),r%from(i)%nn,i,iproc,sendhand(count))
          endif
       enddo
    endif

!<DD>Moved after sends
    !Post receives
    !Probably not much advantage to doing these ASAP as not going to do anything
    !until data is on its way --> Could post sends first so chance of data being
    !available as soon as receive posted.
    if (nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             recvlookup(count)=i
             call nbrecv(r_buff(1:r%to(i)%nn,count),r%to(i)%nn,i,i,recvhand(count))
          endif
       enddo
    endif

    if (nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(nrecv,recvhand)

       !Now unpack data
       do i=1,nrecv
          !Which processor is this data from? Would we be better looping over all procs and just checking
          !if that's one to unpack? Removes memory lookup of array but requires branching and incrementing
          !a counter to work out the equivalent of i here (i.e. which receive we're looking at)         
          ip=recvlookup(i)

          !Could this be vectorised and is it useful to do so?
          !Could we make r_buff a pointer array so that MPI put data directly into
          !to_here array? Would this be beneficial? Would need to be done before RECVs posted
          do j=1,r%to(ip)%nn
             to_here(r%to(ip)%k(j), &
                  r%to(ip)%l(j)) &
                  = r_buff(j,i)
          enddo
       enddo

       deallocate(r_buff)
       deallocate(recvhand)
       deallocate(recvlookup)
    endif

    if (nsend>0) then
       !Shouldn't really need this as we know recvs have worked fine
       call waitall(nsend,sendhand)
       deallocate(s_buff)
       deallocate(sendhand)
    endif

  end subroutine c_redist_32_mpi_copy_nonblock
  !</DD>

  subroutine c_redist_32_mpi_copy_persist_start(r, from_here)
    use mp, only: iproc, nproc
    use mp, only: start_comm
    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here
    integer :: i, count,j

    !Start recv comms
    if (r%nrecv>0) then
       call start_comm(r%recv_hand)
    endif

    !Pack buffer and start communication
    if (r%nsend>0) then
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             !Any benefit in making this a pointer assignment instead?
             !Avoids copying data here but could lead to questionable
             !memory access during the actual send?
             do j=1,r%from(i)%nn
                r%cbuff_send(j,count)=from_here(r%from(i)%k(j), &
                                              r%from(i)%l(j), &
                                              r%from(i)%m(j))
             enddo
             call start_comm(r%send_hand(count))
          endif
       enddo
    endif
  end subroutine c_redist_32_mpi_copy_persist_start

  subroutine c_redist_32_mpi_copy_persist_end(r, to_here)
    use mp, only: waitall
    implicit none
    type (redist_type), intent (in out) :: r
    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here
    integer :: i, ip, j

    !Wait for recv communications to complete and unpack buffer
    if (r%nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(r%nrecv,r%recv_hand)

       !Now unpack data
       do i=1,r%nrecv
          !Which processor is this data from?
          ip=r%recv_ip(i)

          do j=1,r%to(ip)%nn
             to_here(r%to(ip)%k(j), &
                  r%to(ip)%l(j)) &
                  = r%cbuff_recv(j,i)
          enddo
       enddo
    endif

    !Wait for all sends to complete
    if (r%nsend>0) then
       call waitall(r%nsend,r%send_hand)
    endif
  end subroutine c_redist_32_mpi_copy_persist_end
  !</DD>

  subroutine c_redist_32_inv (r, from_here, to_here)

    use job_manage, only: time_message
    use mp, only: proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer, save :: optimised_choice = 0
    real :: time_optimised_loop_1(2), time_optimised_loop_2(2)
    integer :: i

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_32_inv_mpi_copy_persist_start(r, from_here)

!AJ redistribute from local processor to local processor
!AJ The flag opt_local_copy is set by the user in the input 
!AJ file to specify whether optimised local copy routines are used.
!AJ These c_redist_32_inv_*_copy routines are the new local copy 
!AJ functionality where indirect addressing has largely been removed
!AJ There are two different versions of the optimised local copy 
!AJ functionality which are selected at runtime (providing 
!AJ opt_local_copy is true) through a performance measurement process 
!AJ documented below.
    if(opt_local_copy .and. (r%redistname .eq. 'g2x')) then
!AJ Because there are two different optimised local copy routines 
!AJ which provide benefits for different process counts and use cases 
!AJ we use an auto-tuning method to select which routine to use.  This 
!AJ works by timing both of the optimised routines on the first run of 
!AJ this functionality and then using the quickest for the rests of the 
!AJ times this routine is called.  The optimised_choice varible is used 
!AJ to record the choice of routine (it is a SAVE variable so will maintain 
!AJ a value between calls to the routine), it is initialise to 0 when the 
!AJ code first runs.  We run each routine 4 times, only timing the last 
!AJ three runs to avoid any potential initialisation penalties.  We time 
!AJ three rather than one to deal ensure that we collect enough timing 
!AJ data.
       if(optimised_choice .eq. 0) then

          call c_redist_32_inv_new_copy(r, from_here, to_here)

          call time_message(.false.,time_optimised_loop_1,' Optimised Loop 1')
          do i= 1,3
             call c_redist_32_inv_new_copy(r, from_here, to_here)
          end do
          call time_message(.false.,time_optimised_loop_1,' Optimised Loop 1')

          call c_redist_32_inv_new_opt_copy(r, from_here, to_here)

          call time_message(.false.,time_optimised_loop_2,' Optimised Loop 2')
          do i = 1,3
             call c_redist_32_inv_new_opt_copy(r, from_here, to_here)
          end do
          call time_message(.false.,time_optimised_loop_2,' Optimised Loop 2')

          if(time_optimised_loop_1(1) .gt. time_optimised_loop_2(1)) then
             optimised_choice = 2
          else
             optimised_choice = 1
          end if

!AJ This else is encountered once the optimised auto-tuning choice has 
!AJ been calculated above.
       else
          if(optimised_choice .eq. 1) then
!AJ c_redist_32_inv_new_copy is the new local copy functionality where 
!AJ indirect addressing has largely been removed
             call c_redist_32_inv_new_copy(r, from_here, to_here)       
          else
             call c_redist_32_inv_new_opt_copy(r, from_here, to_here)       
          end if
       end if

!AJ This else is encountered when the optimised local copy functionality is 
!AJ not enabled.
    else
!AJ c_redist_32_inv_old_copy is the original local copy functionality
       call c_redist_32_inv_old_copy(r, from_here, to_here)
    end if
       
!AJ c_redist_32_inv_mpi_copy contains all the remote to local 
!AJ copy functionality
    !<DD>
    if(opt_redist_nbk)then
       if(opt_redist_persist)then
          if(.not.opt_redist_persist_overlap) call c_redist_32_inv_mpi_copy_persist_start(r,from_here)
          call c_redist_32_inv_mpi_copy_persist_end(r,to_here)
       else
          call c_redist_32_inv_mpi_copy_nonblock(r,from_here,to_here)
       endif
    else
       call c_redist_32_inv_mpi_copy(r, from_here, to_here)
    endif

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_32_inv

  subroutine c_redist_32_inv_old_copy(r, from_here, to_here)

    use mp, only: iproc

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i

!CMR 
! In the GS2 standard FFT situation this routine maps 
!         xxf(it,ixxf) to g(ig, isgn, iglo) data type 
!         where it is kx (or x) index, ixxf is (y,ig,isgn,"les") 
!         and iglo is ("xyles") 

    do i = 1, r%to(iproc)%nn
!
! redistribute from local processor to local processor
! NB r%to(iproc)%nn is #elements sent by THIS processor to THIS processor
!    In this situation the data at (r%from(iproc)%k(i),r%from(iproc)%l(i),r%from(iproc)%m(i)) 
!    should come from (r%to(iproc)%k(i),r%to(iproc)%l(i)). 
!
! This do loop, in GS2 standard FFT situation, corresponds to:
!    to_here(ig,isgn,iglo)=from_here(it,ixxf)
!
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

  end subroutine c_redist_32_inv_old_copy

  subroutine c_redist_32_inv_new_copy(r, from_here, to_here)
!=====================================================================
!AJ, June 2011: New code from DCSE project on GS2 Indirect Addressing
!=====================================================================
!
! AJ, June 2011:
!  Modified LOCAL COPY part of c_redist_32_inv routine, as used by GS2 to 
!  transform xxf data type to g data type.
!  Here we REDUCE indirect addressing and cache load by EXPLOITING 
!  understanding of g and xxf data types in GS2.
!  These c_redist_32_inv_* routines use exactly the same code created for 
!  the c_redist_32_* functionality, the only difference is in the actual 
!  data copy the to_here_ and from_here indices are swapped.
!

    use mp, only: iproc

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i,k,t2,t1,f3,f2,f1,fhigh,thigh,f2max
    real :: nakyrecip

    i = 1
!AJ We want to be able to divide by naky with floating point 
!AJ arithmetic so we need to convert naky to a float (hence
!AJ the need for nakyrecip).  Also, to remove the necessity 
!AJ to have a division in the body of the loop below we are 
!AJ calculating the reciprocal of naky so that instead of 
!AJ dividing by naky we can multiply  by 1/naky.
    nakyrecip = naky
    nakyrecip = 1/nakyrecip
    f2max = r%from_high(2)
!AJ Loop over all local copies from THIS proc (iproc) to THIS proc
    do while(i .le. r%to(iproc)%nn)
!AJ Initialise the inner loop indices.  As the to_here array has 
!AJ 3 indices (c_redist_32_inv goes from 2 indices to 3 indices in the
!AJ copy) we need two inner loops to be able to move through all 
!AJ the to_here and from_here indices.
!AJ t1 (which equates to the it index) and f3 (which equates to iglo)
!AJ do not change for the iterations of the inner loops.  f2 (isgn) 
!AJ can only be 1 or 2 so the first inner loop is restricted to at 
!AJ most 2 iterations.
       f2 = r%from(iproc)%l(i)
       f3 = r%from(iproc)%m(i)
       t1 = r%to(iproc)%k(i)
!AJ Ensure that isgn (f2) is either 1 or 2.
       do while (f2 .le. f2max)
!AJ Get initial value of f1 (which equates to ig) and t2 (which equates
!AJ to ixxf).
          f1 = r%from(iproc)%k(i)
          t2 = r%to(iproc)%l(i)
!AJ Work out the maximum value ixxf (t2) can have by calculating the range 
!AJ of ixxf that this process owns.  We step through t2 by naky each iteration
!AJ so the line below calculates thigh as the number of ixxf steps for these 
!AJ inner loops of the operation owned by this process.
          thigh = ceiling(((xxf_lo_ulim_proc+1) - t2)*nakyrecip)
!AJ from_here and to_here are calculated for c_redist_32 using ig, isgn, 
!AJ and iglo.  The calculation of the maximum number if ixxf steps above can, 
!AJ therefore, be used to also calculate the maximum number of ig setsp.  THis 
!AJ is what the following line does.
          thigh = thigh + (f1-1)
!AJ Finally check the actual number of inner loop steps by ensuring that the 
!AJ computed maximum bound of ig is not beyond the actual allowed maximum value 
!AJ in r%from_high(1).
          fhigh = min(thigh,r%from_high(1))
          do k = f1,fhigh
             i = i + 1
             to_here(k,f2,f3) = from_here(t1,t2)
             t2 = t2 + naky
          end do
!AJ If at the end of the inner loop we still have theoretical iterations 
!AJ left on ixxf (i.e. the calculated thigh is higher than the allowed 
!AJ maximum of ig) then move to the next isgn.  If there aren't any 
!AJ iterations left then exit the inner loops (the f2 loop).
          if(thigh .gt. r%from_high(1)) then
             f2 = f2 + 1
          else
             f2 = f2max + 1
          end if
       end do
    end do

  end subroutine c_redist_32_inv_new_copy


  subroutine c_redist_32_inv_new_opt_copy(r, from_here, to_here)
!=====================================================================
!AJ, June 2011: New code from DCSE project on GS2 Indirect Addressing
!=====================================================================
!
! AJ, June 2011:
!  Modified LOCAL COPY part of c_redist_32_inv routine, as used by GS2 to 
!  transform xxf data type to g data type.
!  Here we REDUCE indirect addressing and cache load by EXPLOITING 
!  understanding of g and xxf data types in GS2.
!  This functionality is different to c_redist_32_inv_new_copy because 
!  that routine does not always provide optimsed performance 
!  (due to the way memory reads and writes are performed, particularly 
!  the use of a strided write which triggers a write allocate cache 
!  miss).
!  This new functionality moves through the loop indices in a different 
!  order to the other optimised local copy, so instead of following the 
!  original functionality by moving through the to_here and from_here loops
!  using the i variable from i = 1 to to(iproc)%nn we start at i = 1 and 
!  then skip through i's for a certain range, then move back to i = 2 and 
!  follow the same process until all the i range has been processed.
!  This routine uses the same functionality as teh c_redist_32_new_opt_copy
!  routine, the only difference is in the actual copy the indices of the 
!  from_here and to_here loops are swapped.
!

    use mp, only: iproc
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i,t2,t1,f3,f2,f1
    integer :: f3max,f3maxmultiple,f3incr,innermax,iincrem,iglomax,t1test
    integer :: innermaxmultiplier,outerf3limit,startf3
    real :: innermaxrealvalue,tempnaky

!AJ tempnaky is a real version (rather than integer version) of naky
!AJ It is used to enable real arthimetic rather than integer arthimetic.
    tempnaky = naky
!AJ innermaxmultiplier is the upper limit of the xxf data range owned by this process
    innermaxmultiplier = xxf_lo_ulim_proc+1

!AJ This optimised new local copy functionality is reliant on the particular 
!AJ layout used to work out how to increment through the f3 (iglo) indices.
!AJ f3maxmultiple sets the upper bound of the f3 blocks we step through.
!AJ f3incr sets the amount added to f3 at each step of the inner loops.
    select case (layout)
    case ('yxels')
       f3maxmultiple = naky*ntheta0
       f3incr = naky
    case ('yxles')
       f3maxmultiple = naky*ntheta0
       f3incr = naky
    case ('lexys')
       f3maxmultiple = nlambda*negrid*ntheta0
       f3incr = nlambda*negrid
    case ('lxyes')
       f3maxmultiple = nlambda*ntheta0
       f3incr = nlambda
    case ('lyxes')
       f3maxmultiple = nlambda*naky*ntheta0
       f3incr = nlambda*naky
    case('xyles')
       f3maxmultiple = ntheta0
       f3incr = 1
    end select          
    
    i = 1
!AJ iglomax is the maximum iglo (f3) owned by this process
    iglomax = ((iproc+1)*g_lo_blocksize)
!AJ t1test is used handle the it (t1) index which has a 
!AJ range 1 -> (xxf_lo%ntheta0+1)/2) and then 
!AJ (it - xxf_lo%ntheta0 + xxf_lo%nx) -> xxf_lo%nx.  t1test is used to 
!AJ identify when this change in range happens and enable the code to 
!AJ compensate for it.
    t1test  = ((ntheta0+1)/2)+1
!AJ Loop over all local copies from THIS proc (iproc) to THIS proc
    do while(i .le. r%to(iproc)%nn)

!AJ Initial look up of values needed to calculate innermax
!AJ and setup the first iteration of the loop below (do while i .le. innermax)
       f1 = r%from(iproc)%k(i)
       f2 = r%from(iproc)%l(i)
       f3 = r%from(iproc)%m(i)
       t1 = r%to(iproc)%k(i)
       t2 = r%to(iproc)%l(i)

!AJ outerf3limit calculates the number of f3 iterations that can be undertaken 
!AJ before reaching the next f3 point.
       outerf3limit = f3 + f3incr
       do while(f3 .lt. outerf3limit)
          
!AJ iincrem counts how much i needs to be skipped at the end of the final inner loop
!AJ The i loop below moves through the starting i's for these inner loop iterations
!AJ but the innermost loop increments through further i's without actually moving i
!AJ so once we've finished the loop below we need to move i forward to skip all the 
!AJ elements we have just processed.
          iincrem = i
          
!AJ work out the size of the inner loops body 
!AJ innermax is the range of the t2 index (ixxf_ this process 
!AJ owns. 
          innermaxrealvalue = (innermaxmultiplier-t2)/tempnaky
          innermax = ceiling(innermaxrealvalue)-1
!AJ If the isgn (f2) index is at maximum (2) and the ig (f1) index 
!AJ will also be at maximum if we use the innermax value calculated above 
!AJ (ig maximum is r%from_high(1) which currently is ntgrid in the code) 
!AJ then limit innermax to the difference between ig and ntgrid (which is 
!AJ ig max).
          if(f2 .eq. 2 .and. (f1+innermax) .gt. r%from_high(1)) then
             innermax = i + (ntgrid - f1)
!AJ If the isgn (f2) index is not at maximum, but the calculated innermax would 
!AJ take the indexes beyond whaty is allowed for ig (which is the maximum of ig 
!AJ multipled by each increment of ig) then restrict innermax to the maximum 
!AJ number of allowed ig increments.
          else if((f2 .eq. 1 .and. innermax .gt. (((2*ntgrid)+1)+(ntgrid-f1)))) then
             innermax = i + ((2*ntgrid)+1) + (ntgrid - f1)
!AJ If we have reached this else then the calculated innermax does not breach 
!AJ any of the data limits it can be used for the number of iterations of the 
!AJ inner loop. 
          else
             innermax = i + innermax
          end if
          
          do while(i .le. innermax) 
             
!AJ Get the initial address variables.  There are potentially some extra data lookups 
!AJ here so there may be scope to optimise this, although this has not been checked.
             f1 = r%from(iproc)%k(i)
             f2 = r%from(iproc)%l(i)
             f3 = r%from(iproc)%m(i)
             t1 = r%to(iproc)%k(i)
             t2 = r%to(iproc)%l(i)
             
!AJ Record the initial f3 variable value for the inner loops for use in the bounds checking 
!AJ at the end of the inner loops.             
             startf3 = f3
!AJ Calculate the size of the innermost loop by working out the maximum that f3 (iglo) 
!AJ can go to for this process (using the f3maxmultiple as calculated previously that 
!AJ effectively defines the f3 blocksize of this process).
             f3max = ((f3/f3maxmultiple)+1)*f3maxmultiple
!AJ Ensure that the calculate f3max is not larger than the total iglo space this process
!AJ owns.
             f3max = min(f3max,iglomax)
             do while (f3 .lt. f3max)
                to_here(f1,f2,f3) = from_here(t1,t2)                
                f3 = f3 + f3incr
                t1 = t1 + 1
!AJ Deal with the issues with t1 (it) being a split range variable as described in the 
!AJ comment for t1test.
                if(t1 .eq. t1test) then
                   t1 = t1 - ntheta0 + nx
                end if
                iincrem = iincrem + 1
             end do
             
!AJ Increment the outermost inner loop (move forwards in the i loop one step)
             i = i + 1
                          
          end do

!AJ Ensure that we have not finished the i loop altogether
          if(i .lt. r%from(iproc)%nn .and. iincrem .lt. r%from(iproc)%nn) then
!AJ The f1, f2 and t2 lookups are required to ensure that
!AJ the calculation of innermax the next time round this loop
!AJ are correctly performed as we are now moving forward in the loop
!AJ It may be that it is better to re-order the initialisation of values 
!AJ so this happens at the beginning of the loop rather than at this point.
             f1 = r%from(iproc)%k(i)
             f2 = r%from(iproc)%l(i)
             t2 = r%to(iproc)%l(i)     
!AJ Updated to ensure the f3 loop we are in is correctly setup as with the code 
!AJ above we are moving to a new i so we need to get the correct value of f3 for this i     
             f3 =  startf3 + 1
          else
!AJ If we have gone beyond the end of the i loop then force the exit of the inner loops.
!AJ This should never happen as the loop functionality is designed to always perform the 
!AJ correct number of iterations, so this check could be removed from the code.
             f3 = outerf3limit
          end if

       end do
   
!AJ Once all the inner loops have been iterated move the i variable on to the next point 
!AJ point to be considered.    
       i = iincrem
       
    end do

  end subroutine c_redist_32_inv_new_opt_copy

  subroutine c_redist_32_inv_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here


    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_32_inv_mpi_copy


  !<DD>
  subroutine c_redist_32_inv_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    use mp, only: barrier
    use mp, only: waitall, nbsend, nbrecv
    implicit none
    
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here
    complex, dimension(:,:),allocatable :: r_buff,s_buff
    integer,dimension(:),allocatable :: recvhand, sendhand,recvlookup
    integer :: i, ip, nrecv,count,nsend,j

    !Work out how big recv buffer has to be and allocate it
    nrecv=0 !Could really store this in the redist_object
    do i=0,nproc-1
       if (r%from(i)%nn>0.and.i.ne.iproc) nrecv=nrecv+1
    enddo
    if (nrecv>0) then
       allocate(r_buff(SIZE(r%complex_buff),nrecv))
       allocate(recvhand(nrecv))
       allocate(recvlookup(nrecv))
    endif

    !Now work out how big send buffer has to be, allocate it and fill it
    nsend=0
    do i=0,nproc-1
       if(r%to(i)%nn>0.and.i.ne.iproc) nsend=nsend+1
    enddo
    if (nsend>0) then
       allocate(s_buff(SIZE(r%complex_buff),nsend))
       allocate(sendhand(nsend))
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             do j=1,r%to(i)%nn
                s_buff(j,count)=from_here(r%to(i)%k(j), &
                                              r%to(i)%l(j))
             enddo
          endif
       enddo
    endif

    !Post receives
    if (nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             recvlookup(count)=i
            call nbrecv(r_buff(1:r%from(i)%nn,count),r%from(i)%nn,i,i,recvhand(count))
          endif
       enddo
    endif

    !Now send data where required
    count=0
    do i=0,nproc-1
       if(r%to(i)%nn>0.and.i.ne.iproc) then
          !Note we could think about using RSEND here as we should post all receives first
          !IN practice this will require syncing so probaly no good
          count=count+1
          call nbsend(s_buff(1:r%to(i)%nn,count),r%to(i)%nn,i,iproc,sendhand(count))
       endif
    enddo

    if (nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(nrecv,recvhand)

       !Now unpack data
       do i=1,nrecv
          !Which processor is this data from
          ip=recvlookup(i)
          do j=1,r%from(ip)%nn
             to_here(r%from(ip)%k(j), &
                  r%from(ip)%l(j),r%from(ip)%m(j)) &
                  = r_buff(j,i)
          enddo
       enddo

       deallocate(r_buff)
       deallocate(recvhand)
       deallocate(recvlookup)
    endif

    if (nsend>0) then
       !Shouldn't really need this as we know recvs have worked fine
       call waitall(nsend,sendhand)
       deallocate(s_buff)
       deallocate(sendhand)
    endif

  end subroutine c_redist_32_inv_mpi_copy_nonblock
  !</DD>

  subroutine c_redist_32_inv_mpi_copy_persist_start(r, from_here)
    use mp, only: iproc, nproc
    use mp, only: start_comm
    implicit none
    
    type (redist_type), intent (in out) :: r
    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here
    integer :: i, count,j
    !For inverse transform SEND actually refers to received data and
    !RECV refers to data to be sent

    !Post receives
    if (r%nsend>0) then
       call start_comm(r%recv_inv_hand)
    endif

    !Pack send buffer and start sends
    if (r%nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             do j=1,r%to(i)%nn
                r%cbuff_recv(j,count)=from_here(r%to(i)%k(j), &
                     r%to(i)%l(j))
             enddo
             call start_comm(r%send_inv_hand(count))
          endif
       enddo
    endif
  end subroutine c_redist_32_inv_mpi_copy_persist_start

  subroutine c_redist_32_inv_mpi_copy_persist_end(r, to_here)
    use mp, only: waitall
    implicit none
    
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here
    integer :: i, ip, j
    !For inverse transform SEND actually refers to received data and
    !RECV refers to data to be sent

    !Wait for receives to complete and unpack data
    if (r%nsend>0) then
       call waitall(r%nsend,r%recv_inv_hand)

       !Now unpack data
       do i=1,r%nsend
          !Which processor is this data from
          ip=r%recv_inv_ip(i)
          do j=1,r%from(ip)%nn
             to_here(r%from(ip)%k(j), &
                  r%from(ip)%l(j),r%from(ip)%m(j)) &
                  = r%cbuff_send(j,i)
          enddo
       enddo
    endif

    !Wait for sends to complete
    if (r%nrecv>0) then
       call waitall(r%nrecv,r%send_inv_hand)
    endif
  end subroutine c_redist_32_inv_mpi_copy_persist_end

  subroutine c_redist_42 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_42

  subroutine c_redist_42_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    
    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_42_inv

  subroutine c_redist_23 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

  end subroutine c_redist_23

  subroutine c_redist_34 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    
    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i),&
               r%to(iproc)%n(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i), &
                        r%to(ipfrom)%n(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    
  end subroutine c_redist_34

  subroutine r_redist_12 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):), intent (in) :: from_here

    real, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

  end subroutine r_redist_12

  subroutine r_redist_22 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                     r%from_low(2):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine r_redist_22

  subroutine r_redist_22_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine r_redist_22_inv

  subroutine r_redist_32 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine r_redist_32

  subroutine r_redist_32_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                     r%from_low(2):, &
                     r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine r_redist_32_inv

  subroutine r_redist_42 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%real_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine r_redist_42

  subroutine r_redist_42_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    real, dimension (r%to_low(1):, &
                     r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%real_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute') 
  end subroutine r_redist_42_inv

  subroutine i_redist_12 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):), intent (in) :: from_here

    integer, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_12

  subroutine i_redist_22 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_22

  subroutine i_redist_22_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_22_inv

  subroutine i_redist_32 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_32

  subroutine i_redist_32_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_32_inv

  subroutine i_redist_42 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%integer_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%integer_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_42

  subroutine i_redist_42_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%integer_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%integer_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%integer_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine i_redist_42_inv

  subroutine l_redist_12 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):), intent (in) :: from_here

    logical, dimension (r%to_low(1):,r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i), r%to(iproc)%l(i)) &
            = from_here(r%from(iproc)%k(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i)) 
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_12

  subroutine l_redist_22 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_22

  subroutine l_redist_22_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_22_inv

  subroutine l_redist_32 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_32

  subroutine l_redist_32_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    
    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_32_inv

  subroutine l_redist_42 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i), &
                           r%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%logical_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i), &
                                              r%from(ipto)%n(i))
             end do
             call send (r%logical_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do


    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_42

  subroutine l_redist_42_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    type (redist_type), intent (in out) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    
    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i), &
               r%from(iproc)%n(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%logical_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i), &
                        r%from(ipfrom)%n(i)) &
                        = r%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%logical_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i))
             end do
             call send (r%logical_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine l_redist_42_inv

! TT>
  subroutine c_redist_33 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc,proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here

    integer :: i
    !Timing
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_33_mpi_copy_persist_start(r, from_here)


    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(r%to(iproc)%k(i),&
               r%to(iproc)%l(i),&
               r%to(iproc)%m(i)) &
               = from_here(r%from(iproc)%k(i), &
                           r%from(iproc)%l(i), &
                           r%from(iproc)%m(i))
    end do

    !<DD>
    if(opt_redist_nbk)then
       if(opt_redist_persist)then
          if(.not.opt_redist_persist_overlap) call c_redist_33_mpi_copy_persist_start(r, from_here) !DD
          call c_redist_33_mpi_copy_persist_end(r, to_here) !DD
       else
          call c_redist_33_mpi_copy_nonblock(r, from_here, to_here) !DD
       endif
    else
       call c_redist_33_mpi_copy(r, from_here, to_here)
    endif

    !Timing
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_33

  subroutine c_redist_33_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(r%to(ipfrom)%k(i), &
                        r%to(ipfrom)%l(i), &
                        r%to(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                              r%from(ipto)%l(i), &
                                              r%from(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine c_redist_33_mpi_copy

!<DDHACK>
  subroutine c_redist_33_mpi_copy_nonblock(r, from_here, to_here)
    use mp, only: iproc, nproc, send, receive
    use mp, only: barrier
    use mp, only: waitall, nbsend,nbrecv
    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here
    complex, dimension(:,:),allocatable :: r_buff,s_buff
    integer,dimension(:),allocatable :: recvhand, sendhand,recvlookup
    integer :: i, ip, nrecv,count,nsend,j

    !Work out how big recv buffer has to be and allocate it
    nrecv=0 !Could really store this in the redist_object
    do i=0,nproc-1
       if (r%to(i)%nn>0.and.i.ne.iproc) nrecv=nrecv+1
    enddo
    if (nrecv>0) then
       allocate(r_buff(SIZE(r%complex_buff),nrecv))
       allocate(recvhand(nrecv))
       allocate(recvlookup(nrecv))
    endif

    !Now work out how big send buffer has to be, allocate it and fill it
    nsend=0 !Could store this in the redist_object
    do i=0,nproc-1
       if(r%from(i)%nn>0.and.i.ne.iproc) nsend=nsend+1
    enddo
    if (nsend>0) then
       allocate(s_buff(SIZE(r%complex_buff),nsend))
       allocate(sendhand(nsend))
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             !Any benefit in making this a pointer assignment instead?
             !Avoids copying data here but could lead to questionable
             !memory access during the actual send?
             do j=1,r%from(i)%nn
                s_buff(j,count)=from_here(r%from(i)%k(j), &
                                              r%from(i)%l(j), &
                                              r%from(i)%m(j))
             enddo
             call nbsend(s_buff(1:r%from(i)%nn,count),r%from(i)%nn,i,iproc,sendhand(count))
          endif
       enddo
    endif

!<DD>Moved after sends
    !Post receives
    !Probably not much advantage to doing these ASAP as not going to do anything
    !until data is on its way --> Could post sends first so chance of data being
    !available as soon as receive posted.
    if (nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             recvlookup(count)=i
             call nbrecv(r_buff(1:r%to(i)%nn,count),r%to(i)%nn,i,i,recvhand(count))
          endif
       enddo
    endif

    if (nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(nrecv,recvhand)

       !Now unpack data
       do i=1,nrecv
          !Which processor is this data from? Would we be better looping over all procs and just checking
          !if that's one to unpack? Removes memory lookup of array but requires branching and incrementing
          !a counter to work out the equivalent of i here (i.e. which receive we're looking at)         
          ip=recvlookup(i)

          !Could this be vectorised and is it useful to do so?
          !Could we make r_buff a pointer array so that MPI put data directly into
          !to_here array? Would this be beneficial? Would need to be done before RECVs posted
          do j=1,r%to(ip)%nn
             to_here(r%to(ip)%k(j), &
                  r%to(ip)%l(j),r%to(ip)%m(j)) &
                  = r_buff(j,i)
          enddo
       enddo

       deallocate(r_buff)
       deallocate(recvhand)
       deallocate(recvlookup)
    endif

    if (nsend>0) then
       !Shouldn't really need this as we know recvs have worked fine
       call waitall(nsend,sendhand)
       deallocate(s_buff)
       deallocate(sendhand)
    endif

  end subroutine c_redist_33_mpi_copy_nonblock

  subroutine c_redist_33_mpi_copy_persist_start(r, from_here)
    use mp, only: iproc, nproc
    use mp, only: start_comm
    implicit none
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here
    integer :: i, count,j
    !Start recv comms
    if (r%nrecv>0) then
       call start_comm(r%recv_hand)
    endif

    !Pack buffer and start communication
    if (r%nsend>0) then
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             !Any benefit in making this a pointer assignment instead?
             !Avoids copying data here but could lead to questionable
             !memory access during the actual send?
             do j=1,r%from(i)%nn
                r%cbuff_send(j,count)=from_here(r%from(i)%k(j), &
                                              r%from(i)%l(j), &
                                              r%from(i)%m(j))
             enddo
             call start_comm(r%send_hand(count))
          endif
       enddo
    endif
  end subroutine c_redist_33_mpi_copy_persist_start

  subroutine c_redist_33_mpi_copy_persist_end(r, to_here)
    use mp, only: waitall
    implicit none
    type (redist_type), intent (in out) :: r
    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in out) :: to_here
    integer :: i, ip, j

    !Wait for recv communications to complete and unpack buffer
    if (r%nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(r%nrecv,r%recv_hand)

       !Now unpack data
       do i=1,r%nrecv
          !Which processor is this data from?
          ip=r%recv_ip(i)

          do j=1,r%to(ip)%nn
             to_here(r%to(ip)%k(j), &
                  r%to(ip)%l(j),r%to(ip)%m(j)) &
                  = r%cbuff_recv(j,i)
          enddo
       enddo
    endif

    !Wait for all sends to complete
    if (r%nsend>0) then
       call waitall(r%nsend,r%send_hand)
    endif
  end subroutine c_redist_33_mpi_copy_persist_end
  !</DD>

  subroutine c_redist_33_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_33_inv_mpi_copy_persist_start(r, from_here)

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(r%from(iproc)%k(i), &
               r%from(iproc)%l(i), &
               r%from(iproc)%m(i)) &
               = from_here(r%to(iproc)%k(i), &
                           r%to(iproc)%l(i), &
                           r%to(iproc)%m(i))
    end do

    !<DD>
    if(opt_redist_nbk)then
       if(opt_redist_persist)then
          if(.not.opt_redist_persist_overlap) call c_redist_33_inv_mpi_copy_persist_start(r,from_here)
          call c_redist_33_inv_mpi_copy_persist_end(r,to_here)
       else
          call c_redist_33_inv_mpi_copy_nonblock(r,from_here,to_here)
       endif
    else
       call c_redist_33_inv_mpi_copy(r, from_here, to_here)
    endif

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_33_inv
! <TT

  subroutine c_redist_33_inv_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here


    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(r%from(ipfrom)%k(i), &
                        r%from(ipfrom)%l(i), &
                        r%from(ipfrom)%m(i)) &
                        = r%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                              r%to(ipto)%l(i), &
                                              r%to(ipto)%m(i))
             end do
             call send (r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

  end subroutine c_redist_33_inv_mpi_copy


!<DD>
  subroutine c_redist_33_inv_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    use mp, only: barrier
    use mp, only: waitall, nbsend, nbrecv
    implicit none
    
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in) :: from_here
    complex, dimension(:,:),allocatable :: r_buff,s_buff
    integer,dimension(:),allocatable :: recvhand, sendhand,recvlookup
    integer :: i, ip, nrecv,count,nsend,j

    !Work out how big recv buffer has to be and allocate it
    nrecv=0 !Could really store this in the redist_object
    do i=0,nproc-1
       if (r%from(i)%nn>0.and.i.ne.iproc) nrecv=nrecv+1
    enddo
    if (nrecv>0) then
       allocate(r_buff(SIZE(r%complex_buff),nrecv))
       allocate(recvhand(nrecv))
       allocate(recvlookup(nrecv))
    endif

    !Now work out how big send buffer has to be, allocate it and fill it
    nsend=0
    do i=0,nproc-1
       if(r%to(i)%nn>0.and.i.ne.iproc) nsend=nsend+1
    enddo
    if (nsend>0) then
       allocate(s_buff(SIZE(r%complex_buff),nsend))
       allocate(sendhand(nsend))
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             do j=1,r%to(i)%nn
                s_buff(j,count)=from_here(r%to(i)%k(j), &
                     r%to(i)%l(j), &
                     r%to(i)%m(j))
             enddo
          endif
       enddo
    endif

    !Post receives
    if (nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%from(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             recvlookup(count)=i
            call nbrecv(r_buff(1:r%from(i)%nn,count),r%from(i)%nn,i,i,recvhand(count))
          endif
       enddo
    endif

    !Now send data where required
    count=0
    do i=0,nproc-1
       if(r%to(i)%nn>0.and.i.ne.iproc) then
          !Note we could think about using RSEND here as we should post all receives first
          !IN practice this will require syncing so probaly no good
          count=count+1
          call nbsend(s_buff(1:r%to(i)%nn,count),r%to(i)%nn,i,iproc,sendhand(count))
       endif
    enddo

    if (nrecv>0) then
       !Now we don't have anything to do until all the data has arrived (so should we post sends first?)
       !We could do a wait all but may be better to keep on checking each message seperately so that when
       !it arrives we can copy out whilst waiting for others to arrive. 
       !Start with a wait all
       call waitall(nrecv,recvhand)

       !Now unpack data
       do i=1,nrecv
          !Which processor is this data from
          ip=recvlookup(i)
          do j=1,r%from(ip)%nn
             to_here(r%from(ip)%k(j), &
                  r%from(ip)%l(j),r%from(ip)%m(j)) &
                  = r_buff(j,i)
          enddo
       enddo

       deallocate(r_buff)
       deallocate(recvhand)
       deallocate(recvlookup)
    endif

    if (nsend>0) then
       !Shouldn't really need this as we know recvs have worked fine
       call waitall(nsend,sendhand)
       deallocate(s_buff)
       deallocate(sendhand)
    endif

  end subroutine c_redist_33_inv_mpi_copy_nonblock

  subroutine c_redist_33_inv_mpi_copy_persist_start(r, from_here)
    use mp, only: iproc, nproc
    use mp, only: start_comm
    implicit none
    
    type (redist_type), intent (in out) :: r
    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (in) :: from_here
    integer :: i, count,j
    !For inverse transform SEND actually refers to received data and
    !RECV refers to data to be sent

    !Post receives
    if (r%nsend>0) then
       call start_comm(r%recv_inv_hand)
    endif

    !Pack send buffer and start sends
    if (r%nrecv>0) then
       count=0
       do i=0,nproc-1
          if (r%to(i)%nn>0.and.i.ne.iproc) then
             count=count+1
             do j=1,r%to(i)%nn
                r%cbuff_recv(j,count)=from_here(r%to(i)%k(j), &
                     r%to(i)%l(j), &
                     r%to(i)%m(j))
             enddo
             call start_comm(r%send_inv_hand(count))
          endif
       enddo
    endif
  end subroutine c_redist_33_inv_mpi_copy_persist_start
    

  subroutine c_redist_33_inv_mpi_copy_persist_end(r, to_here)
    use mp, only: waitall
    implicit none
    
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here
    integer :: i, ip, j
    !For inverse transform SEND actually refers to received data and
    !RECV refers to data to be sent

    !Wait for receives to complete and unpack data
    if (r%nsend>0) then
       call waitall(r%nsend,r%recv_inv_hand)

       !Now unpack data
       do i=1,r%nsend
          !Which processor is this data from
          ip=r%recv_inv_ip(i)
          do j=1,r%from(ip)%nn
             to_here(r%from(ip)%k(j), &
                  r%from(ip)%l(j),r%from(ip)%m(j)) &
                  = r%cbuff_send(j,i)
          enddo
       enddo
    endif

    !Wait for sends to complete
    if (r%nrecv>0) then
       call waitall(r%nrecv,r%send_inv_hand)
    endif
  end subroutine c_redist_33_inv_mpi_copy_persist_end
  !</DD>

  subroutine c_redist_36 (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc,proc0,barrier
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):, &
                        r%to_low(6):), intent (in out) :: to_here

    integer :: i
    !Timing
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
!    if(opt_redist_persist_overlap) call c_redist_36_mpi_copy_persist_start(r, from_here)
!AJ4
!    call barrier()



    ! redistribute from local processor to local processor
    do i = 1, r%from(iproc)%nn
       to_here(:,:,r%to(iproc)%m(i), &
                   r%to(iproc)%n(i), &
                   r%to(iproc)%o(i), &
                   r%to(iproc)%p(i)) &
                   = from_here(:,:,r%from(iproc)%m(i))
    end do

    if(opt_redist_nbk)then
!       if(opt_redist_persist)then
!          if(.not.opt_redist_persist_overlap) call c_redist_36_mpi_copy_persist_start(r, from_here) !DD
!          call c_redist_36_mpi_copy_persist_end(r, to_here) !DD
!       else
          call c_redist_36_mpi_copy_nonblock(r, from_here, to_here) !DD
!       endif
    else
       call c_redist_36_mpi_copy(r, from_here, to_here)
    endif

    !Timing
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_36

  subroutine c_redist_36_mpi_copy(r, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    use theta_grid, only : ntgrid
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):, &
                        r%to_low(6):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp
    !AJ Added as the r%complex_buf is not designed for sending contigous blocks, rather to pack individual elements
    !AJ into the buffer.  As each element in my redist object is actually (2*ntgrid+1)*2 (sigma) then I need a bigger
    !AJ send and receive buffer.
    complex, dimension(:,:,:), allocatable :: temp_buff

    allocate(temp_buff(-ntgrid:ntgrid, 1:2, size(r%complex_buff)))

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then
          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn 
                temp_buff(:,:,i) = from_here(:,:,r%from(ipto)%m(i))
             end do
             call send (temp_buff(:,:,1:r%from(ipto)%nn), ipto, idp)
          end if

          !AJ FIX TO DO CONTIGIOUS BLOCKS
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (temp_buff(:,:,1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(:,:,r%to(ipfrom)%m(i), &
                            r%to(ipfrom)%n(i), &
                            r%to(ipfrom)%o(i), &
                            r%to(ipfrom)%p(i)) &
                            = temp_buff(:,:,i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%to(ipfrom)%nn > 0) then
             call receive (temp_buff(:,:,1:r%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%to(ipfrom)%nn
                to_here(:,:,r%to(ipfrom)%m(i), &
                            r%to(ipfrom)%n(i), &
                            r%to(ipfrom)%o(i), &
                            r%to(ipfrom)%p(i)) &
                            = temp_buff(:,:,i)
             end do
          end if
          ! send to idpth next processor
          if (r%from(ipto)%nn > 0) then
             do i = 1, r%from(ipto)%nn
                temp_buff(:,:,i) = from_here(:,:,r%from(ipto)%m(i))
             end do
             call send (temp_buff(:,:,1:r%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

    deallocate(temp_buff)

  end subroutine c_redist_36_mpi_copy

  subroutine c_redist_36_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, nbsend, nbrecv, waitany, waitall, status_size, initialise_requests
    use theta_grid, only : ntgrid
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):, &
                        r%to_low(6):), intent (in out) :: to_here

    integer :: i, ip, nsend, nrecv, count, ipfrom
    integer, dimension(status_size) :: status
    !AJ Added as the r%complex_buf is not designed for sending contigous blocks, rather to pack individual elements
    !AJ into the buffer.  As each element in my redist object is actually (2*ntgrid+1)*2 (sigma) then I need a bigger
    !AJ send and receive buffer.
    complex, dimension(:,:,:,:), allocatable :: recv_buff, send_buff
    integer, dimension(:), allocatable :: recv_hand, send_hand, recv_lookup

    nrecv = 0
    nsend = 0
    do ip = 0, nproc-1
       if(r%to(ip)%nn>0.and.ip.ne.iproc) nrecv=nrecv+1
       if(r%from(ip)%nn>0.and.ip.ne.iproc) nsend=nsend+1
    end do

    if(nrecv>0) then
       allocate(recv_buff(-ntgrid:ntgrid,1:2,size(r%complex_buff),nrecv))
       allocate(recv_hand(nrecv))
       allocate(recv_lookup(nrecv))
       call initialise_requests(recv_hand)
       count = 0
       do ip=0,nproc-1
          if(r%to(ip)%nn>0.and.ip.ne.iproc) then
             count = count + 1
             if(count .gt. nrecv) then
                write(*,*) 'error with count',iproc,count,nrecv
             end if
             recv_lookup(count) = ip
             call nbrecv(recv_buff(-ntgrid:ntgrid,:,1:r%to(ip)%nn,count),ip,ip,recv_hand(count))
          end if
       end do

       if(count .ne. nrecv) then
         write(*,*) 'problem with count and nrecv'
       end if

    end if

    if(nsend>0) then
       allocate(send_buff(-ntgrid:ntgrid,1:2,size(r%complex_buff),nsend))
       allocate(send_hand(nsend))       
       call initialise_requests(send_hand)
       count = 0
       do ip=0,nproc-1
          if(r%from(ip)%nn>0.and.ip.ne.iproc) then
             count = count + 1
             do i = 1, r%from(ip)%nn 
                send_buff(:,:,i,count) = from_here(:,:,r%from(ip)%m(i))
             end do
             call nbsend (send_buff(:,:,1:r%from(ip)%nn,count),ip,iproc,send_hand(count))
          end if
       end do
    end if

    if(nrecv>0) then
       call waitall(nrecv, recv_hand)
       do ip = 1,nrecv
!          call waitany(nrecv, recv_hand, count, status)          
!          ipfrom = recv_lookup(count) 
          ipfrom = recv_lookup(ip) 
          do i = 1, r%to(ipfrom)%nn
             to_here(:,:,r%to(ipfrom)%m(i), &
                  r%to(ipfrom)%n(i), &
                  r%to(ipfrom)%o(i), &
                  r%to(ipfrom)%p(i)) &
!                  = recv_buff(:,:,i,count)
                  = recv_buff(:,:,i,ip)
          end do
       end do
    end if
    
    if(nsend>0) then
       call waitall(nsend, send_hand)
    end if


    if(nsend>0) then
       deallocate(send_buff)
       deallocate(send_hand)
    end if

    if(nrecv>0) then
       deallocate(recv_buff)
       deallocate(recv_hand)
       deallocate(recv_lookup)
    end if


  end subroutine c_redist_36_mpi_copy_nonblock


  subroutine c_redist_36_inv (r, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, proc0
    type (redist_type), intent (in out) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):, &
                        r%to_low(6):), intent (in) :: from_here

    integer :: i

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
!    if(opt_redist_persist_overlap) call c_redist_36_inv_mpi_copy_persist_start(r, from_here)

    ! redistribute from local processor to local processor
    do i = 1, r%to(iproc)%nn
       to_here(:,:,r%from(iproc)%m(i)) &
               = from_here(:,:,r%to(iproc)%m(i),&
                               r%to(iproc)%n(i),&
                               r%to(iproc)%o(i),&
                               r%to(iproc)%p(i))
    end do

    !<DD>
    if(opt_redist_nbk)then
!       if(opt_redist_persist)then
!          if(.not.opt_redist_persist_overlap) call c_redist_33_inv_mpi_copy_persist_start(r,from_here)
!          call c_redist_36_inv_mpi_copy_persist_end(r,to_here)
!       else
          call c_redist_36_inv_mpi_copy_nonblock(r,from_here,to_here)
!       endif
    else
       call c_redist_36_inv_mpi_copy(r, from_here, to_here)
    endif

    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_redist_36_inv

  subroutine c_redist_36_inv_mpi_copy(r, from_here, to_here)

    use theta_grid, only : ntgrid
    use mp, only: iproc, nproc, send, receive

    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):, &
                        r%to_low(6):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here


    integer :: i, idp, ipto, ipfrom, iadp
    !AJ Added as the r%complex_buf is not designed for sending contigous blocks, rather to pack individual elements
    !AJ into the buffer.  As each element in my redist object is actually (2*ntgrid+1)*2 (sigma) then I need a bigger
    !AJ send and receive buffer.
    complex, dimension(:,:,:), allocatable :: temp_buff

    allocate(temp_buff(-ntgrid:ntgrid,2,size(r%complex_buff))) 

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                temp_buff(:,:,i) = from_here(:,:,r%to(ipto)%m(i), &
                                                 r%to(ipto)%n(i), &
                                                 r%to(ipto)%o(i), &
                                                 r%to(ipto)%p(i))
             end do
             call send (temp_buff(:,:,1:r%to(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (temp_buff(:,:,1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(:,:,r%from(ipfrom)%m(i)) = temp_buff(:,:,i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (r%from(ipfrom)%nn > 0) then
             call receive (temp_buff(:,:,1:r%from(ipfrom)%nn), ipfrom, idp)
             do i = 1, r%from(ipfrom)%nn
                to_here(:,:,r%from(ipfrom)%m(i)) = temp_buff(:,:,i)
             end do
          end if

          ! send to idpth next processor
          if (r%to(ipto)%nn > 0) then
             do i = 1, r%to(ipto)%nn
                temp_buff(:,:,i) = from_here(:,:,r%to(ipto)%m(i), &
                                                 r%to(ipto)%n(i), &
                                                 r%to(ipto)%o(i), &
                                                 r%to(ipto)%p(i))
             end do
             call send (temp_buff(:,:,1:r%to(ipto)%nn), ipto, idp)
          end if

       end if
    end do

    deallocate(temp_buff)

  end subroutine c_redist_36_inv_mpi_copy

  subroutine c_redist_36_inv_mpi_copy_nonblock(r, from_here, to_here)

    use mp, only: iproc, nproc, nbsend, nbrecv, waitany, waitall, status_size
    use theta_grid, only : ntgrid
    type (redist_type), intent (in out) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):, &
                        r%to_low(6):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in out) :: to_here

    integer :: i, ip, nsend, nrecv, count, ipfrom
    integer, dimension(status_size) :: status
    !AJ Added as the r%complex_buf is not designed for sending contigous blocks, rather to pack individual elements
    !AJ into the buffer.  As each element in my redist object is actually (2*ntgrid+1)*2 (sigma) then I need a bigger
    !AJ send and receive buffer.
    complex, dimension(:,:,:,:), allocatable :: recv_buff, send_buff
    integer, dimension(:), allocatable :: recv_hand, send_hand, recv_lookup

    nrecv = 0
    nsend = 0
    do ip = 0, nproc-1
       if(r%from(ip)%nn>0.and.ip.ne.iproc) nrecv=nrecv+1
       if(r%to(ip)%nn>0.and.ip.ne.iproc) nsend=nsend+1
    end do

    if(nrecv>0) then
       allocate(recv_buff(-ntgrid:ntgrid,1:2,size(r%complex_buff),nrecv))
       allocate(recv_hand(nrecv))
       allocate(recv_lookup(nrecv))

       count = 0
       do ip=0,nproc-1
          if(r%from(ip)%nn>0.and.ip.ne.iproc) then
             count = count + 1
             recv_lookup(count) = ip
             call nbrecv(recv_buff(-ntgrid:ntgrid,:,1:r%from(ip)%nn,count),ip,ip,recv_hand(count))
          end if
       end do
    end if

    if(nsend>0) then
       allocate(send_buff(-ntgrid:ntgrid,1:2,size(r%complex_buff),nsend))
       allocate(send_hand(nsend))       

       count = 0
       do ip=0,nproc-1
          if(r%to(ip)%nn>0.and.ip.ne.iproc) then
             count = count + 1
             do i = 1, r%to(ip)%nn 
                send_buff(:,:,i,count) = from_here(:,:,r%to(ip)%m(i), &
                                                       r%to(ip)%n(i), &
                                                       r%to(ip)%o(i), &
                                                       r%to(ip)%p(i))
             end do
             call nbsend (send_buff(:,:,1:r%to(ip)%nn,count),ip,iproc,send_hand(count))
          end if
       end do
    end if

    if(nrecv>0) then
       do ip = 1,nrecv
          call waitany(nrecv, recv_hand, count, status)          
          ipfrom = recv_lookup(count) 
          do i = 1, r%from(ipfrom)%nn
             to_here(:,:,r%from(ipfrom)%m(i)) = recv_buff(:,:,i,count)
          end do
       end do
    end if
    
    call waitall(nsend, send_hand)

    if(nsend>0) then
       deallocate(send_buff)
       deallocate(send_hand)
    end if

    if(nrecv>0) then
       deallocate(recv_buff)
       deallocate(recv_hand)
       deallocate(recv_lookup)
    end if


  end subroutine c_redist_36_inv_mpi_copy_nonblock


  subroutine c_fill_2 (f, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    implicit none
    type (redist_type), intent (in out) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_fill_2

  subroutine c_fill_3 (f, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, proc0
    implicit none
    type (redist_type), intent (in out) :: f
    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here
    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here
    integer :: i
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')

    !If overlapping then start comms now 
    if(opt_redist_persist_overlap) call c_redist_33_mpi_copy_persist_start(f, from_here)

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do
    
    !Do mpi part of copy
    if(opt_redist_nbk)then
       if(opt_redist_persist)then
          if(.not.opt_redist_persist_overlap) call c_redist_33_mpi_copy_persist_start(f,from_here)
          call c_redist_33_mpi_copy_persist_end(f,to_here)
       else
          call c_fill_3_mpi_copy_nonblock(f,from_here,to_here)
       endif
    else
       call c_fill_3_mpi_copy(f,from_here,to_here)
    endif
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_fill_3

  subroutine c_fill_3_mpi_copy (f, from_here, to_here)
    use mp, only: iproc, nproc, send, receive
    implicit none
    type (redist_type), intent (in out) :: f
    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here
    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here
    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do
  end subroutine c_fill_3_mpi_copy

  !>Non-blocking version of mpi copy
  subroutine c_fill_3_mpi_copy_nonblock (f, from_here, to_here)
    type (redist_type), intent (in out) :: f
    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here
    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here
    call c_redist_33_mpi_copy_nonblock(f,from_here,to_here)
  end subroutine c_fill_3_mpi_copy_nonblock

  subroutine c_fill_4 (f, from_here, to_here)
    use job_manage, only: time_message
    use mp, only: iproc, nproc, send, receive, proc0
    implicit none
    type (redist_type), intent (in out) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%complex_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%complex_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%complex_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%complex_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%complex_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do
    if (proc0.and.(.not.using_measure_scatter)) call time_message(.false.,time_redist,' Redistribute')
  end subroutine c_fill_4

  subroutine r_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_fill_2

  subroutine r_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_fill_3

  subroutine r_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%real_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%real_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%real_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%real_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%real_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine r_fill_4

  subroutine i_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_fill_2

  subroutine i_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_fill_3

  subroutine i_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%integer_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%integer_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%integer_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%integer_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%integer_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine i_fill_4

  subroutine l_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i)) &
                        = f%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_fill_2

  subroutine l_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i)) &
                        = f%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_fill_3

  subroutine l_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, send, receive
    type (redist_type), intent (in out) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (in out) :: to_here

    integer :: i, idp, ipto, ipfrom, iadp

    ! redistribute from local processor to local processor
    do i = 1, f%from(iproc)%nn
       to_here(f%to(iproc)%k(i),&
               f%to(iproc)%l(i), &
               f%to(iproc)%m(i), &
               f%to(iproc)%n(i)) &
               = from_here(f%from(iproc)%k(i), &
                           f%from(iproc)%l(i), &
                           f%from(iproc)%m(i), &
                           f%from(iproc)%n(i))
    end do

    ! redistribute to idpth next processor from idpth preceding processor
    ! or redistribute from idpth preceding processor to idpth next processor
    ! to avoid deadlocks
    do idp = 1, nproc-1
       ipto = mod(iproc + idp, nproc)
       ipfrom = mod(iproc + nproc - idp, nproc)
       iadp = min(idp, nproc - idp)
       ! avoid deadlock AND ensure mostly parallel resolution
       if (mod(iproc/iadp,2) == 0) then

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if

          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%logical_buff(i)
             end do
          end if
       else
          ! receive from idpth preceding processor
          if (f%to(ipfrom)%nn > 0) then
             call receive (f%logical_buff(1:f%to(ipfrom)%nn), ipfrom, idp)
             do i = 1, f%to(ipfrom)%nn
                to_here(f%to(ipfrom)%k(i), &
                        f%to(ipfrom)%l(i), &
                        f%to(ipfrom)%m(i), &
                        f%to(ipfrom)%n(i)) &
                        = f%logical_buff(i)
             end do
          end if

          ! send to idpth next processor
          if (f%from(ipto)%nn > 0) then
             do i = 1, f%from(ipto)%nn
                f%logical_buff(i) = from_here(f%from(ipto)%k(i), &
                                              f%from(ipto)%l(i), &
                                              f%from(ipto)%m(i), &
                                              f%from(ipto)%n(i))
             end do
             call send (f%logical_buff(1:f%from(ipto)%nn), ipto, idp)
          end if
       end if
    end do

  end subroutine l_fill_4

  subroutine report_map_property (r)

    use mp, only: iproc, nproc, proc0, sum_reduce, max_reduce
    type (redist_type), intent (in) :: r
    type :: redist_prp
       integer :: local_max, local_total
       integer :: comm_max, comm_total
       integer :: elm_max, elm_total
    end type redist_prp
    type (redist_prp) :: prp
    integer :: ip, rank_from, rank_to
    integer, dimension(:), allocatable :: lbd_from, lbd_to

    prp%comm_max = 0
    prp%comm_total = 0
    prp%elm_total = 0

    do ip=0, nproc-1
       if (ip == iproc) then
          prp%local_total = r%to(ip)%nn
       else
          if (r%to(ip)%nn > 0) then
             prp%comm_total = prp%comm_total + 1
             prp%elm_total = prp%elm_total + r%to(ip)%nn
          end if
       end if
    end do
    prp%local_max = prp%local_total
    prp%comm_max = prp%comm_total
    prp%elm_max = prp%elm_total

    call max_reduce (prp%local_max, 0)
    call sum_reduce (prp%local_total, 0)
    call max_reduce (prp%comm_max, 0)
    call sum_reduce (prp%comm_total, 0)
    call max_reduce (prp%elm_max, 0)
    call sum_reduce (prp%elm_total, 0)

    if (proc0) then
       rank_from = 1
       if (associated(r%from(0)%l)) rank_from = 2
       if (associated(r%from(0)%m)) rank_from = 3
       if (associated(r%from(0)%n)) rank_from = 4
       rank_to = 1
       if (associated(r%to(0)%l)) rank_to = 2
       if (associated(r%to(0)%m)) rank_to = 3
       if (associated(r%to(0)%n)) rank_to = 4
       allocate (lbd_from(rank_from), lbd_to(rank_to))
       lbd_from = r%from_low(1:rank_from)
       lbd_to = r%to_low(1:rank_to)
       print '(a,i2,a,i2)', 'From rank', rank_from, ' to rank', rank_to
       print '(a,t20,4i10)', 'From lbound (proc0)', r%from_low(1:rank_from)
       print '(a,t20,4i10)', 'To lbound (proc0)',   r%to_low(1:rank_to)
       print '(a,t49,a,t64,a)', '--- Redistribution statistics ---', &
            'max', 'avg'
       print '(a,t40,i12,t55,f15.2)', 'Number of local move elements', &
            prp%local_max, real(prp%local_total)/real(nproc)
       print '(a,t40,i12,t60,f10.2)', &
            'Number of inter-processor communications', &
            prp%comm_max, real(prp%comm_total)/real(nproc)
       print '(a,t40,i12,t55,f15.2)', &
            'Number of inter-processor move elements', &
            prp%elm_max, real(prp%elm_total)/real(nproc)
       print *
    end if

  end subroutine report_map_property

  subroutine measure_gather_32 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_32 (map, gin, gout)
    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    gather_count = gather_count + 1

  end subroutine measure_gather_32

  subroutine measure_scatter_23 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:,:), intent (out) :: gout

    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_32_inv (map, gin, gout)
    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    scatter_count = scatter_count + 1

  end subroutine measure_scatter_23

  subroutine measure_gather_33 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:,:), intent (in) :: gin
    complex, dimension (:,:,:), intent (out) :: gout

    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_33 (map, gin, gout)
    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    gather_count = gather_count + 1

  end subroutine measure_gather_33

  subroutine measure_scatter_33 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:,:), intent (in) :: gin
    complex, dimension (:,:,:), intent (out) :: gout

    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_33_inv (map, gin, gout)
    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    scatter_count = scatter_count + 1

  end subroutine measure_scatter_33

  subroutine measure_gather_22 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_22 (map, gin, gout)
    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    gather_count = gather_count + 1

  end subroutine measure_gather_22

  subroutine measure_scatter_22 (map, gin, gout)

    use job_manage, only: time_message
    use mp, only: proc0

    type (redist_type), intent (in out) :: map
    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    call c_redist_22_inv (map, gin, gout)
    if (proc0.and.using_measure_scatter) call time_message(.false.,time_redist,' Redistribution')
    scatter_count = scatter_count + 1

  end subroutine measure_scatter_22

end module redistribute
