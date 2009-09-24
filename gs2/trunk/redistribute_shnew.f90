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
  private

  public :: index_list_type, delete_list
  public :: redist_type, delete_redist

  public :: init_redist, gather, scatter
  public :: init_fill, fill

  interface gather
     module procedure c_redist_22, r_redist_22, i_redist_22, l_redist_22
     module procedure c_redist_32, r_redist_32, i_redist_32, l_redist_32
     module procedure c_redist_42, r_redist_42, i_redist_42, l_redist_42
     module procedure c_redist_23, r_redist_23, i_redist_23, l_redist_23
     module procedure c_redist_34, r_redist_34, i_redist_34, l_redist_34
  end interface

  interface scatter
     module procedure c_redist_12,     r_redist_12,     i_redist_12,     l_redist_12
     module procedure c_redist_22_inv, r_redist_22_inv, i_redist_22_inv, l_redist_22_inv
     module procedure c_redist_32_inv, r_redist_32_inv, i_redist_32_inv, l_redist_32_inv
     module procedure c_redist_42_inv, r_redist_42_inv, i_redist_42_inv, l_redist_42_inv
  end interface

  interface fill
     module procedure c_fill_2, c_fill_3, c_fill_4
     module procedure r_fill_2, r_fill_3, r_fill_4
     module procedure i_fill_2, i_fill_3, i_fill_4
     module procedure l_fill_2, l_fill_3, l_fill_4
  end interface
!
! For one-sided communication, the j index is the inverse of the i index
!
  type :: index_map
     integer :: ni, nj
     integer, dimension (:), pointer :: i => null ()
     integer, dimension (:), pointer :: j => null ()
  end type index_map

  type :: redist_type
     private
     integer, dimension(4) :: to_low, from_low
     type (index_map), dimension (:), pointer :: to => null ()
     type (index_map), dimension (:), pointer :: from => null ()
     integer :: redist_buff_size
  end type redist_type
  
  type :: index_list_type
     integer, dimension(:), pointer :: first => null ()
     integer, dimension(:), pointer :: second => null ()
     integer, dimension(:), pointer :: third => null ()
     integer, dimension(:), pointer :: fourth => null ()
  end type index_list_type

contains

  subroutine init_redist(r, char, to_low, to_high, to_list, &
       from_low, from_high, from_list, ierr_tmp)

    use mp, only: iproc, nproc, proc0
    implicit none
    type (redist_type), intent (out) :: r
    character(1), intent (in) :: char
    integer, intent (in) :: to_low
    integer, dimension (:), intent (in) :: to_high, from_high
    type (index_list_type), dimension(0:) :: to_list, from_list
    integer, dimension(:), intent (in) :: from_low

    integer, parameter :: one   = 1
    integer, parameter :: two   = 2
    integer, parameter :: three = 3
    integer, parameter :: four  = 4
    
    integer, intent (out), optional :: ierr_tmp
    integer :: j, ip, buff_size, ierr, ni, nj
    integer :: idp, ipfrom, ipto
    integer :: d1, d2, d3

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    allocate (r%to(0:nproc-1), r%from(0:nproc-1))

    ierr = 0
    buff_size = 0

    r%to_low(1) = 1
    r%to_low(2) = to_low

    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          nj = size(to_list(ip)%first)

          r%to(ip)%nj = nj
          r%from(ip)%nj = nj

          allocate (r%to(ip)%j(nj))
          allocate (r%from(ip)%j(nj))
             
          buff_size = max(buff_size, nj)
       else
          r%to(ip)%nj = 0
          r%from(ip)%nj = 0
       endif
    enddo

    do j = 1, size(from_low)
       r%from_low(j) = from_low(j)
    enddo

    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then
          ni = size(from_list(ip)%first)

          r%from(ip)%ni = ni
          r%to(ip)%ni = ni

          allocate (r%from(ip)%i(ni))
          allocate (r%to(ip)%i(ni))
                       
          buff_size = max(buff_size, ni)
       else
          r%from(ip)%ni = 0
          r%to(ip)%ni = 0
       endif
    enddo

    r%redist_buff_size = buff_size

    select case (size(from_high))
       case (one) 
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                r%from(ip)%i = from_list(ip)%first - from_low(1)
             end if
          end do
       case (two) 
          d1 = from_high(1)-from_low(1)+1
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                r%from(ip)%i = from_list(ip)%first  - from_low(1) &
                     + (from_list(ip)%second - from_low(2))*d1
             end if
          end do
       case (three) 
          d1 = from_high(1)-from_low(1)+1
          d2 = d1*(from_high(2)-from_low(2)+1)
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                r%from(ip)%i = from_list(ip)%first  - from_low(1) &
                     + (from_list(ip)%second - from_low(2))*d1 &
                     + (from_list(ip)%third  - from_low(3))*d2
             end if
          end do
       case (four) 
          d1 = from_high(1)-from_low(1)+1
          d2 = d1*(from_high(2)-from_low(2)+1)
          d3 = d2*(from_high(3)-from_low(3)+1)
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                r%from(ip)%i = from_list(ip)%first  - from_low(1) &
                     + (from_list(ip)%second - from_low(2))*d1 &
                     + (from_list(ip)%third  - from_low(3))*d2 &
                     + (from_list(ip)%fourth - from_low(4))*d3
             end if
          end do
       case default
          ierr = 1
    end select

    select case (size(to_high))
       case (one) 
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                r%to(ip)%j = to_list(ip)%first - r%to_low(1)
             end if
          end do
       case (two) 
          d1 = to_high(1)-r%to_low(1)+1
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                r%to(ip)%j = to_list(ip)%first  - r%to_low(1) &
                     + (to_list(ip)%second - r%to_low(2))*d1
             end if
          end do
       case (three) 
! this option is not fully wired up
          ierr = 1
          d1 = to_high(1)-r%to_low(1)+1
          d2 = d1*(to_high(2)-r%to_low(2)+1)
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                r%to(ip)%j = to_list(ip)%first  - r%to_low(1) &
                     + (to_list(ip)%second - r%to_low(2))*d1 &
                     + (to_list(ip)%third  - r%to_low(3))*d2
             end if
          end do
       case (four) 
! this option is not fully wired up
          ierr = 1
          d1 = to_high(1)-r%to_low(1)+1
          d2 = d1*(to_high(2)-r%to_low(2)+1)
          d3 = d2*(to_high(3)-r%to_low(3)+1)
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                r%to(ip)%j = to_list(ip)%first  - r%to_low(1) &
                     + (to_list(ip)%second - r%to_low(2))*d1 &
                     + (to_list(ip)%third  - r%to_low(3))*d2 &
                     + (to_list(ip)%fourth - r%to_low(4))*d3
             end if
          end do
       case default
          ierr = 1
    end select

    if (char == 'c') call reindex (r)

! Now fill in inverse structures

! Loop over PEs
    do idp = 0, nproc-1
       ipto = mod (idp + iproc, nproc)
       ipfrom = mod (iproc + nproc - idp, nproc)

! load local address of data to be fetched
       if (r%to(ipfrom)%nj > 0) iy = loc (r%to(ipfrom)%j)
       call barrier

       if (r%to(ipto)%ni > 0) then

! get address of data on remote PE
          call shmem_get (ix, iy, 1, ipto)
! get data from remote PE
          call shmem_get (r%to(ipto)%i, x, r%to(ipto)%ni, ipto)

       end if
       call barrier
    end do

! Loop over PEs
    do idp = 0, nproc-1
       ipto = mod (idp + iproc, nproc)
       ipfrom = mod (iproc + nproc - idp, nproc)

! load local address of data to be fetched
       if (r%from(ipfrom)%ni > 0) iy = loc (r%from(ipfrom)%i)
       call barrier

       if (r%from(ipto)%nj > 0) then

! get address of data on remote PE
          call shmem_get (ix, iy, 1, ipto)
! get data from remote PE
          call shmem_get (r%from(ipto)%j, x, r%from(ipto)%nj, ipto)

       end if
       call barrier
    end do

    if (present(ierr_tmp)) ierr_tmp = ierr

  end subroutine init_redist

  subroutine init_fill (f, char, to_low, to_high, to_list, &
       from_low, from_high, from_list, ierr_tmp)

    use mp, only: nproc, proc0, iproc
    implicit none
    type (redist_type), intent (out) :: f
    character(1), intent (in) :: char
    type (index_list_type), dimension(0:) :: to_list, from_list
    integer, dimension(:), intent (in) :: to_low,  from_low, &
                                          to_high, from_high

    type (index_map), dimension(0:nproc-1) :: to_tmp

    integer, parameter :: one   = 1
    integer, parameter :: two   = 2
    integer, parameter :: three = 3
    integer, parameter :: four  = 4
    
    integer, intent (out), optional :: ierr_tmp
    integer :: j, ip, n_to, n_from, buff_size, ierr, ni, nj
    integer :: idp, ipfrom, ipto
    integer :: d1, d2, d3

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    do j = 1, size(to_low)
       f%to_low(j) = to_low(j)
    enddo

    do j = 1, size(from_low)
       f%from_low(j) = from_low(j)
    enddo

    allocate (f%to(0:nproc-1), f%from(0:nproc-1))

    ierr = 0
    buff_size = 0

    do ip = 0, nproc - 1
       if (associated(to_list(ip)%first)) then
          nj = size(to_list(ip)%first)

          f%to(ip)%nj = nj
          f%from(ip)%nj = nj

          allocate (f%to(ip)%j(nj))
          allocate (f%from(ip)%j(nj))

          buff_size = max(buff_size, nj)
       else
          f%to(ip)%nj = 0
          f%from(ip)%nj = 0
       endif
    enddo

    do ip = 0, nproc - 1
       if (associated(from_list(ip)%first)) then
          ni = size(from_list(ip)%first)

          f%from(ip)%ni = ni
          f%to(ip)%ni = ni

          allocate (f%from(ip)%i(ni))
          allocate (f%to(ip)%i(ni))
          
          buff_size = max(buff_size, ni)
       else
          f%from(ip)%ni = 0
          f%to(ip)%ni = 0
       endif
    enddo

    f%redist_buff_size = buff_size
   
    select case (size(from_high))
       case (one) 
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                f%from(ip)%i = from_list(ip)%first - from_low(1)
             end if
          end do
       case (two) 
          d1 = from_high(1)-from_low(1)+1
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                f%from(ip)%i = from_list(ip)%first  - from_low(1) &
                     + (from_list(ip)%second - from_low(2))*d1
             end if
          end do
       case (three) 
          d1 = from_high(1)-from_low(1)+1
          d2 = d1*(from_high(2)-from_low(2)+1)
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then
                f%from(ip)%i = from_list(ip)%first  - from_low(1) &
                     + (from_list(ip)%second - from_low(2))*d1 &
                     + (from_list(ip)%third  - from_low(3))*d2
             end if
          end do
       case (four) 
          d1 = from_high(1)-from_low(1)+1
          d2 = d1*(from_high(2)-from_low(2)+1)
          d3 = d2*(from_high(3)-from_low(3)+1)
          do ip = 0, nproc-1 
             if (associated (from_list(ip)%first)) then             
                f%from(ip)%i = from_list(ip)%first  - from_low(1) &
                     + (from_list(ip)%second - from_low(2))*d1 &
                     + (from_list(ip)%third  - from_low(3))*d2 &
                     + (from_list(ip)%fourth - from_low(4))*d3
             end if
          end do
       case default
          ierr = 1
    end select

    select case (size(to_high))
       case (one) 
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                f%to(ip)%j = to_list(ip)%first - to_low(1)
             end if
          end do
       case (two) 
          d1 = to_high(1)-to_low(1)+1
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                f%to(ip)%j = to_list(ip)%first  - to_low(1) &
                     + (to_list(ip)%second - to_low(2))*d1
             end if
             end do
       case (three) 
          d1 = to_high(1)-to_low(1)+1
          d2 = d1*(to_high(2)-to_low(2)+1)
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                f%to(ip)%j = to_list(ip)%first  - to_low(1) &
                     + (to_list(ip)%second - to_low(2))*d1 &
                     + (to_list(ip)%third  - to_low(3))*d2
             end if
          end do
       case (four) 
          d1 = to_high(1)-to_low(1)+1
          d2 = d1*(to_high(2)-to_low(2)+1)
          d3 = d2*(to_high(3)-to_low(3)+1)
          do ip = 0, nproc-1 
             if (associated (to_list(ip)%first)) then
                f%to(ip)%j = to_list(ip)%first  - to_low(1) &
                     + (to_list(ip)%second - to_low(2))*d1 &
                     + (to_list(ip)%third  - to_low(3))*d2 &
                     + (to_list(ip)%fourth - to_low(4))*d3
             end if
          end do
       case default
          ierr = 1
    end select

    if (char == 'c') call reindex (f)

! Now fill in f%to structure

! Loop over PEs
    do idp = 0, nproc-1
       ipto = mod (idp + iproc, nproc)
       ipfrom = mod (iproc + nproc - idp, nproc)

! load local address of data to be fetched
       if (f%to(ipfrom)%nj > 0) iy = loc (f%to(ipfrom)%j)
       call barrier

       if (f%to(ipto)%ni > 0) then

! get address of data on remote PE
          call shmem_get (ix, iy, 1, ipto)
! get data from remote PE
          call shmem_get (f%to(ipto)%i, x, f%to(ipto)%ni, ipto)

       endif
       call barrier
    end do


! Loop over PEs
    do idp = 0, nproc-1
       ipto = mod (idp + iproc, nproc)
       ipfrom = mod (iproc + nproc - idp, nproc)

! load local address of data to be fetched
       if (f%from(ipfrom)%ni > 0) iy = loc (f%from(ipfrom)%i)
       call barrier

       if (f%from(ipto)%nj > 0) then

! get address of data on remote PE
          call shmem_get (ix, iy, 1, ipto)
! get data from remote PE
          call shmem_get (f%from(ipto)%j, x, f%from(ipto)%nj, ipto)

       end if 
       call barrier
    end do

    if (present(ierr_tmp)) ierr_tmp = ierr

  end subroutine init_fill

  subroutine delete_redist(r)

    use mp, only: nproc
    type (redist_type), intent (in out) :: r
    
    integer :: n

    if (associated(r%to)) then
       do n = 0, nproc-1
          if (associated (r%to(n)%i)) deallocate (r%to(n)%i)
          if (associated (r%to(n)%j)) deallocate (r%to(n)%j)
       enddo
       deallocate (r%to)
    endif

    if (associated(r%from)) then
       do n = 0, nproc-1
          if (associated (r%from(n)%i)) deallocate (r%from(n)%i)
          if (associated (r%from(n)%j)) deallocate (r%from(n)%j)
       enddo
       deallocate (r%from)
    endif

  end subroutine delete_redist

  subroutine delete_list(list)
    use mp, only: nproc
    type (index_list_type), dimension(0:) :: list

    integer :: ip

    do ip = 0, nproc-1
       if(associated(list(ip)%first))  deallocate(list(ip)%first)
       if(associated(list(ip)%second)) deallocate(list(ip)%second)
       if(associated(list(ip)%third))  deallocate(list(ip)%third)
       if(associated(list(ip)%fourth)) deallocate(list(ip)%fourth)
    enddo

  end subroutine delete_list

  subroutine reindex (r)
!
! SHMEM limitation for complex data: work around with re-indexing of data
! Use to_tmp k index to store r%from%i; re-index data as if real
!
    use mp, only: nproc
    type (redist_type) :: r
    integer, dimension(:), allocatable :: tmp
    integer :: ip, ni, nj, n

    do ip = 0, nproc-1
! Do from structure first
       ni = r%from(ip)%ni
       if (ni > 0) then
! make tmp space
          allocate (tmp(ni))
! copy data to tmp space
          tmp = r%from(ip)%i
! reallocate main index array
          deallocate (r%from(ip)%i)
          allocate (r%from(ip)%i(2*ni))
! change indices
          do n = 1, ni
             r%from(ip)%i(2*n-1) = 2 * tmp(n)
             r%from(ip)%i(2*n  ) = 2 * tmp(n) + 1
          end do
! delete tmp space
          deallocate (tmp)
! double number of indices to be shipped around later
          r%from(ip)%ni = 2 * ni 
! take care of inverse indices
          deallocate (r%to(ip)%i)
          r%to(ip)%ni = 2 * ni
          allocate (r%to(ip)%i(2*ni))
       end if

! Now do to structure
       nj = r%to(ip)%nj
       if (nj > 0) then
! make tmp space
          allocate (tmp(nj))
! copy data to tmp space
          tmp = r%to(ip)%j
! reallocate main index array
          deallocate (r%to(ip)%j)
          allocate (r%to(ip)%j(2*nj))
! change indices
          do n = 1, nj
             r%to(ip)%j(2*n-1) = 2 * tmp(n)
             r%to(ip)%j(2*n  ) = 2 * tmp(n) + 1
          end do
! delete tmp space
          deallocate (tmp)
! double number of indices to be shipped around later
          r%to(ip)%nj = 2 * nj
! take care of inverse indices
          deallocate (r%from(ip)%j)
          r%from(ip)%nj = 2 * nj
          allocate (r%from(ip)%j(2*nj))
       endif
    end do
      
  end subroutine reindex

  subroutine c_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%from_low(1):), intent (in) :: from_here

    complex, dimension (r%to_low(1):,r%to_low(2):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_12

  subroutine c_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_22

  subroutine c_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, complex_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_22_inv

  subroutine c_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_32

  subroutine c_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, complex_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_32_inv

  subroutine c_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_42

  subroutine c_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, complex_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_42_inv

  subroutine c_redist_23 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_23

  subroutine c_redist_34 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    complex, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    complex, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (out) :: to_here

    complex, dimension (r%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_redist_34

  subroutine c_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (out) :: to_here

    complex, dimension (f%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_fill_2

  subroutine c_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (out) :: to_here

    complex, dimension (f%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_fill_3

  subroutine c_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    complex, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    complex, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (out) :: to_here

    complex, dimension (f%redist_buff_size) :: complex_buff
    integer :: ip, idp

    pointer (ix, x)
    complex x
    
    pointer (iy, y)
    complex y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (complex_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, complex_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine c_fill_4

  subroutine r_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%from_low(1):), intent (in) :: from_here

    real, dimension (r%to_low(1):,r%to_low(2):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_12

  subroutine r_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_22

  subroutine r_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, real_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_22_inv

  subroutine r_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_32

  subroutine r_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, real_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_32_inv

  subroutine r_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_42

  subroutine r_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, real_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_42_inv

  subroutine r_redist_23 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_23

  subroutine r_redist_34 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    real, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    real, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (out) :: to_here

    real, dimension (r%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       call shmem_get (ix, iy, 1, idp)
       if (r%from(idp)%ni > 0) then
          call shmem_ixget (real_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_redist_34

  subroutine r_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (out) :: to_here

    real, dimension (f%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_fill_2

  subroutine r_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (out) :: to_here

    real, dimension (f%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_fill_3

  subroutine r_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    real, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    real, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (out) :: to_here

    real, dimension (f%redist_buff_size) :: real_buff
    integer :: ip, idp

    pointer (ix, x)
    real x
    
    pointer (iy, y)
    real y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (real_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, real_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine r_fill_4

  subroutine i_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%from_low(1):), intent (in) :: from_here

    integer, dimension (r%to_low(1):,r%to_low(2):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_12

  subroutine i_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_22

  subroutine i_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, integer_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_22_inv

  subroutine i_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_32

  subroutine i_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, integer_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_32_inv

  subroutine i_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_42

  subroutine i_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, integer_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_42_inv

  subroutine i_redist_23 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_23

  subroutine i_redist_34 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    integer, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    integer, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (out) :: to_here

    integer, dimension (r%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_redist_34

  subroutine i_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (out) :: to_here

    integer, dimension (f%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_fill_2

  subroutine i_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (out) :: to_here

    integer, dimension (f%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_fill_3

  subroutine i_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    integer, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    integer, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (out) :: to_here

    integer, dimension (f%redist_buff_size) :: integer_buff
    integer :: ip, idp

    pointer (ix, x)
    integer x
    
    pointer (iy, y)
    integer y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (integer_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, integer_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine i_fill_4

  subroutine l_redist_12 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%from_low(1):), intent (in) :: from_here

    logical, dimension (r%to_low(1):,r%to_low(2):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_12

  subroutine l_redist_22 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_22

  subroutine l_redist_22_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, logical_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_22_inv

  subroutine l_redist_32 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_32

  subroutine l_redist_32_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, logical_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_32_inv

  subroutine l_redist_42 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_42

  subroutine l_redist_42_inv (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):), intent (in) :: from_here

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):, &
                        r%from_low(4):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%nj > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%to(idp)%j, r%to(idp)%nj, iproc)
          call shmem_ixput (x, logical_buff, r%from(idp)%j, r%from(idp)%nj, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_42_inv

  subroutine l_redist_23 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_23

  subroutine l_redist_34 (r, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: r

    logical, dimension (r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent (in) :: from_here

    logical, dimension (r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):), intent (out) :: to_here

    logical, dimension (r%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (r%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, r%from(idp)%i, r%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, r%to(idp)%i, r%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_redist_34

  subroutine l_fill_2 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):), intent (out) :: to_here

    logical, dimension (f%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_fill_2

  subroutine l_fill_3 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):), intent (out) :: to_here

    logical, dimension (f%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_fill_3

  subroutine l_fill_4 (f, from_here, to_here)

    use mp, only: iproc, nproc, barrier
    type (redist_type), intent (in) :: f

    logical, dimension (f%from_low(1):, &
                        f%from_low(2):, &
                        f%from_low(3):, &
                        f%from_low(4):), intent (in) :: from_here

    logical, dimension (f%to_low(1):, &
                        f%to_low(2):, &
                        f%to_low(3):, &
                        f%to_low(4):), intent (out) :: to_here

    logical, dimension (f%redist_buff_size) :: logical_buff
    integer :: ip, idp

    pointer (ix, x)
    logical x
    
    pointer (iy, y)
    logical y       

    iy = loc (to_here)    
    call barrier

    do ip = 0, nproc-1
       idp = mod(ip + iproc, nproc)

       if (f%from(idp)%ni > 0) then
          call shmem_get (ix, iy, 1, idp)
          call shmem_ixget (logical_buff, from_here, f%from(idp)%i, f%from(idp)%ni, iproc)
          call shmem_ixput (x, logical_buff, f%to(idp)%i, f%to(idp)%ni, idp)
       end if

    end do

    call barrier

  end subroutine l_fill_4

end module redistribute
