    type(shm_node_pointers_t), pointer :: next => null()
    integer i, ierr
    integer, allocatable :: ashape(:)
    logical found

    found=.false.

    next => shm_ptr_head
     do while (associated(next))
        !write(0,*) shm_info%id, counter, next%nd(shm_info%id), c_loc(pin)
        found = c_associated(next%nd(0), c_loc(a))
        if(found) then
           !write(0,*) '## -> got win flush'
           call MPI_WIN_FENCE(0, next%win, ierr)
           exit
        else
           next => next%next
        endif
     end do

     if (.not. found) then
        call error_abort("shm_flush: no shared block found for input array")
     endif
