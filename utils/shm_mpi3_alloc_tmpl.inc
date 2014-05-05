  !  not working with the current installation of mpi3 on mac  
  !call MPI_SIZEOF(x, disp, ierr)
  disp = storage_size(x)/8 ! this is 2008 standard

    allocate (ashape(ndim))
    do i = 1, ndim
       ashape(i) = lubd(2*i) - lubd(2*i-1) +1
    enddo

    asize = disp * product(ashape)
    if (debug) write(0,*) "alloc size ", shm_info%id, asize	

    call mpi_win_allocate_shared(asize, disp, info_noncontig, shm_info%comm, base_ptr, win, ierr)

    select case (ndim)
    case(1)
       call c_f_pointer(base_ptr, a1, ashape)
       a1 => remap_bounds(lubd(1), a1)
    case(2)
       call c_f_pointer(base_ptr, a2, ashape)
       a2 =>  remap_bounds(lubd(1), lubd(3), a2)
    case(3)
       call c_f_pointer(base_ptr, a3, ashape)
       a3 =>  remap_bounds(lubd(1), lubd(3), lubd(5), a3)
    case default
        call error_abort("shm_alloc: wrong ndim")
     end select

    ! store shared memory info in the linked list
    if (.not.associated(shm_pointers)) then 
       allocate(shm_pointers)
       shm_ptr_head => shm_pointers
    else
       allocate(shm_pointers%next)
       shm_pointers => shm_pointers%next
    end if
    shm_pointers%next => null()

    shm_pointers%win = win
    id = shm_info%id
    shm_pointers%ndim = ndim
    allocate(shm_pointers%se(2*ndim, 0:shm_info%size-1))
    allocate(shm_pointers%nd(0 : shm_info%size - 1))
    shm_pointers%nd(id) = base_ptr 
    ! what if base_ptr is c_null_ptr ???
    counter = counter + 1
    shm_pointers%tag=counter
    if (present(tag)) then
       tag=counter
    endif
    if(present(label)) then
       shm_pointers%label=label
    endif

    do i = 0, shm_info%size -1
       if (i == id) cycle
       call mpi_win_shared_query(win, i, asize, disp, base_ptr, ierr)
       shm_pointers%nd(i) = base_ptr
    end do

    call mpi_allgather(lubd,2*ndim,MPI_INTEGER,shm_pointers%se, &
         2*ndim,MPI_INTEGER,shm_info%comm,ierr)
