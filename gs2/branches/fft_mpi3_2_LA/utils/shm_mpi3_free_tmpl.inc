 integer ierr
    logical done
    type(shm_node_pointers_t), pointer :: next=>null()

    done = .false.
    next => shm_ptr_head

    do while(associated(next))
        !write(*,*) 'shm_free ', shm_info%id, next%nd(shm_info%id), c_loc(a), next%se
       if (c_associated(next%nd(shm_info%id),c_loc(a))) then
          !write(*,*) 'shm_free win free'
          call mpi_win_free(next%win, ierr)
          done = .true.
          exit
       else
          next => next%next
       end if
    enddo

    if( .not. done) then 
       write(0,*) " Warning: shm_free failed to find a window !!!"
    endif