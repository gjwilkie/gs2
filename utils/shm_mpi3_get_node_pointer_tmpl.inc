    type(shm_node_pointers_t), pointer :: next => null()

    integer i, jd, ierr
    integer, allocatable :: ashape(:)
    logical found

    if ( present(a1)) a1 =>null()
    if ( present(a2)) a2 =>null()
    if ( present(a3)) a3 =>null()

    ! id = -1 is for pointer from node to whole array
    ! a shm symbolic name must be used 

    if ( id == -1) then 
       jd = 0
    else
       jd = id
    endif

     next => shm_ptr_head
     do while (associated(next))
        if (present(tag)) then
           found = tag == next%tag
        else
           found = c_associated(next%nd(shm_info%id), c_loc(pin))
        endif
        if(found) then
           if ( present(a1) .and. next%ndim /= 1 .or. &
                present(a2) .and. next%ndim /= 2 .or. &
                present(a3) .and. next%ndim /= 3 ) then
              call error_abort("shm_get_node_pointer: inconsistent out pointer and dimension")
           endif
           allocate(ashape(next%ndim))
           do i = 1, next%ndim
              ashape(i) =next%se(2*i,jd) - next%se(2*i-1,jd) +1
           enddo

           if ( id == -1) then 
              i = next%ndim
              ashape(i) = next%se(2*i,shm_info%size -1) - next%se(2*i-1,0) +1
           endif
           !write(0,*) 'ashape', shm_info%id, ashape, next%se
           
           if ( id == -1 ) then
              select case (next%ndim)
              case(1)
                 call c_f_pointer(next%nd(jd), a1, ashape)
              case(2)
                 call c_f_pointer(next%nd(jd), a2, ashape)
              case(3)
                 call c_f_pointer(next%nd(jd), a3, ashape)
              case default
                 call error_abort(" shm_get_node_pointer : wrong ndim")
              end select
           else
              select case (next%ndim)
              case(1)
                 call c_f_pointer(next%nd(jd), a1, ashape)
                 a1 => remap_bounds(next%se(1,jd), a1)
              case(2)
                 call c_f_pointer(next%nd(jd), a2, ashape)
                 a2 =>  remap_bounds(next%se(1,jd), next%se(3,jd), a2)
              case(3)
                 call c_f_pointer(next%nd(jd), a3, ashape)
                 a3 =>  remap_bounds(next%se(1,jd), next%se(3,jd), next%se(5,jd), a3)
              case default
                 call error_abort(" shm_get_node_pointer : wrong ndim")
              end select
           endif
           exit
        else
           next => next%next
        endif
     end do

! ATTENTION this works only pin is allocated, i..e. all e >= s

!!$    pointer_array: do j = 1, size(shm_info%g_lo_ptr)
!!$       do i = 0, shm_info%size -1
!!$          s = shm_info%g_lo_se(1, i)
!!$          e = shm_info%g_lo_se(2, i)
!!$          aux => shm_info%g_lo_ptr(j)%p(:,:, s:e)
!!$         !print*,'get_node_pointer', iproc, i,j,s,e
!!$          if ( associated(aux, pin)) then
!!$             get_node_pointer => shm_info%g_lo_ptr(j)%p
!!$             !print'(a,7(i5,x))','get_node_pointer', iproc, lbound(get_node_pointer), ubound(get_node_pointer)
!!$             exit pointer_array
!!$          endif
!!$       end do
!!$ enddo pointer_array

 if ( present(a1) .and. .not. associated(a1) .or. &
      present(a2) .and. .not. associated(a2) .or. &
      present(a3) .and. .not. associated(a3) ) then
    call error_abort('shm_get_node_pointer fail to find a node pointer')
 endif
