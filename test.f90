program test
  use mp, only: init_mp, finish_mp, proc0, nproc, broadcast, iproc
  use mp, only: send, receive, barrier

  implicit none
  integer, dimension (:,:), allocatable :: npass
  integer :: nsize, i, ip, ns2, j
  
  call init_mp

  nsize = 500000
  ns2 = 450
  if (.not. proc0) then
     allocate (npass(nsize, ns2))
  end if

  if (iproc == 1) then
     do j=1,ns2
        npass(:,j) = (/ (i-j, i=1,nsize) /)
     end do
  end if

  i=0
  do 
     i=i+1
     do ip=1,nproc-1
        if (iproc == ip) call send (npass(:,mod(i-1,ns2)+1), mod(ip,nproc-1)+1)
        if (iproc == mod(ip,nproc-1)+1) call receive (npass(:,mod(i-1,ns2)+1), ip)
        call barrier
        if (proc0) write(*,*) ip,' to ',mod(ip,nproc-1)+1,' worked'
     end do
     if (proc0) write(*,*) i
  end do

  if (iproc == nproc-1) then
     write (*,*) iproc, npass(nsize,ns2)
  end if

  call finish_mp

end program test
