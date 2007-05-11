program test
  use mp, only: init_mp, finish_mp, proc0, nproc, broadcast, iproc
  use mp, only: send, receive, barrier

  implicit none
  real, dimension (:,:), allocatable :: rpass
  integer :: nsize, i, ip, j, ns2
  real :: a, b
  
  call init_mp

  nsize = 1000000
  ns2 = 100
  if (.not. proc0) then
     allocate (rpass(nsize,ns2))
  end if

  if (iproc == 1) then
     call random_number (rpass)
  end if

  i=0
  do 
     i=i+1
     do ip=1,nproc-1
        j = mod(i,ns2)+1
        if (iproc == ip) call send (rpass(:,j), mod(ip,nproc-1)+1)
        if (iproc == mod(ip,nproc-1)+1) then
           call receive (rpass(:,j), ip)
           rpass(:,j) = rpass(:,j) + real(ip)*rpass(:,j)/sum(rpass(:,j))
        end if
        call barrier
!        if (iproc == ip) call send (rpass(2897,j), 0)
!        if (proc0) call receive (a, ip)
!        call barrier
!        if (iproc == mod(ip,nproc-1)+1) call send (rpass(2897,j), 0)
!        if (proc0) call receive (b, mod(ip,nproc-1)+1)
!        call barrier
!        if (proc0) write(*,*) a,' = ',b,' => ',ip,' to ',mod(ip,nproc-1)+1,' ok'
        if (proc0) write(*,*) ip,' to ',mod(ip,nproc-1)+1,' ok'
     end do

     if (proc0) write(*,*) i
  end do

  call finish_mp

end program test
