!>A module to hold different routines which can be used to invert a matrix
module mat_inv

  !/Serial
!  public :: inverse_lapack
  public :: inverse_gj
  public :: inverse_gj_t
  !/OpenMP
  public :: inverse_gj_omp
  !/MPI

!NOTE: No attempt has been made to optimise any of these routines yet.
 
  private

contains
  
  !>Serial
  ! !>A routine to invert a matrix using lapack routines
  ! !Will probably have to remove this before merging back as it introduces
  ! !new dependencies.
  ! subroutine inverse_lapack(a,n)
  !   integer,intent(in)::n
  !   complex,dimension(:,:),intent(inout)::a
  !   complex,dimension(:),allocatable::work
  !   integer :: lwork
  !   integer :: info
  !   integer,dimension(:),allocatable::ipiv

  !   !Setup
  !   lwork=n*n

  !   !Allocates
  !   allocate(work(lwork),ipiv(n))

  !   !LU decompose
  !   CALL ZGETRF( N, N, a, N, IPIV, INFO )
  !   if(info.ne.0)then
  !      write(6,'("Error during lapack LU factorization.")')
  !   endif

  !   !Invert
  !   CALL ZGETRI(N, a, N, IPIV, WORK, LWORK, INFO)
  !   if(info.ne.0)then
  !      write(6,'("Error during lapack invert")')
  !   endif

  !   !Cleanup
  !   deallocate(work,ipiv)
  ! end subroutine inverse_lapack

  !>Serial Gauss-Jordan elimination
  subroutine inverse_gj(a,n)
    implicit none 
    integer,intent(in):: n
    !Should we make this assumed shape?
    !<DD>Tagged
    complex,dimension(n,n),intent(inout) :: a

    !Is it better to make these allocatable?
    !<DD>Tagged
    complex :: tmp,fac
    integer i, j, k

    do i=1,n
       fac=1/a(i,i) !This would become inverse if done on blocks
       a(i,i)=1
       a(:,i)=a(:,i)*fac
       !NOTE:These loops don't shrink as i increases so load should
       !be balanced?
       do k=1,i-1
          tmp=a(i,k) 
          a(i,k)=0
          a(:,k)=a(:,k)-a(:,i)*tmp
       enddo
       do k=i+1,n
          tmp=a(i,k)
          a(i,k)=0
          a(:,k)=a(:,k)-a(:,i)*tmp
       enddo
    enddo
  end subroutine inverse_gj

  !>Serial Gauss-Jordan elimination with transposed addressing
  !This has much worse memory access than non-transposed.
  subroutine inverse_gj_t(a,n)
    implicit none 
    integer,intent(in):: n
    !Should we make this assumed shape?
    !<DD>Tagged
    complex,dimension(n,n),intent(inout) :: a

    !Is it better to make these allocatable?
    !<DD>Tagged
    complex :: tmp,fac
    integer i, j, k

    do i=1,n
       fac=1/a(i,i) !This would become inverse if done on blocks
       a(i,i)=1
       a(i,:)=a(i,:)*fac
       !NOTE:These loops don't shrink as i increases so load should
       !be balanced?
       do k=1,i-1
          tmp=a(k,i) 
          a(k,i)=0
          a(k,:)=a(k,:)-a(i,:)*tmp
       enddo
       do k=i+1,n
          tmp=a(k,i)
          a(k,i)=0
          a(k,:)=a(k,:)-a(i,:)*tmp
       enddo
    enddo
  end subroutine inverse_gj_t

  !>openmp Gauss-Jordan elimination
  subroutine inverse_gj_omp(a,n)
    implicit none 
    integer,intent(in):: n
    !Should we make this assumed shape?
    !<DD>Tagged
    complex,dimension(n,n),intent(inout) :: a

    !Is it better to make these allocatable?
    !<DD>Tagged
    complex :: tmp,fac
    integer i, j, k

    !Can't omp at top level currently due to data dependencies
    do i=1,n
       fac=1/a(i,i) !This would become inverse if done on blocks
       a(i,i)=1
       a(:,i)=a(:,i)*fac
       !NOTE: Slight difference to serial version to avoid
       !having two separate parallel regions
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(tmp)
       do k=1,n
          if(k.eq.i)cycle
          tmp=a(i,k)
          a(i,k)=0
          a(:,k)=a(:,k)-a(:,i)*tmp
       enddo
       !$OMP END PARALLEL DO
    enddo

  end subroutine inverse_gj_omp
end module mat_inv
