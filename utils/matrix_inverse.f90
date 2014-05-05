!>A module to hold different routines which can be used to invert a matrix
module mat_inv

!>Really we should provide an interface to different versions of these routines for
!different data types, doing things in-place or not and with assumed arrays and not etc.

!>If our code were needs to do lots of inverts then it may make sense to create a small
!timing routine which can be used to select the most efficient method for a given size on
!the current machine etc. This could then be used with a generic invert routine, the basic
!idea being:
!subroutine invert(a)
!:
!!If this is the first call then do some timing to pick best method
!if(first) call pick_best_method
!!This will set a module level variable identifying the best routine
!select case(best_method)
!case(0)
! call invert_type_1(a)
!:
!end select
!:
!end subroutine 
!Of course the best method may change with matrix size so we may need to calculate
!best_method for a number of sizes and then pick the one corresponding to the closest
!sized array to the one we have.

  !/Serial
!  public :: inverse_lapack
  public :: inverse_gj
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
    integer i, k

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
    integer i, k

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
