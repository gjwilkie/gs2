module centering

  implicit none

  public :: get_cell_value, invert_matrix

  interface get_cell_value
     module procedure get_cell_value_1d
     module procedure get_cell_value_2d
     module procedure get_cell_value_1d_complex
     module procedure get_cell_value_2d_complex
  end interface

contains

  ! takes an array of grid values and calculates a cell
  ! value, with the location in the cell determined by impfac
  ! as cell = impfac*grid_{i+1} + (1.0-impfac)*grid_{i};
  ! e.g., if impfac=0.5, then cell is the value at the cell center
  subroutine get_cell_value_1d (impfac, grid, cell, l)

    implicit none

    real, intent (in) :: impfac
    real, dimension (l:), intent (in) :: grid
    real, dimension (l:), intent (out) :: cell
    integer, intent (in) :: l

    integer :: u

    u = -l-1

    cell(l:u) = impfac*grid(l+1:) + (1.0-impfac)*grid(:u)

  end subroutine get_cell_value_1d

  ! takes an array of grid values in 2d and calculates cell
  ! values, with the location in the cell determined by th_imp and vp_imp
  subroutine get_cell_value_2d (th_imp, vp_imp, grid, cell, l1, l2)

    implicit none

    real, intent (in) :: th_imp, vp_imp
    real, dimension (l1:,l2:), intent (in) :: grid
    real, dimension (l1:,l2:), intent (out) :: cell
    integer, intent (in) :: l1, l2

    integer :: u1, u2, ig, iv

    ! upper limits for first and second indices of grid and cell
    u1 = -l1
    u2 = -l2

    ! first do centering in vpar
    ! for theta<0
    do ig = l1, -1 
       call get_cell_value (vp_imp,grid(ig,:),cell(ig,:),l2)
    end do
    ! for theta>0
    do ig = 0, u1-1
       call get_cell_value (1.0-vp_imp,grid(ig,:),cell(ig,:),l2)
    end do
    
    ! next do centering in theta
    ! for vpa<0
    do iv = l2, -1
       call get_cell_value (1.0-th_imp,grid(:,iv),cell(:,iv),l1)
    end do
    ! for vpa>0
    do iv = 0, u2-1
       call get_cell_value (th_imp,grid(:,iv),cell(:,iv),l1)
    end do
    
  end subroutine get_cell_value_2d

  ! takes an array of grid values and calculates a cell
  ! value, with the location in the cell determined by impfac
  ! as cell = impfac*grid_{i+1} + (1.0-impfac)*grid_{i};
  ! e.g., if impfac=0.5, then cell is the value at the cell center
  subroutine get_cell_value_1d_complex (impfac, grid, cell, l)

    implicit none

    real, intent (in) :: impfac
    complex, dimension (l:), intent (in) :: grid
    complex, dimension (l:), intent (out) :: cell
    integer, intent (in) :: l

    integer :: u

    u = -l-1

    cell(l:u) = impfac*grid(l+1:) + (1.0-impfac)*grid(:u)

  end subroutine get_cell_value_1d_complex

  ! takes an array of grid values in 2d and calculates cell
  ! values, with the location in the cell determined by th_imp and vp_imp
  subroutine get_cell_value_2d_complex (th_imp, vp_imp, grid, cell, l1, l2)

    implicit none

    real, intent (in) :: th_imp, vp_imp
    complex, dimension (l1:,l2:), intent (in) :: grid
    complex, dimension (l1:,l2:), intent (out) :: cell
    integer, intent (in) :: l1, l2

    integer :: u1, u2, ig, iv

    ! upper limits for first and second indices of grid and cell
    u1 = -l1
    u2 = -l2

    ! first do centering in vpar
    ! for theta<0
    do ig = l1, -1 
       call get_cell_value (vp_imp,grid(ig,:),cell(ig,:),l2)
    end do
    ! for theta>0
    do ig = 0, u1-1
       call get_cell_value (1.0-vp_imp,grid(ig,:),cell(ig,:),l2)
    end do
    
    ! next do centering in theta
    ! for vpa<0
    do iv = l2, -1
       call get_cell_value (1.0-th_imp,grid(:,iv),cell(:,iv),l1)
    end do
    ! for vpa>0
    do iv = 0, u2-1
       call get_cell_value (th_imp,grid(:,iv),cell(:,iv),l1)
    end do
    
  end subroutine get_cell_value_2d_complex

  subroutine invert_matrix (mat)

    implicit none

    real, dimension (:,:), intent (in out) :: mat

    integer :: i, n
    real :: sgn
    integer, dimension (:), allocatable :: indx
    real, dimension (:), allocatable :: rhs
    real, dimension (:,:), allocatable :: lu

    n = size(mat,1)

    allocate (indx(n))
    allocate (rhs(n))
    allocate (lu(n,n))

    lu = mat

    ! get LU decomposition of lu=mat and overwrite it
    call ludcmp (lu, indx, sgn)

    do i = 1, size(mat,1)
       rhs = 0.0 ; rhs(i) = 1.0
       ! get solution to Ax=b
       ! with A = mat, b = input rhs, x = output rhs
       call lubksb (lu, indx, rhs)
       ! construct the inverse of A from the solution columns
       mat(:,i) = rhs
    end do

    deallocate (indx, rhs, lu)

  contains

    subroutine ludcmp (a, indx, d)
      
      implicit none
      
      real, dimension (:,:), intent (in out) :: a
      integer, dimension (:), intent (out) :: indx
      real, intent (out) :: d
      
      real, dimension (size(a,1)) :: vv
      real, parameter :: tiny = 1.0e-20
      integer :: j, n, imax
      
      n = assert_eq(size(a,1),size(a,2), size(indx), 'ludcmp')
      d = 1.0
      vv = maxval(abs(a),dim=2)
      if (any(vv == 0.0)) then
         write (*,*) 'singular matrix in ludcmp'
      end if
      vv = 1.0/vv
      do j = 1, n
         imax = (j-1) + imaxloc(vv(j:n)*abs(a(j:n,j)))
         if (j /= imax) then
            call swap(a(imax,:),a(j,:))
            d = -d
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (a(j,j) == 0.0) a(j,j) = tiny
         a(j+1:n,j) = a(j+1:n,j) / a(j,j)
         a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - spread(a(j+1:n,j),dim=2,ncopies=n-j) * &
              spread(a(j,j+1:n),dim=1,ncopies=n-j)
      end do
      
    end subroutine ludcmp

    subroutine lubksb (a, indx, b)
      
      implicit none
      
      real, dimension (:,:), intent (in) :: a
      integer, dimension (:), intent (in) :: indx
      real, dimension (:), intent (in out) :: b
      
      integer :: i, n, ii, ll
      real :: summ
      
      n = assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
      ii = 0
      do i = 1, n
         ll = indx(i)
         summ = b(ll)
         b(ll) = b(i)
         if (ii /= 0) then
            summ = summ - dot_product(a(i,ii:i-1),b(ii:i-1))
         else if (summ /= 0.0) then
            ii = i
         end if
         b(i) = summ
      end do
      do i = n, 1, -1
         b(i) = (b(i) - dot_product(a(i,i+1:n),b(i+1:n))) / a(i,i)
      end do
      
    end subroutine lubksb
    
    subroutine swap (a,b)
      real, dimension (:), intent (in out) :: a, b
      real, dimension (size(a)) :: dum
      dum=a
      a=b
      b=dum
    end subroutine swap
    
    function assert_eq (n1,n2,n3,string)
      character (*), intent (in) :: string
      integer, intent (in) :: n1,n2,n3
      integer :: assert_eq
      if (n1 == n2 .and. n2 == n3) then
         assert_eq=n1
      else
         write (*,*) 'error: an assert_eq failed with this tag:', &
              string
         stop 'program terminated by assert_eq3'
      end if
    end function assert_eq
    
    function imaxloc (arr)
      real, dimension (:), intent (in) :: arr
      integer :: imaxloc
      integer, dimension (1) :: imax
      imax = maxloc(arr(:))
      imaxloc = imax(1)
    end function imaxloc
    
  end subroutine invert_matrix

end module centering
