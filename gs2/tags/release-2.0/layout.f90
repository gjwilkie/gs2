module layout
  implicit none

  integer, parameter :: maxdim = 7

  type layout_type
     integer :: nprocs, nblocksize, ndim
     integer, dimension (maxdim) :: llim, ulim, size, skip
  end type layout_type

  interface make_new_layout
     module_procedure make_new_layout_from_sizes
     module_procedure make_new_layout_from_nprocs_and_sizes
     module_procedure make_new_layout_from_bounds
     module_procedure make_new_layout_from_nprocs_and_bounds
  end interface

contains

  subroutine init_layout (lo)
    implicit none
    type (layout_type), intent (in out) :: lo
    integer :: i

    lo%size = lo%ulim - lo%llim
    lo%skip(1) = 1
    do i = 2, lo%ndim
       lo%skip(i) = lo%skip(i-1)*lo%size(i-1)
    end do
    lo%nblocksize = (layout_size(lo) - 1)/lo%nprocs + 1
  end subroutine init_layout

  subroutine make_new_layout_from_nprocs_and_sizes (lo, nprocs, sizes)
    implicit none
    type (layout_type), intent (out) :: lo
    integer, intent (in) :: nprocs
    integer, dimension (:), intent (in) :: sizes

    lo%nprocs = nprocs
    lo%ndim = size(sizes)
    lo%llim = 1
    lo%ulim = 0
    lo%ulim(:lo%ndim) = sizes
    call init_layout (lo)
  end subroutine make_new_layout_from_nprocs_and_sizes

  subroutine make_new_layout_from_sizes (lo, sizes)
    implicit none
    type (layout_type), intent (out) :: lo
    integer, dimension (:), intent (in) :: sizes
    include 'mpif.h'
    integer :: nprocs, ierror

    call mpi_comm_size (mpi_comm_world, nprocs, ierror)
    call make_new_layout_from_nprocs_and_sizes (lo, nprocs, sizes)
  end subroutine make_new_layout_from_sizes

  subroutine make_new_layout_from_bounds (lo, llim, ulim)
    implicit none
    type (layout_type), intent (out) :: lo
    integer, dimension (:), intent (in) :: llim, ulim
    include 'mpif.h'
    integer :: nprocs, ierror

    call mpi_comm_size (mpi_comm_world, nprocs, ierror)
    call make_new_layout_from_nprocs_and_bounds (lo, nprocs, llim, ulim)
  end subroutine make_new_layout_from_sizes

  subroutine make_new_layout_from_nprocs_and_bounds (lo, nprocs, llim, ulim)
    implicit none
    type (layout_type), intent (out) :: lo
    integer, intent (in) :: nprocs
    integer, dimension (:), intent (in) :: llim, ulim

    lo%nprocs = nprocs
    lo%ndim = size(llim)
    lo%llim = 1
    lo%ulim = 0
    lo%llim(:lo%ndim) = llim
    lo%ulim(:lo%ndim) = ulim
    call init_layout (lo)
  end subroutine make_new_layout_from_sizes

  pure function layout_size (lo)
    integer :: layout_size
    type (layout_type), intent (in) :: lo
    layout_size = product(lo%size(:lo%ndim))
  end function layout_size

  pure function layout_linear_index (idx, lo)
    integer :: layout_index
    integer, dimension (:), intent (in) :: idx
    type (layout_type), intent (in) :: lo

    layout_index = sum(idx(:lo%ndim)*lo%skip(lo%ndim))
  end function layout_linear_index

  pure function layout_home (idx, lo)
    integer :: layout_home
    integer, dimension (:), intent (in) :: idx
    type (layout_type), intent (in) :: lo

    layout_home = layout_linear_index(idx,lo)/lo%nblocksize
  end function layout_home

  pure function layout_index (linear, dim, lo)
    integer :: layout_index
    integer, intent (in) :: linear, dim
    type (layout_type), intent (in) :: lo

    layout_index = linear/lo%skip(dim) + lo%llim(dim)
  end function layout_index

end module layout
