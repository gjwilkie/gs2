module layouts

  implicit none
  private

  type :: x_layout_type 
     integer :: nblock
     integer :: iproc
     integer :: llim_world, ulim_world
     integer :: llim_proc, ulim_proc, ulim_alloc
     integer :: nx, ny, nz
  end type x_layout_type

  type :: y_layout_type 
     integer :: nblock
     integer :: iproc
     integer :: llim_world, ulim_world
     integer :: llim_proc, ulim_proc, ulim_alloc
     integer :: nx, ny, nz
  end type y_layout_type

  type :: layout_type
     type (x_layout_type) :: x
     type (y_layout_type) :: y
     integer :: nx, ny, nz
  end type layout_type

  type (layout_type) :: lo

  interface proc_id
     module procedure proc_id_x, proc_id_y
  end interface

  interface idx
     module procedure idx_x, idx_y
  end interface

  interface idx_local
     module procedure idx_local_x
     module procedure idx_local_y
  end interface

contains

  subroutine init_layouts (nx, ny, nzp)
    
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: nx, ny
    integer, intent (in), optional :: nzp
    logical, save :: initialized = .false.
    integer :: nyblock, nxblock, nz

    if (initialized) return
    initialized = .true.

    if (present(nzp)) then
       nz = nzp
    else
       nz = 0
    end if

    lo%nx = nx
    lo%ny = ny
    lo%nz = nz

    nyblock = nx/nproc

    lo % y % iproc = iproc
    lo % y % nblock = nyblock
    lo % y % llim_world = 1
    lo % y % ulim_world = nx
    lo % y % llim_proc = iproc*nyblock + 1
    lo % y % ulim_proc = min (nx, lo % y % llim_proc + nyblock - 1)
    lo % y % ulim_alloc = max (lo % y % llim_proc, lo % y % ulim_proc)
    lo % y % nx = nx
    lo % y % ny = ny
    lo % y % nz = nz


    nxblock = ny/nproc

    lo % x % iproc = iproc
    lo % x % nblock = nxblock
    lo % x % llim_world = 1
    lo % x % ulim_world = ny
    lo % x % llim_proc = iproc*nxblock + 1
    lo % x % ulim_proc = min (ny, lo % x % llim_proc + nxblock - 1)
    lo % x % ulim_alloc = max (lo % x % llim_proc, lo % x % ulim_proc)
    lo % x % nx = nx
    lo % x % ny = ny
    lo % x % nz = nz


  end subroutine init_layouts

  pure function proc_id_x (lo, i)
    implicit none
    integer :: proc_id_x
    type (x_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_x = (i-1)/lo%nblock

  end function proc_id_x

  pure function idx_x (lo, j, k)
    implicit none
    integer :: idx_x
    type (x_layout_type), intent (in) :: lo
    integer, intent (in) :: j
    integer, intent (in), optional :: k
    
    if (present(k)) then
       idx_x = k + lo%nz*(j-1)
    else
       idx_x = j
    end if

  end function idx_x

  pure function idx_local_x (lo, j)
    implicit none
    logical :: idx_local_x
    type (x_layout_type), intent (in) :: lo
    integer, intent (in) :: j
    
    idx_local_x = lo%iproc == proc_id (lo, j)

  end function idx_local_x

  pure function proc_id_y (lo, i)
    implicit none
    integer :: proc_id_y
    type (y_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_y = (i-1)/lo%nblock

  end function proc_id_y

  pure function idx_y (lo, i)
    implicit none
    integer :: idx_y
    type (y_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
    idx_y = i

  end function idx_y

  pure function idx_local_y (lo, i)
    implicit none
    logical :: idx_local_y
    type (y_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
    idx_local_y = lo%iproc == proc_id (lo, i)

  end function idx_local_y

end module layouts
