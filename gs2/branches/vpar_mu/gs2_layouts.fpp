# include "define.inc"

module gs2_layouts

  use layouts_type, only: g_layout_type

  implicit none
  private

  public :: layout, finish_layouts, local_field_solve
  
  public :: init_dist_fn_layouts, init_gs2_layouts
  public :: g_lo

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_jfields_layouts
  public :: jf_lo, jf_layout_type
  public :: mj, ij, dj

  public :: init_x_transform_layouts, init_y_transform_layouts
  public :: xxf_lo, xxf_layout_type, yxf_lo, yxf_layout_type
  public :: gidx2xxfidx, xxfidx2yxfidx, yxfidx2xxfidx
  public :: xxf_ky_is_zero

  public :: ig_idx, ik_idx, it_idx, iv_idx, imu_idx, is_idx, if_idx
  public :: im_idx, in_idx, ij_idx, ifield_idx
  public :: idx, proc_id, idx_local

  logical :: initialized_x_transform = .false.
  logical :: initialized_y_transform = .false.

  logical :: local_field_solve
  character (len=4) :: layout
  logical :: exist

  type (g_layout_type) :: g_lo

  type :: f_layout_type
     integer :: iproc
     integer :: nidx, nfield, ntgrid, nindex, naky, ntheta0, M_class, N_class, i_class
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
     integer, dimension(:,:), pointer :: ik => null ()
     integer, dimension(:,:), pointer :: it => null ()
  end type f_layout_type

  type (f_layout_type), dimension(:), allocatable :: f_lo

  type :: jf_layout_type
     integer :: iproc
     integer :: nindex, naky, ntheta0, ntgrid, nfield
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type jf_layout_type

  type (jf_layout_type) :: jf_lo

  integer, dimension(:), allocatable :: ij, mj
  integer, dimension(:,:), allocatable :: dj

  type :: xxf_layout_type
     integer :: iproc
     integer :: ntgrid, nvgrid, naky, ntheta0, nx, nadd, nmu, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
  end type xxf_layout_type

  type (xxf_layout_type) :: xxf_lo

  type :: yxf_layout_type
     integer :: iproc
     integer :: ntgrid, nvgrid, naky, ny, ntheta0, nx, nmu, nspec, ntgridtotal, nvgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
     integer :: it_ord, ig_ord, iv_ord, imu_ord, is_ord
     integer :: it_comp, ig_comp, iv_comp, imu_comp, is_comp
     integer, dimension (5) :: compound_count
  end type yxf_layout_type

  type (yxf_layout_type) :: yxf_lo

  interface ij_idx
     module procedure ij_idx_f
     module procedure ij_idx_jf
  end interface

  interface im_idx
     module procedure im_idx_f
  end interface

  interface in_idx
     module procedure in_idx_f
  end interface

  interface ifield_idx
     module procedure ifield_idx_f
  end interface

  interface if_idx
     module procedure if_idx_f
     module procedure if_idx_jf
  end interface

  interface ig_idx
     module procedure ig_idx_xxf
     module procedure ig_idx_yxf
     module procedure ig_idx_f
  end interface

  interface iv_idx
     module procedure iv_idx_xxf
     module procedure iv_idx_yxf
  end interface

  interface ik_idx
     module procedure ik_idx_g
     module procedure ik_idx_jf
     module procedure ik_idx_xxf
  end interface

  interface it_idx
     module procedure it_idx_jf
     module procedure it_idx_yxf
  end interface

  interface imu_idx
     module procedure imu_idx_g
     module procedure imu_idx_xxf
     module procedure imu_idx_yxf
  end interface

  interface is_idx
     module procedure is_idx_g
     module procedure is_idx_xxf
     module procedure is_idx_yxf
  end interface

  interface proc_id
     module procedure proc_id_g
     module procedure proc_id_f
     module procedure proc_id_jf
     module procedure proc_id_xxf
     module procedure proc_id_yxf
  end interface

  interface idx
     module procedure idx_g
     module procedure idx_f
     module procedure idx_jf
     module procedure idx_xxf
     module procedure idx_yxf
  end interface

  interface idx_local
     module procedure idx_local_g,      ig_local_g
     module procedure idx_local_f,      ig_local_f
     module procedure idx_local_jf,     ig_local_jf
     module procedure idx_local_xxf,    ig_local_xxf
     module procedure idx_local_yxf,    ig_local_yxf
  end interface

contains

  subroutine init_gs2_layouts
    
    use mp, only: proc0
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    if (proc0) call read_parameters
    call broadcast_results

  end subroutine init_gs2_layouts

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist, error_unit
    implicit none
    integer :: in_file
    namelist /layouts_knobs/ layout, local_field_solve
    local_field_solve = .false.
    layout = 'xyms'
    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)
    if (layout.ne.'xyms' .and. layout.ne.'xmys' .and. layout.ne.'yxms' &
         .and. layout.ne.'mxys' .and. layout.ne.'ymxs' .and. layout.ne.'myxs') &
    then
       write(6,*) "gs2_layouts: read_parameters finds illegal layout=",layout," =>stop"
       stop
    endif

  end subroutine read_parameters
    
  subroutine broadcast_results
    use mp, only: broadcast
    implicit none

    call broadcast (layout)
    call broadcast (local_field_solve)

  end subroutine broadcast_results

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dist_fn_layouts &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec)

    use mp, only: iproc, nproc
    use file_utils, only: error_unit

    implicit none

    integer, intent (in) :: ntgrid, naky, ntheta0, nvgrid, nmu, nspec
    logical, save :: initialized = .false.

# ifdef USE_C_INDEX
    integer :: ierr
    interface
       function init_indices_glo_c (layout)
         integer :: init_indices_glo_c
         character(*) :: layout
       end function init_indices_glo_c
    end interface
# endif

    if (initialized) return
    initialized = .true.
   
    g_lo%iproc = iproc
    g_lo%ntgrid = ntgrid
    g_lo%naky = naky
    g_lo%ntheta0 = ntheta0
    g_lo%nvgrid = nvgrid
    g_lo%nmu = nmu
    g_lo%nspec = nspec
    g_lo%llim_world = 0
    g_lo%ulim_world = naky*nmu*nspec - 1
      
    g_lo%blocksize = g_lo%ulim_world/nproc + 1
    g_lo%llim_proc = g_lo%blocksize*iproc
    g_lo%ulim_proc = min(g_lo%ulim_world, g_lo%llim_proc + g_lo%blocksize - 1)
    g_lo%ulim_alloc = max(g_lo%llim_proc, g_lo%ulim_proc)

# ifdef USE_C_INDEX
    ierr = init_indices_glo_c (layout)
    if (ierr /= 0) &
         & write (error_unit(),*) 'ERROR: layout not found: ', trim(layout)
# endif

  end subroutine init_dist_fn_layouts

# ifdef USE_C_INDEX
  function is_idx_g (lo, i)
# else
  elemental function is_idx_g (lo, i)
# endif

    implicit none
    integer :: is_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

# ifdef USE_C_INDEX
    interface
       function is_idx_g_c (lo,num)
         use layouts_type, only: g_layout_type
         integer :: is_idx_g_c
         type (g_layout_type) :: lo
         integer :: num
       end function is_idx_g_c
    end interface
    is_idx_g = is_idx_g_c (lo,i)
# else
    ! TT: the order of the division does not matter, so no need for branching
    is_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%nmu, lo%nspec)
# endif

  end function is_idx_g

# ifdef USE_C_INDEX
  function imu_idx_g (lo, i)
# else
  elemental function imu_idx_g (lo, i)
# endif

    implicit none

    integer :: imu_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

# ifdef USE_C_INDEX
    interface
       function imu_idx_g_c (lo,num)
         use layouts_type, only: g_layout_type
         integer :: imu_idx_g_c
         type (g_layout_type) :: lo
         integer :: num
       end function imu_idx_g_c
    end interface
    imu_idx_g = imu_idx_g_c (lo,i)
# else
    select case (layout)
    case ('yxms')
       imu_idx_g = 1 + mod((i - lo%llim_world)/lo%naky, lo%nmu)
    case ('xyms')
       imu_idx_g = 1 + mod((i - lo%llim_world)/lo%naky, lo%nmu)
    case ('mxys')
       imu_idx_g = 1 + mod((i - lo%llim_world), lo%nmu)
    case ('myxs')
       imu_idx_g = 1 + mod((i - lo%llim_world), lo%nmu)
    end select
# endif

  end function imu_idx_g

! TT>
# ifdef USE_C_INDEX
  function ik_idx_g (lo, i)
# else
  elemental function ik_idx_g (lo, i)
# endif
! <TT
    implicit none
    integer :: ik_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

# ifdef USE_C_INDEX
    interface
       function ik_idx_g_c (lo,num)
         use layouts_type, only: g_layout_type
         integer :: ik_idx_g_c
         type (g_layout_type) :: lo
         integer :: num
       end function ik_idx_g_c
    end interface
    ik_idx_g = ik_idx_g_c (lo,i)
# else
    select case (layout)
    case ('xyms')
       ik_idx_g = 1 + mod((i - lo%llim_world), lo%naky)
    case ('xmys')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nmu, lo%naky)
    case ('yxms')
       ik_idx_g = 1 + mod(i - lo%llim_world, lo%naky)
    case ('mxys')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nmu, lo%naky)
    case ('myxs')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nmu, lo%naky)
    case ('ymxs')
       ik_idx_g = 1 + mod(i - lo%llim_world, lo%naky)
    end select
# endif

  end function ik_idx_g

  elemental function proc_id_g (lo, i)
    implicit none
    integer :: proc_id_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_g = i/lo%blocksize

  end function proc_id_g

# ifdef USE_C_INDEX
  function idx_g (lo, ik, imu, is)
# else
  elemental function idx_g (lo, ik, imu, is)
# endif

    implicit none

    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, imu, is

# ifdef USE_C_INDEX
    interface
       function idx_g_c (lo,ik,imu,is)
         use layouts_type, only: g_layout_type
         integer :: idx_g_c
         type (g_layout_type) :: lo
         integer :: ik,imu,is
       end function idx_g_c
    end interface
    idx_g = idx_g_c (lo, ik, imu, is)
# else
    select case (layout)
    case ('xyms')
       idx_g = ik-1 + lo%naky*(imu-1 + lo%nmu*(is-1))
    case ('xmys')
       idx_g = imu-1 + lo%nmu*(ik-1 + lo%naky*(is-1))
    case ('yxms')
       idx_g = ik-1 + lo%naky*(imu-1 + lo%nmu*(is-1))
    case ('mxys')
       idx_g = imu-1 + lo%nmu*(ik-1 + lo%naky*(is-1))
    case ('myxs')
       idx_g = imu-1 + lo%nmu*(ik-1 + lo%naky*(is-1))
    case ('ymxs')
       idx_g = ik-1 + lo%naky*(imu-1 + lo%nmu*(is-1))
    end select
# endif

  end function idx_g

! TT>
# ifdef USE_C_INDEX
  function idx_local_g (lo, ik, imu, is)
# else
  elemental function idx_local_g (lo, ik, imu, is)
# endif
! <TT
    implicit none
    logical :: idx_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, imu, is

    idx_local_g = idx_local(lo, idx(lo, ik, imu, is))
  end function idx_local_g

  elemental function ig_local_g (lo, ig)
    implicit none
    logical :: ig_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_g = lo%iproc == proc_id(lo, ig)
  end function ig_local_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Field layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! M_class(i) = number of i_th supercell
! N_class(i) = size of i_th supercell
! i_class = number of classes of (different sized) supercells.  

  subroutine init_fields_layouts (nfield, nindex, naky, ntheta0, &
       M_class, N_class, i_class)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: nfield, nindex, naky, ntheta0
    integer, dimension(:), intent (in) :: M_class, N_class
    integer, intent (in) :: i_class
    logical, save :: initialized = .false.
    integer :: i, utmp, btmp

    if (initialized) return
    initialized = .true.
    
    allocate (f_lo(i_class))
    do i = 1, i_class
       allocate (f_lo(i)%ik(M_class(i), N_class(i)))
       allocate (f_lo(i)%it(M_class(i), N_class(i)))
    end do

    do i = 1, i_class
       f_lo(i)%iproc = iproc
       f_lo(i)%nfield = nfield
       f_lo(i)%nidx = nindex                 ! cell size
       f_lo(i)%ntgrid = (nindex/nfield-1)/2
       f_lo(i)%nindex = nindex*N_class(i)    ! supercell size
       f_lo(i)%naky = naky
       f_lo(i)%ntheta0 = ntheta0
       f_lo(i)%M_class = M_class(i)
       f_lo(i)%N_class = N_class(i)
       f_lo(i)%i_class = i_class
       f_lo(i)%llim_world = 0
       f_lo(i)%ulim_world = f_lo(i)%nindex*M_class(i) - 1
       if (local_field_solve) then    ! guarantee local operations for matrix inversion
          utmp = M_class(i) - 1
          btmp = utmp/nproc + 1
          f_lo(i)%blocksize = f_lo(i)%nindex*btmp
       else
          f_lo(i)%blocksize = f_lo(i)%ulim_world/nproc + 1
       end if
       f_lo(i)%llim_proc = f_lo(i)%blocksize*iproc
       f_lo(i)%ulim_proc = min(f_lo(i)%ulim_world, f_lo(i)%llim_proc + f_lo(i)%blocksize - 1)
       f_lo(i)%ulim_alloc = max(f_lo(i)%llim_proc, f_lo(i)%ulim_proc)
    end do

  end subroutine init_fields_layouts

  function im_idx_f (lo, i)
    implicit none
    integer :: im_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    im_idx_f = 1 + mod((i - lo%llim_world)/lo%nindex, lo%M_class)
  end function im_idx_f

  function if_idx_f (lo, i)
    implicit none
    integer :: if_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    if_idx_f = 1 + mod(i - lo%llim_world, lo%nindex)
  end function if_idx_f

  function idx_f (lo, if, im)
    implicit none
    integer :: idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: if, im
    idx_f = if-1 + lo%nindex*(im-1)
  end function idx_f

  function proc_id_f (lo, i)
    implicit none
    integer :: proc_id_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_f = i/lo%blocksize
  end function proc_id_f

  function idx_local_f (lo, if, im)
    implicit none
    logical :: idx_local_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: if, im

    idx_local_f = idx_local(lo, idx(lo, if, im))
  end function idx_local_f

  function ig_local_f (lo, ig)
    implicit none
    logical :: ig_local_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_f = lo%iproc == proc_id(lo, ig)
  end function ig_local_f

  function in_idx_f (lo, i)
    implicit none
    integer :: in_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    in_idx_f = 1 + mod((i - lo%llim_world)/lo%nidx, lo%N_class)
  end function in_idx_f

  function ig_idx_f (lo, i)
    implicit none
    integer :: ig_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_f = -lo%ntgrid + mod((i - lo%llim_world)/lo%nfield, (2*lo%ntgrid+1))
  end function ig_idx_f

  function ij_idx_f (lo, ig, if, n)
    implicit none
    integer :: ij_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, if, n
    ij_idx_f = ig+lo%ntgrid + (2*lo%ntgrid+1)*(if-1 + lo%nfield*(n-1))
  end function ij_idx_f

  function ifield_idx_f (lo, i)
    implicit none
    integer :: ifield_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ifield_idx_f = 1 + mod((i - lo%llim_world), lo%nfield)
  end function ifield_idx_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fast field layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_jfields_layouts (nfield, nindex, naky, ntheta0, i_class)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: nfield, nindex, naky, ntheta0, i_class
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    jf_lo%iproc = iproc
    jf_lo%nindex = nindex
    jf_lo%ntgrid = (nindex/nfield-1)/2
    jf_lo%nfield = nfield
    jf_lo%naky = naky
    jf_lo%ntheta0 = ntheta0
    jf_lo%llim_world = 0
    jf_lo%ulim_world = nindex*ntheta0*naky - 1
    jf_lo%blocksize = jf_lo%ulim_world/nproc + 1
    jf_lo%llim_proc = jf_lo%blocksize*iproc
    jf_lo%ulim_proc = min(jf_lo%ulim_world, jf_lo%llim_proc + jf_lo%blocksize - 1)
    jf_lo%ulim_alloc = max(jf_lo%llim_proc, jf_lo%ulim_proc)

    allocate (ij(jf_lo%llim_proc:jf_lo%ulim_alloc))
    allocate (mj(jf_lo%llim_proc:jf_lo%ulim_alloc))
    allocate (dj(i_class,jf_lo%llim_proc:jf_lo%ulim_alloc))
    ij = 1  ; mj = 1;  dj = 0
    
  end subroutine init_jfields_layouts

  function ik_idx_jf (lo, i)
    implicit none
    integer :: ik_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_jf = 1 + mod((i - lo%llim_world)/lo%nindex/lo%ntheta0, lo%naky)
  end function ik_idx_jf

  function it_idx_jf (lo, i)
    implicit none
    integer :: it_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_jf = 1 + mod((i - lo%llim_world)/lo%nindex, lo%ntheta0)
  end function it_idx_jf

  function if_idx_jf (lo, i)
    implicit none
    integer :: if_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    if_idx_jf = 1 + mod(i - lo%llim_world, lo%nindex)
  end function if_idx_jf

  function idx_jf (lo, if, ik, it)
    implicit none
    integer :: idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: if, ik, it
    idx_jf = if-1 + lo%nindex*(it-1 + lo%ntheta0*(ik-1))
  end function idx_jf

  function ij_idx_jf (lo, ig, if, ik, it)
    implicit none
    integer :: ij_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, if, ik, it
    ij_idx_jf = ig+lo%ntgrid + (2*lo%ntgrid+1)*(if-1+lo%nfield*(it-1+lo%ntheta0*(ik-1)))
  end function ij_idx_jf

  function proc_id_jf (lo, i)
    implicit none
    integer :: proc_id_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_jf = i/lo%blocksize
  end function proc_id_jf

  function idx_local_jf (lo, if, ik, it)
    implicit none
    logical :: idx_local_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: if, ik, it

    idx_local_jf = idx_local(lo, idx(lo, if, ik, it))
  end function idx_local_jf

  function ig_local_jf (lo, ig)
    implicit none
    logical :: ig_local_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_jf = lo%iproc == proc_id(lo, ig)
  end function ig_local_jf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_x_transform_layouts &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx)
    use mp, only: iproc, nproc, barrier
    implicit none
    integer, intent (in) :: ntgrid, nvgrid, naky, ntheta0, nmu, nspec, nx

    if (initialized_x_transform) return
    initialized_x_transform = .true.

    xxf_lo%iproc = iproc
    xxf_lo%ntgrid = ntgrid
    xxf_lo%ntgridtotal = (2*ntgrid+1)
    xxf_lo%nvgrid = nvgrid
    xxf_lo%naky = naky
    xxf_lo%ntheta0 = ntheta0
    if (nx > ntheta0) then
       xxf_lo%nx = nx
    else
       xxf_lo%nx = (3*ntheta0+1)/2
    end if
    xxf_lo%nadd = xxf_lo%nx - ntheta0
    xxf_lo%nmu = nmu
    xxf_lo%nspec = nspec
    xxf_lo%llim_world = 0
    xxf_lo%ulim_world = naky*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec - 1

    xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
    xxf_lo%llim_proc = xxf_lo%blocksize*iproc
    xxf_lo%ulim_proc &
         = min(xxf_lo%ulim_world, xxf_lo%llim_proc + xxf_lo%blocksize - 1)
    xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)

  end subroutine init_x_transform_layouts

  elemental function is_idx_xxf (lo, i)
    implicit none
    integer :: is_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ! TT: the order of the division does not matter, so no need for branching
    is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/(2*lo%nvgrid+1)/lo%nmu, lo%nspec)

  end function is_idx_xxf

  elemental function imu_idx_xxf (lo, i)
    implicit none
    integer :: imu_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    imu_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/(2*lo%nvgrid + 1), lo%nmu)
  end function imu_idx_xxf

  elemental function iv_idx_xxf (lo, i)
    implicit none
    integer :: iv_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    iv_idx_xxf = -lo%nvgrid + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), 2*lo%nvgrid + 1)
  end function iv_idx_xxf

  elemental function ig_idx_xxf (lo, i)
    implicit none
    integer :: ig_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
  end function ig_idx_xxf

  elemental function ik_idx_xxf (lo, i)
    implicit none
    integer :: ik_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
  end function ik_idx_xxf

  elemental function idx_xxf (lo, ig, iv, ik, imu, is)
    implicit none
    integer :: idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, ik, imu, is

    idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(iv+lo%nvgrid &
         + (2*lo%nvgrid+1)*(imu-1 + lo%nmu*(is-1))))
  end function idx_xxf

  elemental function proc_id_xxf (lo, i)
    implicit none
    integer :: proc_id_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
    proc_id_xxf = i/lo%blocksize

  end function proc_id_xxf

  elemental function idx_local_xxf (lo, ig, iv, ik, imu, is)
    implicit none
    logical :: idx_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, ik, imu, is
    idx_local_xxf = idx_local (lo, idx(lo, ig, iv, ik, imu, is))
  end function idx_local_xxf

  elemental function ig_local_xxf (lo, i)
    implicit none
    logical ig_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_xxf = lo%iproc == proc_id(lo, i)
  end function ig_local_xxf

  elemental function xxf_ky_is_zero (lo, i)
    implicit none
    logical xxf_ky_is_zero
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    xxf_ky_is_zero = 0 == mod(i, lo%naky)

  end function xxf_ky_is_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Y-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_y_transform_layouts &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, nvgrid, naky, ntheta0, nmu, nspec
    integer, intent (in) :: nx, ny
    integer :: nnx, nny

    if (initialized_y_transform) return
    initialized_y_transform = .true.

    if (nx > ntheta0) then
       nnx = nx
    else
       nnx = (3*ntheta0+1)/2
    end if
    if (ny > naky) then
       nny = ny
    else
       nny = 3*naky
    end if

    yxf_lo%iproc = iproc
    yxf_lo%ntgrid = ntgrid
    yxf_lo%ntgridtotal = (2*ntgrid+1)
    yxf_lo%nvgrid = nvgrid
    yxf_lo%nvgridtotal = (2*nvgrid+1)
    yxf_lo%naky = naky
    yxf_lo%ny = nny
    yxf_lo%ntheta0 = ntheta0
    yxf_lo%nx = nnx
    yxf_lo%nmu = nmu
    yxf_lo%nspec = nspec
    yxf_lo%llim_world = 0
    yxf_lo%ulim_world = nnx*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec - 1

    yxf_lo%it_ord=1
    yxf_lo%ig_ord=2
    yxf_lo%iv_ord=3
    yxf_lo%imu_ord=4
    yxf_lo%is_ord=5

    yxf_lo%compound_count(1)=1
    yxf_lo%compound_count(2)=yxf_lo%nx
    yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*yxf_lo%ntgridtotal
    yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*yxf_lo%nvgridtotal
    yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nmu

    yxf_lo%ig_comp=yxf_lo%compound_count(yxf_lo%ig_ord)
    yxf_lo%iv_comp=yxf_lo%compound_count(yxf_lo%iv_ord)
    yxf_lo%it_comp=yxf_lo%compound_count(yxf_lo%it_ord)
    yxf_lo%imu_comp=yxf_lo%compound_count(yxf_lo%imu_ord)
    yxf_lo%is_comp=yxf_lo%compound_count(yxf_lo%is_ord)

    yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
    yxf_lo%llim_proc = yxf_lo%blocksize*iproc
    yxf_lo%ulim_proc &
         = min(yxf_lo%ulim_world, yxf_lo%llim_proc + yxf_lo%blocksize - 1)
    yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)
    
  end subroutine init_y_transform_layouts

  elemental function is_idx_yxf (lo, i)
    implicit none
    integer :: is_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_yxf = 1 + mod((i-lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_yxf

  elemental function imu_idx_yxf (lo, i)
    implicit none
    integer :: imu_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    imu_idx_yxf = 1 + mod((i-lo%llim_world)/lo%imu_comp, lo%nmu)
  end function imu_idx_yxf

  elemental function iv_idx_yxf (lo, i)
    implicit none
    integer :: iv_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    iv_idx_yxf = -lo%nvgrid + mod((i-lo%llim_world)/lo%iv_comp, lo%nvgridtotal)
  end function iv_idx_yxf

  elemental function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_yxf = -lo%ntgrid + mod((i-lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_yxf

  elemental function it_idx_yxf (lo, i)
    implicit none
    integer :: it_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_yxf = 1 + mod((i-lo%llim_world)/lo%it_comp, lo%nx)
  end function it_idx_yxf

  elemental function idx_yxf (lo, ig, iv, it, imu, is)
    implicit none
    integer :: idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, it, imu, is

    idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(iv+lo%nvgrid &
         + (2*lo%nvgrid+1)*(imu-1 + lo%nmu*(is-1))))
  end function idx_yxf

  elemental function proc_id_yxf (lo, i)
    implicit none
    integer :: proc_id_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_yxf = i/lo%blocksize

  end function proc_id_yxf

  elemental function idx_local_yxf (lo, ig, iv, it, imu, is)
    implicit none
    logical :: idx_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, iv, it, imu, is
    idx_local_yxf = idx_local (lo, idx(lo, ig, iv, it, imu, is))
  end function idx_local_yxf

  elemental function ig_local_yxf (lo, i)
    implicit none
    logical ig_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_yxf = lo%iproc == proc_id(lo, i)
  end function ig_local_yxf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transformation subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# ifdef USE_C_INDEX
  subroutine gidx2xxfidx (ig, iv, it_in, iglo, g_lo, xxf_lo, it_out, ixxf)
# else
  elemental subroutine gidx2xxfidx (ig, iv, it_in, iglo, g_lo, xxf_lo, it_out, ixxf)
# endif

    implicit none
    integer, intent (in) :: ig, iv, it_in, iglo
    type (g_layout_type), intent (in) :: g_lo
    type (xxf_layout_type), intent (in) :: xxf_lo
    integer, intent (out) :: ixxf, it_out

    if (it_in > (xxf_lo%ntheta0+1)/2) then
       it_out = it_in - xxf_lo%ntheta0 + xxf_lo%nx
    else
       it_out = it_in
    end if

    ixxf = idx(xxf_lo, ig, iv, ik_idx(g_lo,iglo), &
         imu_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2xxfidx

  elemental subroutine xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
    implicit none
    integer, intent (in) :: it, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (yxf_layout_type), intent (in) :: yxf_lo
    integer, intent (out) :: ik, iyxf

    ik = ik_idx(xxf_lo,ixxf)
    iyxf = idx(yxf_lo, ig_idx(xxf_lo,ixxf), iv_idx(xxf_lo,ixxf), &
         it, imu_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
  end subroutine xxfidx2yxfidx

  elemental subroutine yxfidx2xxfidx (ik, iyxf, yxf_lo, xxf_lo, it, ixxf)
    implicit none
    integer, intent (in) :: ik, iyxf
    type (yxf_layout_type), intent (in) :: yxf_lo
    type (xxf_layout_type), intent (in) :: xxf_lo
    integer, intent (out) :: it, ixxf
    integer :: ik0

    ik0 = ik
    if (ik0 > xxf_lo%naky) then
       it = -999999
       ixxf = -999999
       return
    end if

    it = it_idx(yxf_lo,iyxf)
    ixxf = idx(xxf_lo, ig_idx(yxf_lo,iyxf), iv_idx(yxf_lo,iyxf), &
         ik0, imu_idx(yxf_lo,iyxf), is_idx(yxf_lo,iyxf))
  end subroutine yxfidx2xxfidx

  ! subroutine factors (n, j, div)
  !   integer, intent (in) :: n
  !   integer, intent (out) :: j
  !   integer, dimension (:), intent (out) :: div
  !   integer :: i, imax

  !   do i=2,n
  !      if (mod(n,i)==0) exit
  !   end do
  !   imax = n/i
  !   j=1
  !   do i=1,imax
  !      if (mod(n,i)==0) then
  !         div(j) = i
  !         j=j+1
  !      end if
  !   end do
  !   div(j) = n
  ! end subroutine factors

  subroutine finish_layouts

    implicit none

    deallocate (f_lo, ij, mj, dj)

  end subroutine finish_layouts

end module gs2_layouts
