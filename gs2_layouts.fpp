! Modifications for unbalanced decomposition functionality:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2012
! on behalf of EPSRC for the HECToR project
# include "define.inc"

module gs2_layouts

  use layouts_type, only: g_layout_type

  implicit none
  private

  public :: layout, finish_layouts, local_field_solve
  public :: is_kx_local
  
  public :: init_dist_fn_layouts, init_gs2_layouts
  public :: wnml_gs2_layouts
  public :: g_lo

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_jfields_layouts
  public :: jf_lo, jf_layout_type
  public :: mj, ij, dj

  public :: init_x_transform_layouts, init_y_transform_layouts
  public :: calculate_unbalanced_x, calculate_unbalanced_y, calculate_idle_processes
  public :: xxf_lo, xxf_layout_type, yxf_lo, yxf_layout_type
  public :: gidx2xxfidx, xxfidx2yxfidx, yxfidx2xxfidx, xxfidx2gidx
  public :: xxf_ky_is_zero

  public :: ig_idx, ik_idx, it_idx, iv_idx, imu_idx, is_idx, if_idx
  public :: im_idx, in_idx, ij_idx, ifield_idx
  public :: idx, proc_id, idx_local

  logical :: initialized_x_transform = .false.
  logical :: initialized_y_transform = .false.

  logical :: local_field_solve, unbalanced_xxf, unbalanced_yxf
  real :: max_unbalanced_xxf, max_unbalanced_yxf
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
     integer :: ntgrid, nvgrid, naky, ny, ntheta0, nx, nmu, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
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
     module procedure it_idx_g
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
  subroutine wnml_gs2_layouts(unit)
    implicit none
    integer :: unit
    if (.not. exist) return
       write (unit, *)
       write (unit, fmt="(' &',a)") "layouts_knobs"
       write (unit, fmt="(' layout = ',a)") '"'//trim(layout)//'"'
       write (unit, fmt="(' local_field_solve = ',L1)") local_field_solve
       write (unit, fmt="(' /')")
  end subroutine wnml_gs2_layouts

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
    namelist /layouts_knobs/ layout, local_field_solve, unbalanced_xxf, &
         max_unbalanced_xxf, unbalanced_yxf, max_unbalanced_yxf

    local_field_solve = .false.
    unbalanced_xxf = .false.
    unbalanced_yxf = .false.
    max_unbalanced_xxf = 0.0
    max_unbalanced_yxf = 0.0
    layout = 'xyms'
    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)
    if (layout.ne.'xyms' .and. layout.ne.'xmys' .and. layout.ne.'yxms' &
         .and. layout.ne.'mxys' .and. layout.ne.'ymxs' .and. layout.ne.'myxs') &
    then
       write(6,*) "gs2_layouts: read_parameters finds illegal layout=",layout," =>stop"
       stop
    endif

! max_unbalanced_xxf and max_unbalanced_yxf have a maximum range of 0.0 - 1.0
! so check here whether the user has specified the range correctly and adjust if 
! they have not.
    if(max_unbalanced_xxf .gt. 1.0) then
       max_unbalanced_xxf = 1.0
       write(*,*) 'max_unbalanced_xxf too large, setting to 1.0'
    else if(max_unbalanced_xxf .lt. 0.0) then
       max_unbalanced_xxf = 0.0
       write(*,*) 'max_unbalanced_xxf too small, setting to 0.0'
    end if
    if(max_unbalanced_yxf .gt. 1.0) then
       max_unbalanced_yxf = 1.0
       write(*,*) 'max_unbalanced_yxf too large, setting to 1.0'
    else if(max_unbalanced_yxf .lt. 0.0) then
       max_unbalanced_yxf = 0.0
       write(*,*) 'max_unbalanced_yxf too small, setting to 0.0'
    end if



  end subroutine read_parameters
    
  subroutine broadcast_results
    use mp, only: broadcast
    implicit none

    call broadcast (layout)
    call broadcast (local_field_solve)
    call broadcast (unbalanced_xxf)
    call broadcast (unbalanced_yxf)
    call broadcast (max_unbalanced_xxf)
    call broadcast (max_unbalanced_yxf)

  end subroutine broadcast_results

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dist_fn_layouts &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec)

    use mp, only: iproc, nproc, proc0
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
    g_lo%ulim_world = naky*ntheta0*nmu*nspec - 1
      
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

  subroutine is_kx_local(nspec, nmu, naky, ntheta0, kx_local)  

    use mp, only: nproc
    implicit none
    integer, intent (in) :: nspec, nmu, naky, ntheta0
    logical, intent (out) :: kx_local
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nmfacs, nyfacs, nxfacs
    integer :: i

    kx_local = .true.

    allocate (facs(max(nspec,nmu,naky)/2+1,4))

    select case (layout)

    case ('xyms')
       
       call factors (nspec, nsfacs, facs(:,1))
       call factors (nmu, nmfacs, facs(:,2))
       call factors (naky, nyfacs, facs(:,3))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nmfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nyfacs-1
          if (nproc == facs(i,3)*naky*nspec) goto 100
       end do

    case ('xmys')
       
       call factors (nspec, nsfacs, facs(:,1))
       call factors (naky, nyfacs, facs(:,2))
       call factors (nmu, nmfacs, facs(:,3))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nmfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nyfacs-1
          if (nproc == facs(i,3)*nmu*nspec) goto 100
       end do

    case ('mxys')
       
       call factors (nspec, nsfacs, facs(:,1))
       call factors (naky, nyfacs, facs(:,2))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nyfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

    case ('yxms')
       
       call factors (nspec, nsfacs, facs(:,1))
       call factors (nmu, nmfacs, facs(:,2))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nmfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

    case ('myxs')
       
       call factors (nspec, nsfacs, facs(:,1))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do

    case ('ymxs')
       
       call factors (nspec, nsfacs, facs(:,1))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do

    end select

    kx_local = .false.    
                       
100 deallocate (facs)

  end subroutine is_kx_local

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
    ! TT: the order of the division doesn't matter, so no need for branching
    is_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%nmu, lo%nspec)
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
       imu_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nmu)
    case ('xyms')
       imu_idx_g = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%nmu)
    case ('mxys')
       imu_idx_g = 1 + mod((i - lo%llim_world), lo%nmu)
    case ('myxs')
       imu_idx_g = 1 + mod((i - lo%llim_world), lo%nmu)
    end select
# endif

  end function imu_idx_g

# ifdef USE_C_INDEX
  function it_idx_g (lo, i)
# else
  elemental function it_idx_g (lo, i)
# endif

    implicit none

    integer :: it_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
# ifdef USE_C_INDEX
    interface
       function it_idx_g_c (lo,num)
         use layouts_type, only: g_layout_type
         integer :: it_idx_g_c
         type (g_layout_type) :: lo
         integer :: num
       end function it_idx_g_c
    end interface
    it_idx_g = it_idx_g_c (lo,i)
# else
    select case (layout)
    case ('xyms')
       it_idx_g = 1 + mod(i - lo%llim_world, lo%ntheta0)
    case ('xmys')
       it_idx_g = 1 + mod(i - lo%llim_world, lo%ntheta0)
    case ('yxms')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('mxys')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%nmu, lo%ntheta0)
    case ('myxs')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%nmu/lo%naky, lo%ntheta0)
    case ('ymxs')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%nmu, lo%ntheta0)
    end select
# endif

  end function it_idx_g

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
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
    case ('xmys')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%nmu, lo%naky)
    case ('yxms')
       ik_idx_g = 1 + mod(i - lo%llim_world, lo%naky)
    case ('mxys')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nmu/lo%ntheta0, lo%naky)
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
  function idx_g (lo, ik, it, imu, is)
# else
  elemental function idx_g (lo, ik, it, imu, is)
# endif

    implicit none

    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, imu, is

# ifdef USE_C_INDEX
    interface
       function idx_g_c (lo,ik,it,imu,is)
         use layouts_type, only: g_layout_type
         integer :: idx_g_c
         type (g_layout_type) :: lo
         integer :: ik,it,imu,is
       end function idx_g_c
    end interface
    idx_g = idx_g_c (lo, ik, it, imu, is)
# else
    select case (layout)
    case ('xyms')
       idx_g = it-1 + lo%ntheta0*(ik-1 + lo%naky*(imu-1 + lo%nmu*(is-1)))
    case ('xmys')
       idx_g = it-1 + lo%ntheta0*(imu-1 + lo%nmu*(ik-1 + lo%naky*(is-1)))
    case ('yxms')
       idx_g = ik-1 + lo%naky*(it-1 + lo%ntheta0*(imu-1 + lo%nmu*(is-1)))
    case ('mxys')
       idx_g = imu-1 + lo%nmu*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
    case ('myxs')
       idx_g = imu-1 + lo%nmu*(ik-1 + lo%naky*(it-1 + lo%ntheta0*(is-1)))
    case ('ymxs')
       idx_g = ik-1 + lo%naky*(imu-1 + lo%nmu*(it-1 + lo%ntheta0*(is-1)))
    end select
# endif

  end function idx_g

! TT>
# ifdef USE_C_INDEX
  function idx_local_g (lo, ik, it, imu, is)
# else
  elemental function idx_local_g (lo, ik, it, imu, is)
# endif
! <TT
    implicit none
    logical :: idx_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, imu, is

    idx_local_g = idx_local(lo, idx(lo, ik, it, imu, is))
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
    use mp, only: iproc, nproc, proc0, barrier
    implicit none
    integer, intent (in) :: ntgrid, nvgrid, naky, ntheta0, nmu, nspec, nx
    integer :: nprocset, ngroup, ip, nblock
    real :: unbalanced_amount

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

    ! AJ November 2011
    ! unbalanced_xxf is a variable initialised in the input file
    ! which if set to .true. will enable the code below to
    ! investigate whether an unbalanced decomposition can be
    ! constructed.  Whether the unbalanced decomposition is used
    ! or not is also dependent on the variable max_unbalanced_xxf
    ! which is also set in the input file and defines the maximum
    ! amount of imbalance permitted.
    ! The code that constructs the unbalance decomposition works
    ! through the different indicies that are used in the xxf_lo
    ! to construct the full xxf_lo data space, namely:
    ! naky*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec
    ! Note: We precalculate 2*ntgrid+1 as ntgridtotal.
    if (unbalanced_xxf) then
       
       call calculate_unbalanced_x(nproc, iproc, unbalanced_amount)		
       
    end if
    
    ! If we are not using the unbalanced code, either because the
    ! user has chosen not to do this in the code or because the
    ! calculated imbalance was too large then use the original
    ! decomposition.
    if (.not. unbalanced_xxf) then
       
       xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
       xxf_lo%llim_proc = xxf_lo%blocksize*iproc
       xxf_lo%ulim_proc &
            = min(xxf_lo%ulim_world, xxf_lo%llim_proc + xxf_lo%blocksize - 1)
       xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)
       
    else
       
       if(proc0) then	
          write(*, fmt="('Using unbalanced decomposition for xxf. '&
               'Unbalanced fraction',F6.2)") unbalanced_amount
       end if
       
    end if

  end subroutine init_x_transform_layouts

  subroutine calculate_unbalanced_x(nproc, iproc, unbalanced_amount)
  !====================================================================
  ! AJ, November 2011: New code from DCSE project
  ! This subroutine (calculate_unbalanced_x) wraps up all the 
  ! functionality required to calculate the xxf blocksizes and 
  ! work out if and unbalanced decomposition should be used, and if so 
  ! the size of the unbalanced blocks.
  !
  ! This routine takes the number of processes to be used and the rank 
  ! of the calling process and returns the amount of unbalance in the 
  ! calculate decomposition (in a real variable with a range between 
  ! 0.0 and 1.0 where 1.0 is a 100% difference between the small and 
  ! large block size and 0.0 is no difference between the small and 
  ! large block sizes).
  !
  ! This routine can be called from routines within the gs2_layouts
  ! module (specifically init_x_transform_layouts) or from outside 
  ! this module (specifically it is currently called from ingen to 
  ! provide suggestions for unbalanced decompositions and process 
  ! counts).  However, this routine relies on init_x_transform_layouts 
  ! having been called to initialise some of the values it used 
  ! (those in xxf_lo%) which is why there is a check whether the 
  ! initialize_x_transform variable has been set to true.  If not this 
  ! routine will not run and will notify the user of the problem.
  !====================================================================

    implicit none
    integer, intent (in) :: nproc, iproc
    real, intent (out) :: unbalanced_amount
    integer :: level_proc_num, i, j, k, m

    ! Ensure that the xxf_lo% data has been properly initialized as this 
    ! routine relies on some data from that data structure.  If it has not 
    ! then abort this routine.
    if(.not. initialized_x_transform) then
       write(*,*) 'X Transform data structures not initialized so calculate_unbalanced_x subroutine will not operate correctly'
       write(*,*) 'Aborting subroutine calculate_unbalanced_x'
       return
    end if
    
    ! The xxf_lo%blocksize is initialised to the value that
    ! is used if this unbalanced code is not used (the
    ! original code).
    xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
    xxf_lo%small_block_balance_factor = 1
    xxf_lo%large_block_balance_factor = 1                                                        
    
    ! Initialise the number of processors/cores used in the
    ! decomposition as the total number available.
    level_proc_num = nproc
    
    ! The first factor to be considered is nspec (number
    ! of species).
    k = xxf_lo%nspec
    
    ! Calculate the factors of nspec and the number of
    ! processors/cores we have left from the decomposition.
    ! Further details of i, m, and j are provided in the
    ! comments for the subroutine calculate_block_breakdown
    ! but briefly described i is the remainder left from k/level_proc_num,
    ! m is the remainder of level_proc_num/k and j is the integer
    ! result of k/level_proc_num.
    call calculate_block_breakdown(k, i, m, j, level_proc_num)
    
    ! If j = 0 this means that level_proc_num is larger than k
    ! the factor we are looking to divide so we have reached the
    ! point where we can't divide anymore and can now build the
    ! decomposition.
    if(j .eq. 0) then
       
       ! If we have got to this point then k is larger than the
       ! number of processes we have so we can try and split it
       ! up over the processes.
       ! If m = 0 then we know that k exactly divides into the
       ! number of processes we have at this point.  This means
       ! we can eliminate this factor from the decomposition
       ! and reduce the number of processes we have available
       ! by that factor.   Therefore we divide the number
       ! of processes we have by k and set the factor we need
       ! to decompose to 1.
       ! If m does not equal zero then we keep the current number
       ! of processes and pass the factor (k) on to the next
       ! level of the decomposition.
       if(m .eq. 0) then
          level_proc_num = level_proc_num/k
          k = 1
       end if
       
       ! Multiple the previous factor by the next factor to be
       ! be used.  If the previous factor divided exactly by the
       ! processes available then k will be 1 so we will only be
       ! considering this factor.
       k = xxf_lo%nmu * k
       
       ! The next code follows the pattern described above
       ! moving through all the other factor, i.e.:
       ! nmu*nvgridtotal*ntgridtotal*naky
       ! until it gets to a point where j = 0 (i.e. level_proc_num
       ! is greater than the factors being decomposed), or we
       ! have run out of factors to be decomposed.
       ! Once either of those points is reached then the
       ! code moves on to calculate the blocksizes in the
       ! else clauses of this code below.
       ! This part of the code is commented in the first occurence
       ! of "if(i .ne. 0) then" below.
       call calculate_block_breakdown(k, i, m, j, level_proc_num)
       
       if(j .eq. 0) then
          
          if(m .eq. 0) then
             level_proc_num = level_proc_num/k
             k = 1
          end if
          
          k = xxf_lo%nmu * k
          
          call calculate_block_breakdown(k, i, m, j, level_proc_num)
          
          if(j .eq. 0) then
             
             if(m .eq. 0) then
                level_proc_num = level_proc_num/k
                k = 1
             end if
             
             k = (2*xxf_lo%nvgrid+1) * k
             call calculate_block_breakdown(k, i, m, j, level_proc_num)
             
             if(j .eq. 0) then
                
                if(m .eq. 0) then
                   level_proc_num = level_proc_num/k
                   k = 1
                end if
                
                k = xxf_lo%ntgridtotal * k
                call calculate_block_breakdown(k, i, m, j, level_proc_num)
                
                if(j .eq. 0) then
                   
                   if(m .eq. 0) then
                      level_proc_num = level_proc_num/k
                      k = 1
                   end if
                   
                   k = xxf_lo%naky * k
                   call calculate_block_breakdown(k, i, m, j, level_proc_num)
                   
                   ! At this point we have run out of factors to
                   ! decompose.  Now we check whether i is not
                   ! equal to 0.  If i is equal to 0 (remember that
                   ! i is the remainder left when k/level_proc_num)
                   ! then this is not an unbalanced decomposition so
                   ! we do not have to calculate the unbalanced decompostion.
                   ! If i is not equal to 0 then k cannot be exactly divided
                   ! by level_proc_num.  This means we cannot evenly divide
                   ! the decomposition over the processes we have available so
                   ! we need to find a way to create a different decomposition.
                   if(i .ne. 0) then
                      
                      ! Calculate the least unbalanced split of k over level_proc_num and also
                      ! how what factor of processes are assigned to each block size.
                      call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
                      ! Calculate the actual block sizes using the factors calculated above and the remaining data space to be
                      ! decomposed.  In this instance we are at the lowest level of the decomposition so there is nothing left to decompose
                      ! so the remaining data space is 1.  For other levels of the decomposition this 1 is replaced by the part
                      ! of the data space that has not been split up yet.
                      call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, nproc, &
                           1, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, &
                           xxf_lo%ulim_alloc)
                      
                   end if
                   
                else
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)     
                      ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                      ! naky for this level of the decomposition.
                      call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, nproc, &
                           xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple,  xxf_lo%llim_proc, xxf_lo%ulim_proc, &
                           xxf_lo%ulim_alloc)
                      
                   end if
                   
                end if
                
             else
                
                if(i .ne. 0) then
                   
                   call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)                         
                   ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                   ! naky,ntgridtotal for this level of the decomposition.
                   call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, nproc, &
                        xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, &
                        xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
                   
                end if
                
             end if
             
             
          else
             
             if(i .ne. 0) then
                
                call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nvgridtotal for this level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, nproc, &
                     (2*xxf_lo%nvgrid+1)*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, &
                     xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
                
             end if
             
          end if
          
       else
          
          if(i .ne. 0) then
             
             call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
             
             ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
             ! naky,ntgridtotal,2*nvgrid+1,nmu for this layout and level of the decomposition.
             call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, nproc, &
                  xxf_lo%nmu*(2*xxf_lo%nvgrid+1)*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, &
                  xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
             
          end if
          
       end if
       
    else
       
       if(i .ne. 0) then
          
          call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
          ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
          ! naky,ntgridtotal,2*nvgrid+1,nmu for this level of the decomposition.
          call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, nproc, &
               xxf_lo%nmu*(2*xxf_lo%nvgrid+1)*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, &
               xxf_lo%large_block_size, xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
          
       end if
       
    end if
    
    
    ! If the decomposition has calculated that the large block and small blocks are equal then
    ! we don't have an unbalanced decomposition so set the unbalanced_amount to 0.
    ! large_block_balance_factor and small_block_balance_factor are initialised to 1 at
    ! the start of this decomposition code so if they haven't been changed we already
    ! have a balanced decomposition.
    if(xxf_lo%large_block_balance_factor .eq. 1 .and. xxf_lo%small_block_balance_factor .eq. 1) then
       
       unbalanced_amount = 0
       
    else 
       
       ! If there is an unbalanced decomposition work out the percentage of
       ! difference between the two blocks.  This is used to ensure that
       ! we don't create decompositions that have significant differences
       ! between the two block sizes which would impact the amount of
       ! computation the different groups of processes have to perform.
       unbalanced_amount = real(xxf_lo%large_block_size)/real(xxf_lo%small_block_size)
       unbalanced_amount = unbalanced_amount - 1
       
    end if
    
    ! If we calculate that there is not an unbalanced decomposition or that the
    ! amount of unbalance is larger than a integer the user sets in the input file
    ! called: max_unbalanced_xxf ; then do not use the unbalanced decomposition.
    if (unbalanced_amount .gt. max_unbalanced_xxf .or. unbalanced_amount .eq. 0) then
       
       unbalanced_xxf = .false.       
       
    end if
    
    
  end subroutine calculate_unbalanced_x
  
  subroutine calculate_block_breakdown (k, i, m, j, level_proc_num)
  !====================================================================
  ! AJ, November 2011: New code from DCSE project
  ! This subroutine (calculate_block_breakdown) is used in the code
  ! to create unbalanced xxf and yxf layouts to address MPI
  ! communication overheads in the redistribution of data for the
  ! FFTs used in the nonlinear calculations to move from Fourier
  ! space to real space and back again.
  !
  ! The routine calculates the factor k divided by a provided number
  ! of processors/cores and returns it in integer j.  As j is an
  ! integer the number returned will be more than zero if k is larger
  ! than level_proc_num and zero if k is smaller than level_proc_num.
  !
  ! The integer i that is returned is the remained when k is divided
  ! by level_proc_num.  This is used in the unbalanced decomposition
  ! code to work out whether the factor divides exactly by the number
  ! of procs  (in level_proc_num).  If it does then it is not an
  ! unbalanced decomposition (so the original decomposition can be
  ! used.
  !
  ! The integer m which is returned by this subroutine is the remainder
  ! of level_proc_num divided by k.  This is used in the unbalanced
  ! decomposition code to work out if the factor k exactly divides
  ! the number of processors/cores in level_proc_num.  If it does
  ! (i.e. m is zero) then the decomposition code can perform this
  ! division and move on to the next factor.  If it does not then
  ! the code keeps the factor k and multiplies it by the next factor
  ! and checks that one (and so on).
  !====================================================================
    implicit none
   
    integer, intent(in) :: k, level_proc_num
    integer, intent(out) :: i, m, j

    i = mod(k, level_proc_num)
    m = mod(level_proc_num, k)
    j = k/level_proc_num

  end subroutine calculate_block_breakdown


  subroutine calculate_unbalanced_decomposition(k, j, i, numsmall, numlarge, localprocs)
  !====================================================================
  ! AJ, November 2011: New code from DCSE project
  ! This subroutine (calculate_unbalanced_decomposition) is used in the code
  ! to create unbalanced xxf and yxf layouts to address MPI
  ! communication overheads in the redistribution of data for the
  ! FFTs used in the nonlinear calculations to move from Fourier
  ! space to real space and back again.
  !
  ! k is an integer input which is the factor that needs to be decomposed
  ! across the number of processes (the input integer localprocs)
  ! that we have at this point in the decomposition.
  !
  ! When this subroutine is called k is larger than
  ! localprocs (this is enforced through the if and else statements
  ! in the code that calls this routine).  It is also only called when
  ! localprocs does not completely/equally divide k (so k/localprocs
  ! will not given a round number).
  !
  ! The subroutine returns four integers; i, j, numlarge, and numsmall
  !
  ! j is calculated as the integer division of k/localprocs.  As this
  ! is an integer division the result is rounded down (remembering
  ! that k/localprocs will not be an equal division).  j is used by the
  ! small block factor (the relative size of the small block).
  !
  ! i is the large block factor.  As j is the rounded down division,
  ! i = j + 1.  This gives two different block sizes (when these are
  ! later calculated).
  !
  ! numlarge is used to work out how many processes will use the large
  ! block size in the decomposition.
  !
  ! numsmall is used to work out how many processes will use the small
  ! block size in the decomposition.
  !
  !====================================================================
    implicit none

    integer, intent(in) :: k, localprocs
    integer, intent(out) :: i, j, numlarge, numsmall

    ! To calculate the number of processes to give the large block size
    ! we work out the remainder left from dividing the input factor
    ! but the processes we have left.  This remainder is proportionate
    ! to the difference in the two block sizes.  For example, if k = 62
    ! and localprocs = 3 the code will work out i to be 21 and j to be
    ! 20, although k/localprocs is actually equal to 20 and 2/3.
    ! The code will also work out numlarge to be 2 and numsmall to
    ! be 1, so of the total processes in the simulatino 1/3 will have
    ! the small block and 2/3 the large block.  This maps to the 2/3
    ! in the actual division as opposed to the integer division (i.e
    ! the factor we have removed and replaced with small and large
    ! sizes, so replacing 20 2/3 by 1/3 x 20 and 2/3 * 21.
    numlarge = mod(k,localprocs)
    j  = k/localprocs
    i = j + 1  
    numsmall = localprocs - numlarge

  end subroutine calculate_unbalanced_decomposition


  subroutine calculate_block_size(iproc, numsmall, numlarge, smalldecomp, largedecomp, nproc, sizeblock, blocksize, smallblocksize, largeblocksize, block_multiple, llim, ulim, ulim_alloc)
  !====================================================================
  ! AJ, November 2011: New code from DCSE project
  ! This subroutine (calculate_block_size) is used in the code
  ! to create unbalanced xxf and yxf layouts to address MPI
  ! communication overheads in the redistribution of data for the
  ! FFTs used in the nonlinear calculations to move from Fourier
  ! space to real space and back again.
  !
  ! Input iproc is this processes MPI id.
  !
  ! Input numsmall is the number of processes that will use the small
  ! block size (calculated in subroutine
  ! calculate_unbalanced_decomposition).
  !
  ! Input numlarge is the number of processes that will use the large
  ! block size (calculated in suboutine
  ! calcuate_unbalanced_decomposition).
  !
  ! Input smalldecomp is the small_block_balance_factor of the
  ! decomposition (calculated in subroutine
  ! calculate_unbalanced_decomposition where it is referred to as j.
  !
  ! Input largedecomp is the large_block_balance_factor of the
  ! decomposition (calculated in subroutine
  ! calculate_unbalanced_decomposition where it is referred to as i.
  !
  ! Input nproc is the total number of MPI processes that this program
  ! is using.
  !
  ! Input sizeblock is the factor that needs to be split over the
  ! cores/processors we have (i.e. everything that has not yet been
  ! decomposed).
  !
  ! Output blocksize is the blocksize assigned to this actual process
  ! (to be used throughout the code of the decompositions).
  !
  ! Output smallblocksize is the small block size.
  !
  ! Output largeblocksize is the large block size.
  !
  ! Output block_multiple is the size of the small block and large block
  ! added together.  This is used to calculate which process owns a particular
  ! data point at later points in the code (particular in the subroutine 
  ! proc_id).
  ! 
  ! Output llim is the lower limit of the data block for this process.
  !
  ! Output ulim is the upper limit of the data block for this process.
  !
  ! Output ulim_alloc is the upper allocation limit for this processes block.
  !
  !====================================================================
    implicit none

    integer, intent(in) :: iproc, numsmall, numlarge, smalldecomp, largedecomp, nproc, sizeblock
    integer, intent(out) :: blocksize, smallblocksize, largeblocksize, block_multiple, llim, ulim, ulim_alloc
    integer :: modproc, procfactors

    ! Small blocksize is simply calculated by taking the factors left to be
    ! distributed and multiplying it by the small_block_balance_factor.
    smallblocksize = sizeblock * smalldecomp
    ! Likewise large blocksize is the factors multiplied by the
    ! large_block_balance_factor.
    largeblocksize = sizeblock * largedecomp

    ! The block multiple is the chunks that decomposition blocks are arranged
    ! in.  For instance if we have a decomposition with 3 blocks, one small 
    ! of size 640 elements and two large of 672 elements then block_multiple
    ! will be equal to 1*640+2*672.  This is then used in the proc_id code 
    ! to isolate which block an element id belongs to by factorising the id 
    ! with this multiple to isolate the id to a small set of points within 
    ! a single chunk (i.e. set of small and large blocks).  So, for instance, 
    ! if the id 15646 is presented to proc_id then you can work out that it is 
    ! in a large block 
    ! (by doing 15646 - ((15646/block_multiple)*block_multiple) in integer 
    ! arithmetic) which would give 1758, which is in the 3 block of that chunk 
    ! of blocks.  With this information it is possible to work out the
    ! particular proc that owns that chunk.  For more information the
    ! subroutines proc_id_xxf and proc_id_yxf
    
    block_multiple = (smallblocksize * numsmall) + (largeblocksize * numlarge)

    ! This is also used in the proc_id functionality as well as working out 
    ! whether this process has a small or large block.
    procfactors = numsmall + numlarge

    ! This is used to calculate whether this process has a small or large block.
    modproc = mod(iproc, procfactors)

    ! Set this processes blocksize depending if it is to use the smallblock
    ! or large block.
    if(modproc .lt. numsmall) then
       blocksize = smallblocksize
    else
       blocksize = largeblocksize
    end if

    ! Calculate the lower limit to the block this process owns
    llim = ((iproc / procfactors) * block_multiple)
    if(modproc .ne. 0) then
       if(modproc .lt. numsmall) then
!AJ	  llim = llim + smallblocksize * (modproc - 1)
          llim = llim + smallblocksize * (modproc)
       else
          llim = llim + (smallblocksize * numsmall) + (largeblocksize  * (modproc - numsmall))
       end if
    end if

    ! The upper block limit is the lower limit plus the blocksize
    ulim = llim + blocksize  - 1
    ! The allocation upper limit is the upper block limit unless this process 
    ! has no elements of this index (in which situation the ulim will be less
    ! than or equal to the llim.  This ensures that ulim_alloc = llim for zero 
    ! sized blocks.
    ulim_alloc = max(llim, ulim)

  end subroutine calculate_block_size

  subroutine calculate_idle_processes(nprocs, idle_percentage)
  !====================================================================
  ! AJ, November 2011: New code from DCSE project
  ! This subroutine (calculate_idle_processes) is used to calculate the 
  ! difference between the number of processes used in the xxf and yxf
  ! data layouts.  This is important as it can affect the amount of 
  ! communication that the code has to undertake when moving between 
  ! linear and non-linear calculations.  
  !
  ! This routine is used by ingen when it is suggesting optimal process
  ! counts for users to flag up when suggested process counts will 
  ! results in there being a significant difference in the processes 
  ! used in the two layouts, and therefore a significant communication 
  ! overhead moving between the two layouts.
  !====================================================================

    implicit none

    integer, intent(in) :: nprocs
    real, intent(out) :: idle_percentage
    integer :: xxf_blocksize, yxf_blocksize
    real :: xxf_usedprocs, xxf_idleprocs
    real :: yxf_usedprocs, yxf_idleprocs
    real :: delta_idle_procs

    ! Ensure that the xxf_lo% and yxf_lo%data has been properly initialized as this 
    ! routine relies on some data from those data structures.  If it has not 
    ! then abort this routine.
    if(.not. initialized_x_transform .and. .not. initialized_y_transform) then
       write(*,*) 'X and/or Y transform data structures not initialized so calculate_idle_processes will not operate correctly'
       write(*,*) 'Aborting subroutine calculate_idle_processes'
       return
    end if

    if(nprocs .lt. 1) then
       write(*,*) 'nprocs value in calculate_idle_processes subroutine is less than 1 which is incorrect.'
       write(*,*) 'calculate_idle_processes aborting.'
       return
    end if

    ! Calculate the standard xxf_blocksize
    xxf_blocksize = xxf_lo%ulim_world/nprocs + 1
    ! Use the blocksize calculated above to calculate how many processes the
    ! xxf space maps to using this block size
    xxf_usedprocs = (xxf_lo%ulim_world+1)/real(xxf_blocksize)  
    ! Now work out how many processes do not have any xxf data space assigned
    ! to them.  This is calculated using real arthimetic so it will also 
    ! include partial data spaces (so for instance it will calculate where 
    ! a process only has half a block assigned to it).
    xxf_idleprocs = nprocs - xxf_usedprocs
 
    ! Calculate the standard yxf_blocksize
    yxf_blocksize = yxf_lo%ulim_world/nprocs + 1
    ! Use the blocksize calculated above to calculate how many processes the
    ! yxf space maps to using this block size    
    yxf_usedprocs = (yxf_lo%ulim_world+1)/real(yxf_blocksize)
    ! Now work out how many processes do not have any yxf data space assigned
    ! to them.  This is calculated using real arthimetic so it will also 
    ! include partial data spaces (so for instance it will calculate where 
    ! a process only has half a block assigned to it).
    yxf_idleprocs = nprocs - yxf_usedprocs

    ! Calculate the difference between the idle processes in the yxf and xxf 
    ! decompositions.  A high delta_idle_procs will cause high communication
    ! costs in the transform routines.
    delta_idle_procs = abs(yxf_idleprocs - xxf_idleprocs)
 
    ! Roughly calculate the percentage of data to be transferred in the
    ! transform between the xxf and yxf data spaces using the delta_idle_procs
    ! variable calculated above.
    if ( delta_idle_procs .le. 1 ) then
       idle_percentage = 0.5d0 * delta_idle_procs
    else
       idle_percentage = (1.0d0 - 1.0d0/(2.0d0 * delta_idle_procs))
    end if

  end subroutine calculate_idle_processes


  elemental function is_idx_xxf (lo, i)
    implicit none
    integer :: is_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ! TT: the order of the division doesn't matter, so no need for branching
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
    integer :: block_offset, offset_block_number, j, k, tempi
    
    !AJ This code has been added to deal with the unbalanced decomposition functionality.
    !AJ If an unbalanced xxf decomposition is being used then the proc_id cannot
    !AJ use a simple lo%blocksize as there will be two separate block
    !AJ sizes used so we have to work out which block this i lives in and 
    !AJ therefore which process it belongs to.
    if (unbalanced_xxf) then
       !AJ block_offset works out how many groups of blocks this i point is
       !AJ after (the small and large blocks used in the decomposition are
       !AJ grouped together).
       block_offset = (i / lo%block_multiple)          
       !AJ j represents how many blocks are in each group of blocks
       j = lo%num_small + lo%num_large
       !AJ offset_block_number is the number of blocks up to the start of the 
       !AJ group of blocks we are considering.
       offset_block_number = block_offset * j
       !AJ tempi represents where this index is inside the group of blocks
       !AJ this i point sits.
       tempi = i - (block_offset * lo%block_multiple)
       !AJ Work through each block in the group of blocks and see if this i 
       !AJ is within that block.  If it is set the proc_id_xxf as this block 
       !AJ owner.          
       do k=1,j
          if(k .le. lo%num_small) then
             !AJ TODO: The if-else construct used below could potentially be rationalised 
             !AJ TODO: to a more efficient formula where:
             !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%small_block_size
             !AJ TODO: Although a method for selecting if tempi is in the small or large 
             !AJ TODO: blocks would have to be provided.
             if(tempi .lt. lo%small_block_size) then
                !AJ (k -1) is the number of blocks that we have already considered 
                !AJ within this group of blocks we have selected. 
                proc_id_xxf =  offset_block_number + (k - 1)
                exit 
             else
                !AJ If the index is not in this block then reduce tempi by a small
                !AJ block size and move on to the next block.
                tempi = tempi - lo%small_block_size
             end if
          else
             !AJ TODO: The if-else construct used below could potentially be rationalised 
             !AJ TODO: to a more efficient formula where:
             !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%large_block_size
             !AJ TODO: Although a method for selecting if tempi is in the small or large 
             !AJ TODO: blocks would have to be provided.
             if(tempi .lt. lo%large_block_size) then
                proc_id_xxf = offset_block_number + (k - 1)
                exit 
             else
                !AJ If the index is not in this block then reduce tempi by a large
                !AJ block size and move on to the next block.
                tempi = tempi - lo%large_block_size
             end if
          end if
       end do
    else
       !AJ This code is called if the unbalanced decomposition is not being used.
       proc_id_xxf = i/lo%blocksize
    end if

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
    use mp, only: iproc, nproc, proc0
    implicit none
    integer, intent (in) :: ntgrid, nvgrid, naky, ntheta0, nmu, nspec
    integer, intent (in) :: nx, ny
    integer :: nnx, nny, ngroup, nprocset, nblock
    real :: unbalanced_amount

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
    yxf_lo%naky = naky
    yxf_lo%ny = nny
    yxf_lo%ntheta0 = ntheta0
    yxf_lo%nx = nnx
    yxf_lo%nmu = nmu
    yxf_lo%nspec = nspec
    yxf_lo%llim_world = 0
    yxf_lo%ulim_world = nnx*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec - 1

    ! AJ November 2011
    ! unbalanced_yxf is a variable initialised in the input file
    ! which if set to .true. will enable the code below to
    ! investigate whether an unbalanced decomposition can be
    ! constructed.  Whether the unbalanced decomposition is used
    ! or not is also dependent on the variable max_unbalanced_yxf
    ! which is also set in the input file and defines the maximum
    ! amount of imbalance permitted.
    ! The code that constructs the unbalance decomposition works
    ! through the different indicies that are used in the yxf_lo
    ! to construct the full xxf_lo data space, namely:
    ! nxx*(2*ntgrid+1)*(2*nvgrid+1)*nmu*nspec
    ! Note: We precalculate 2*ntgrid+1 as ntgridtotal.
    ! This functionality is the same as the functionality in
    ! init_x_transform_layouts which is extensively commented.
    ! Please see init_x_transform_layouts for more details on the 
    ! functionality below.
    if (unbalanced_yxf) then
       
       call calculate_unbalanced_y(nproc, iproc, unbalanced_amount)
       
    end if
    
    if (.not. unbalanced_yxf) then
       
       yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
       yxf_lo%llim_proc = yxf_lo%blocksize*iproc
       yxf_lo%ulim_proc &
            = min(yxf_lo%ulim_world, yxf_lo%llim_proc + yxf_lo%blocksize - 1)
       yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)
       
    else
       
       if(proc0) then	
          write(*, fmt="('Using unbalanced decomposition for yxf. '&
               'Unbalanced fraction',F6.2)") unbalanced_amount
       end if
       
    end if
    
  end subroutine init_y_transform_layouts

  subroutine calculate_unbalanced_y(nproc, iproc,  unbalanced_amount)

    implicit none
    integer, intent (in) :: nproc, iproc
    real, intent (out) :: unbalanced_amount
    integer :: level_proc_num, i, j, k, m

    if(.not. initialized_y_transform) then
       write(*,*) 'Y Transform data structures not initialized so calculate_unbalanced_y subroutine will not operate correctly'
       write(*,*) 'Aborting subroutine calculate_unbalanced_y'
       return
    end if

    yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
    yxf_lo%small_block_balance_factor = 1
    yxf_lo%large_block_balance_factor = 1                                                        
    
    level_proc_num = nproc

    k = yxf_lo%nspec

    call calculate_block_breakdown(k, i, m, j, level_proc_num)
    
    if(j .eq. 0) then
       
       if(m .eq. 0) then
          level_proc_num = level_proc_num/k
          k = 1
       end if
       
       k = yxf_lo%nmu * k

       call calculate_block_breakdown(k, i, m, j, level_proc_num)
       
       if(j .eq. 0) then
          
          if(m .eq. 0) then
             level_proc_num = level_proc_num/k
             k = 1
          end if

          k = yxf_lo%nmu * k
          
          call calculate_block_breakdown(k, i, m, j, level_proc_num)
          
          if(j .eq. 0) then
             
             if(m .eq. 0) then
                level_proc_num = level_proc_num/k
                k = 1
             end if
             
             k = (2*yxf_lo%nvgrid+1) * k
             call calculate_block_breakdown(k, i, m, j, level_proc_num)
             
             if(j .eq. 0) then
                
                if(m .eq. 0) then
                   level_proc_num = level_proc_num/k
                   k = 1
                end if
                
                k = yxf_lo%ntgridtotal * k
                call calculate_block_breakdown(k, i, m, j, level_proc_num)
                
                if(j .eq. 0) then
                   
                   if(m .eq. 0) then
                      level_proc_num = level_proc_num/k
                      k = 1
                   end if
                   
                   k = yxf_lo%nx * k
                   call calculate_block_breakdown(k, i, m, j, level_proc_num)
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
                      call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, nproc, &
                           1, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                      
                   end if
                   
                else
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)     
                      call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, nproc, &
                           yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                      
                   end if
                   
                end if
                
             else
                
                if(i .ne. 0) then
                   
                   call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)                         
                   call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, nproc, &
                        yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, &
                        yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                   
                end if
                
             end if
             
             
          else
             
             if(i .ne. 0) then
                
                call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, nproc, &
                     (2*yxf_lo%nvgrid+1)*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, &
                     yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                
             end if
             
          end if
          
       else
          
          if(i .ne. 0) then
             
             call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)

             call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, nproc, &
                  yxf_lo%nmu*(2*yxf_lo%nvgrid+1)*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, &
                  yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)

          end if
          
       end if
       
    else
       
       if(i .ne. 0) then
          
          call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
          call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, nproc, &
               yxf_lo%nmu*(2*yxf_lo%nvgrid+1)*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, &
               yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
          
       end if
       
    end if
    
    
    if(yxf_lo%large_block_balance_factor .eq. 1 .and. yxf_lo%small_block_balance_factor .eq. 1) then
       
       unbalanced_amount = 0
       
    else 
       
       unbalanced_amount = real(yxf_lo%large_block_size)/real(yxf_lo%small_block_size)
       unbalanced_amount = unbalanced_amount - 1
              
    end if
    
    if (unbalanced_amount .gt. max_unbalanced_yxf .or. unbalanced_amount .eq. 0) then
       
       unbalanced_yxf = .false.
       
    end if
    

  end subroutine calculate_unbalanced_y

  elemental function is_idx_yxf (lo, i)
    implicit none
    integer :: is_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ! TT: the order of the division doesn't matter, so no need for branching
    is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/(2*lo%nvgrid+1)/lo%nmu, lo%nspec)

  end function is_idx_yxf

  elemental function imu_idx_yxf (lo, i)
    implicit none
    integer :: imu_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    imu_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/(2*lo%nvgrid + 1), lo%nmu)
  end function imu_idx_yxf

  elemental function iv_idx_yxf (lo, i)
    implicit none
    integer :: iv_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    iv_idx_yxf = -lo%nvgrid + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), 2*lo%nvgrid+1)
  end function iv_idx_yxf

  elemental function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
  end function ig_idx_yxf

  elemental function it_idx_yxf (lo, i)
    implicit none
    integer :: it_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
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
    integer :: block_offset, offset_block_number, j, k, tempi

    !AJ This code has been added to deal with the unbalanced decomposition functionality.
    !AJ If an unbalanced yxf decomposition is being used then the proc_id cannot
    !AJ use a simple lo%blocksize as there will be two separate block
    !AJ sizes used so we have to work out which block this i lives in and 
    !AJ therefore which process it belongs to.
    if (unbalanced_yxf) then
       !AJ block_offset works out how many groups of blocks this i point is
       !AJ after (the small and large blocks used in the decomposition are
       !AJ grouped together).
       block_offset = (i / lo%block_multiple)
       !AJ j represents how many blocks are in each group of blocks
       j = lo%num_small + lo%num_large
       !AJ offset_block_number is the number of blocks up to the start of the 
       !AJ group of blocks we are considering.
       offset_block_number = block_offset * j
       !AJ tempi represents where this index is inside the group of blocks
       !AJ this i point sits.
       tempi = i - (block_offset * lo%block_multiple)
       !AJ Work through each block in the group of blocks and see if this i 
       !AJ is within that block.  If it is set the proc_id_xxf as this block 
       !AJ owner.          
       do k=1,j
          if(k .le. lo%num_small) then
             !AJ TODO: The if-else construct used below could potentially be rationalised 
             !AJ TODO: to a more efficient formula where:
             !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%small_block_size
             !AJ TODO: Although a method for selecting if tempi is in the small or large 
             !AJ TODO: blocks would have to be provided.
             if(tempi .lt. lo%small_block_size) then
                !AJ (k -1) is the number of blocks that we have already considered 
                !AJ within this group of blocks we have selected. 
                proc_id_yxf = offset_block_number + (k - 1)
                exit 
             else
                !AJ If the index is not in this block then reduce tempi by a small
                !AJ block size and move on to the next block.
                tempi = tempi - lo%small_block_size
             end if
          else
             !AJ TODO: The if-else construct used below could potentially be rationalised 
             !AJ TODO: to a more efficient formula where:
             !AJ TODO: proc_id_xxf = offset_block_number + tempi/lo%large_block_size
             !AJ TODO: Although a method for selecting if tempi is in the small or large 
             !AJ TODO: blocks would have to be provided.
             if(tempi .lt. lo%large_block_size) then
                proc_id_yxf = (block_offset * j) + (k - 1)
                exit 
             else
                !AJ If the index is not in this block then reduce tempi by a large
                !AJ block size and move on to the next block.
                tempi = tempi - lo%large_block_size
             end if
          end if
       end do
    else
       proc_id_yxf = i/lo%blocksize
    end if

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
  subroutine gidx2xxfidx (ig, iv, iglo, g_lo, xxf_lo, it, ixxf)
# else
  elemental subroutine gidx2xxfidx (ig, iv, iglo, g_lo, xxf_lo, it, ixxf)
# endif

    implicit none
    integer, intent (in) :: ig, iv, iglo
    type (g_layout_type), intent (in) :: g_lo
    type (xxf_layout_type), intent (in) :: xxf_lo
    integer, intent (out) :: it, ixxf

    it = it_idx(g_lo,iglo)
    if (it > (xxf_lo%ntheta0+1)/2) then
       it = it - xxf_lo%ntheta0 + xxf_lo%nx
    end if

    ixxf = idx(xxf_lo, ig, iv, ik_idx(g_lo,iglo), &
         imu_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2xxfidx

! TT>
# ifdef USE_C_INDEX
  subroutine xxfidx2gidx (it, ixxf, xxf_lo, g_lo, ig, iv, iglo)
# else
  elemental subroutine xxfidx2gidx (it, ixxf, xxf_lo, g_lo, ig, iv, iglo)
# endif
! <TT
    implicit none
    integer, intent (in) :: it, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (g_layout_type), intent (in) :: g_lo
    integer, intent (out) :: ig, iv, iglo
    integer :: it0

    it0 = it
    if (it0 > (xxf_lo%ntheta0+1)/2) then
       it0 = it0 + xxf_lo%ntheta0 - xxf_lo%nx
       if (it0 <= (xxf_lo%ntheta0+1)/2) then
          ig = -999999
          iv = -999999
          iglo = -999999
          return
       end if
    end if

    ig = ig_idx(xxf_lo,ixxf)
    iv = iv_idx(xxf_lo,ixxf)
    iglo = idx(g_lo, ik_idx(xxf_lo,ixxf), it0, &
         imu_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
  end subroutine xxfidx2gidx

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

  subroutine factors (n, j, div)
    integer, intent (in) :: n
    integer, intent (out) :: j
    integer, dimension (:), intent (out) :: div
    integer :: i, imax

    do i=2,n
       if (mod(n,i)==0) exit
    end do
    imax = n/i
    j=1
    do i=1,imax
       if (mod(n,i)==0) then
          div(j) = i
          j=j+1
       end if
    end do
    div(j) = n
  end subroutine factors

  subroutine finish_layouts

    implicit none

    deallocate (f_lo, ij, mj, dj)

  end subroutine finish_layouts

end module gs2_layouts
