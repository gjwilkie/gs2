module gs2_layouts
  implicit none
  private

  public :: pe_layout

  public :: init_dist_fn_layouts
  public :: g_lo, g_layout_type
  public :: gint_lo, gint_layout_type
  public :: geint_lo, geint_layout_type

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_jfields_layouts
  public :: jf_lo, jf_layout_type
  public :: mj, ij, dj

!  public :: lambda_lo, lambda_layout_type
!  public :: init_lambda_layouts

  public :: init_lorentz_layouts
  public :: lz_lo, lz_layout_type
  public :: gidx2lzidx !, lzidx2gidx
!  public :: gidx2lamidx, lamidx2gintidx
  public :: gidx2gintidx, gintidx2geidx

  public :: init_x_transform_layouts, init_y_transform_layouts
  public :: xxf_lo, xxf_layout_type, yxf_lo, yxf_layout_type
  public :: gidx2xxfidx, xxfidx2yxfidx, yxfidx2xxfidx, xxfidx2gidx
  public :: xxf_ky_is_zero

  public :: init_accel_transform_layouts
  public :: accel_lo, accel_layout_type, dealiasing
  public :: accelx_lo, accelx_layout_type

  public :: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx, if_idx, isign_idx
  public :: im_idx, in_idx, ij_idx, ifield_idx
  public :: idx, proc_id, idx_local

  type :: g_layout_type
     integer :: iproc
     integer :: naky, ntheta0, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type g_layout_type

  type :: gint_layout_type
     integer :: iproc
     integer :: naky, ntheta0, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type gint_layout_type

  type :: geint_layout_type
     integer :: iproc
     integer :: naky, ntheta0, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type geint_layout_type

  type (g_layout_type) :: g_lo
  type (gint_layout_type) :: gint_lo
  type (geint_layout_type) :: geint_lo

  type :: f_layout_type
     integer :: iproc
     integer :: nidx, nfield, ntgrid, nindex, naky, ntheta0, M_class, N_class, i_class
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
     integer, dimension(:,:), pointer :: ik !=> null()
     integer, dimension(:,:), pointer :: it !=> null()
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

  type :: lz_layout_type
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type lz_layout_type

  type (lz_layout_type) :: lz_lo

!  type :: lambda_layout_type
!     integer :: iproc
!     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2
!     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!  end type lambda_layout_type

!  type (lambda_layout_type) :: lambda_lo

  type :: xxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ntheta0, nx, nadd, negrid, nlambda, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type xxf_layout_type

  type (xxf_layout_type) :: xxf_lo

  type :: yxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ny, ntheta0, nx, negrid, nlambda, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type yxf_layout_type

  type (yxf_layout_type) :: yxf_lo

  type :: accel_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ndky, ny, ntheta0
     integer :: nx, nxny, nxnky, negrid, nlambda, nspec, nia
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type accel_layout_type

  type (accel_layout_type) :: accel_lo

  type :: accelx_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ny, ntheta0, nx, nxny, negrid, nlambda, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type accelx_layout_type

  type (accelx_layout_type) :: accelx_lo

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
     module procedure ig_idx_lz
     module procedure ig_idx_xxf
     module procedure ig_idx_yxf
     module procedure ig_idx_f
  end interface

  interface ik_idx
     module procedure ik_idx_g
     module procedure ik_idx_gint
     module procedure ik_idx_geint
     module procedure ik_idx_jf
     module procedure ik_idx_lz
!     module procedure ik_idx_lambda
     module procedure ik_idx_xxf
     module procedure ik_idx_accel
  end interface

  interface it_idx
     module procedure it_idx_g
     module procedure it_idx_gint
     module procedure it_idx_geint
     module procedure it_idx_jf
     module procedure it_idx_lz
!     module procedure it_idx_lambda
     module procedure it_idx_yxf
     module procedure it_idx_accel
  end interface

  interface il_idx
     module procedure il_idx_g
     module procedure il_idx_xxf
     module procedure il_idx_yxf
  end interface

  interface ie_idx
     module procedure ie_idx_g
     module procedure ie_idx_gint
     module procedure ie_idx_lz
!     module procedure ie_idx_lambda
     module procedure ie_idx_xxf
     module procedure ie_idx_yxf
  end interface

  interface is_idx
     module procedure is_idx_g
     module procedure is_idx_gint
     module procedure is_idx_geint
     module procedure is_idx_lz
!     module procedure is_idx_lambda
     module procedure is_idx_xxf
     module procedure is_idx_yxf
  end interface

  interface isign_idx
     module procedure isign_idx_xxf
     module procedure isign_idx_yxf
  end interface

  interface proc_id
     module procedure proc_id_g
     module procedure proc_id_gint
     module procedure proc_id_geint
     module procedure proc_id_f
     module procedure proc_id_jf
     module procedure proc_id_lz
!     module procedure proc_id_lambda
     module procedure proc_id_xxf
     module procedure proc_id_yxf
  end interface

  interface idx
     module procedure idx_g
     module procedure idx_gint
     module procedure idx_geint
     module procedure idx_f
     module procedure idx_jf
     module procedure idx_lz
     module procedure idx_xxf
     module procedure idx_yxf
!     module procedure idx_lambda
  end interface

  interface idx_local
     module procedure idx_local_g,      ig_local_g
     module procedure idx_local_gint,   ig_local_gint
     module procedure idx_local_geint,  ig_local_geint
     module procedure idx_local_f,      ig_local_f
     module procedure idx_local_jf,     ig_local_jf
     module procedure idx_local_lz,     ig_local_lz
     module procedure idx_local_xxf,    ig_local_xxf
     module procedure idx_local_yxf,    ig_local_yxf
!     module procedure idx_local_lambda, ig_local_lambda
 end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_dist_fn_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    g_lo%iproc = iproc
    g_lo%naky = naky
    g_lo%ntheta0 = ntheta0
    g_lo%nlambda = nlambda
    g_lo%negrid = negrid
    g_lo%nspec = nspec
    g_lo%llim_world = 0
    g_lo%ulim_world = nlambda*negrid*ntheta0*naky*nspec - 1
    g_lo%blocksize = g_lo%ulim_world/nproc + 1
    g_lo%llim_proc = g_lo%blocksize*iproc
    g_lo%ulim_proc = min(g_lo%ulim_world, g_lo%llim_proc + g_lo%blocksize - 1)
    g_lo%ulim_alloc = max(g_lo%llim_proc, g_lo%ulim_proc)

    gint_lo%iproc = iproc
    gint_lo%naky = naky
    gint_lo%ntheta0 = ntheta0
    gint_lo%negrid = negrid
    gint_lo%nspec = nspec
    gint_lo%llim_world = 0
    gint_lo%ulim_world = negrid*ntheta0*naky*nspec - 1
    gint_lo%blocksize = gint_lo%ulim_world/nproc + 1
    gint_lo%llim_proc = gint_lo%blocksize*iproc
    gint_lo%ulim_proc &
         = min(gint_lo%ulim_world, gint_lo%llim_proc + gint_lo%blocksize - 1)
    gint_lo%ulim_alloc = max(gint_lo%llim_proc, gint_lo%ulim_proc)
    
    geint_lo%iproc = iproc
    geint_lo%naky = naky
    geint_lo%ntheta0 = ntheta0
    geint_lo%nspec = nspec
    geint_lo%llim_world = 0
    geint_lo%ulim_world = ntheta0*naky*nspec - 1
    geint_lo%blocksize = geint_lo%ulim_world/nproc + 1
    geint_lo%llim_proc = geint_lo%blocksize*iproc
    geint_lo%ulim_proc &
         = min(geint_lo%ulim_world, geint_lo%llim_proc + geint_lo%blocksize -1)
    geint_lo%ulim_alloc = max(geint_lo%llim_proc, geint_lo%ulim_proc)
  end subroutine init_dist_fn_layouts

  pure function is_idx_g (lo, i)
    implicit none
    integer :: is_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
  end function is_idx_g

  pure function ik_idx_g (lo, i)
    implicit none
    integer :: ik_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%negrid/lo%ntheta0, lo%naky)
  end function ik_idx_g

  pure function it_idx_g (lo, i)
    implicit none
    integer :: it_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%negrid, lo%ntheta0)
  end function it_idx_g

  pure function ie_idx_g (lo, i)
    implicit none
    integer :: ie_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda, lo%negrid)
  end function ie_idx_g

  pure function il_idx_g (lo, i)
    implicit none
    integer :: il_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_g = 1 + mod(i - lo%llim_world, lo%nlambda)
  end function il_idx_g

  pure function proc_id_g (lo, i)
    implicit none
    integer :: proc_id_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_g = i/lo%blocksize
  end function proc_id_g

  pure function idx_g (lo, ik, it, il, ie, is)
    implicit none
    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    idx_g = il-1 + lo%nlambda*(ie-1 + lo%negrid*(it-1 + lo%ntheta0*(ik-1 &
         + lo%naky*(is-1))))
  end function idx_g

  pure function idx_local_g (lo, ik, it, il, ie, is)
    implicit none
    logical :: idx_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is

    idx_local_g = idx_local(lo, idx(lo, ik, it, il, ie, is))
  end function idx_local_g

  pure function ig_local_g (lo, ig)
    implicit none
    logical :: ig_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_g = lo%iproc == proc_id(lo, ig)
  end function ig_local_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Once-integrated distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function is_idx_gint (lo, i)
    implicit none
    integer :: is_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_gint = 1 + mod((i - lo%llim_world)/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
  end function is_idx_gint

  pure function ik_idx_gint (lo, i)
    implicit none
    integer :: ik_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_gint = 1 + mod((i - lo%llim_world)/lo%negrid/lo%ntheta0, lo%naky)
  end function ik_idx_gint

  pure function it_idx_gint (lo, i)
    implicit none
    integer :: it_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_gint = 1 + mod((i - lo%llim_world)/lo%negrid, lo%ntheta0)
  end function it_idx_gint

  pure function ie_idx_gint (lo, i)
    implicit none
    integer :: ie_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_gint = 1 + mod(i - lo%llim_world, lo%negrid)
  end function ie_idx_gint

  pure function proc_id_gint (lo, i)
    implicit none
    integer :: proc_id_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_gint = i/lo%blocksize
  end function proc_id_gint

  pure function idx_gint (lo, ik, it, ie, is)
    implicit none
    integer :: idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, ie, is
    idx_gint = ie-1 + lo%negrid*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
  end function idx_gint

  pure function idx_local_gint (lo, ik, it, ie, is)
    implicit none
    logical :: idx_local_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, ie, is

    idx_local_gint = idx_local(lo, idx(lo, ik, it, ie, is))
  end function idx_local_gint

  pure function ig_local_gint (lo, ig)
    implicit none
    logical :: ig_local_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_gint = lo%iproc == proc_id(lo, ig)
  end function ig_local_gint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Twice-integrated distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function is_idx_geint (lo, i)
    implicit none
    integer :: is_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%nspec)
  end function is_idx_geint

  pure function ik_idx_geint (lo, i)
    implicit none
    integer :: ik_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
  end function ik_idx_geint

  pure function it_idx_geint (lo, i)
    implicit none
    integer :: it_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_geint = 1 + mod(i - lo%llim_world, lo%ntheta0)
  end function it_idx_geint

  pure function proc_id_geint (lo, i)
    implicit none
    integer :: proc_id_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_geint = i/lo%blocksize
  end function proc_id_geint

  pure function idx_geint (lo, ik, it, is)
    implicit none
    integer :: idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, is
    idx_geint = it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1))
  end function idx_geint

  pure function idx_local_geint (lo, ik, it, is)
    implicit none
    logical :: idx_local_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, is

    idx_local_geint = idx_local(lo, idx(lo, ik, it, is))
  end function idx_local_geint

  pure function ig_local_geint (lo, ig)
    implicit none
    logical :: ig_local_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_geint = lo%iproc == proc_id(lo, ig)
  end function ig_local_geint

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
    integer :: i

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
       f_lo(i)%blocksize = f_lo(i)%ulim_world/nproc + 1
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
! Lorentz layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_lorentz_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    lz_lo%iproc = iproc
    lz_lo%ntgrid = ntgrid
    lz_lo%naky = naky
    lz_lo%ntheta0 = ntheta0
    lz_lo%negrid = negrid
    lz_lo%nspec = nspec
    lz_lo%ng2 = ng2
    lz_lo%llim_world = 0
    lz_lo%ulim_world = (2*ntgrid+1)*negrid*ntheta0*naky*nspec - 1
    lz_lo%blocksize = lz_lo%ulim_world/nproc + 1
    lz_lo%llim_proc = lz_lo%blocksize*iproc
    lz_lo%ulim_proc &
         = min(lz_lo%ulim_world, lz_lo%llim_proc + lz_lo%blocksize - 1)
    lz_lo%ulim_alloc = max(lz_lo%llim_proc, lz_lo%ulim_proc)
  end subroutine init_lorentz_layouts

  pure function is_idx_lz (lo, i)
    implicit none
    integer :: is_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_lz = 1 + mod((i - lo%llim_world) &
         /(2*lo%ntgrid + 1)/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
  end function is_idx_lz

  pure function ik_idx_lz (lo, i)
    implicit none
    integer :: ik_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_lz = 1 + mod((i - lo%llim_world) &
         /(2*lo%ntgrid + 1)/lo%negrid/lo%ntheta0, lo%naky)
  end function ik_idx_lz

  pure function it_idx_lz (lo, i)
    implicit none
    integer :: it_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_lz = 1 + mod((i - lo%llim_world) &
         /(2*lo%ntgrid + 1)/lo%negrid, lo%ntheta0)
  end function it_idx_lz

  pure function ie_idx_lz (lo, i)
    implicit none
    integer :: ie_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%negrid)
  end function ie_idx_lz

  pure function ig_idx_lz (lo, i)
    implicit none
    integer :: ig_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
  end function ig_idx_lz

  pure function idx_lz (lo, ig, ik, it, ie, is)
    implicit none
    integer :: idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is
    idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ie-1 + lo%negrid*(it-1 &
         + lo%ntheta0*(ik-1 + lo%naky*(is-1))))
  end function idx_lz

  pure function proc_id_lz (lo, i)
    implicit none
    integer :: proc_id_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_lz = i/lo%blocksize
  end function proc_id_lz

  pure function idx_local_lz (lo, ig, ik, it, ie, is)
    implicit none
    logical :: idx_local_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is

    idx_local_lz = idx_local(lo, idx(lo, ig, ik, it, ie, is))
  end function idx_local_lz

  pure function ig_local_lz (lo, ig)
    implicit none
    logical :: ig_local_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_lz = lo%iproc == proc_id(lo, ig)
  end function ig_local_lz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lambda layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine init_lambda_layouts &
!       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)       
!    use mp, only: iproc, nproc
!    implicit none
!    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2
!    logical :: initialized = .false.
!    
!    if (initialized) return
!    initialized = .true.
!
!    lambda_lo%iproc = iproc
!    lambda_lo%ntgrid = ntgrid
!    lambda_lo%naky = naky
!    lambda_lo%ntheta0 = ntheta0
!    lambda_lo%negrid = negrid
!    lambda_lo%nspec = nspec
!    lambda_lo%ng2 = ng2
!    lambda_lo%llim_world = 0
!    lambda_lo%ulim_world = negrid*ntheta0*naky*nspec - 1
!    lambda_lo%blocksize = lambda_lo%ulim_world/nproc + 1
!    lambda_lo%llim_proc = lambda_lo%blocksize*iproc
!    lambda_lo%ulim_proc &
!         = min(lambda_lo%ulim_world, lambda_lo%llim_proc + lambda_lo%blocksize - 1)
!    lambda_lo%ulim_alloc = max(lambda_lo%llim_proc, lambda_lo%ulim_proc)
!
!  end subroutine init_lambda_layouts
!
!  pure function is_idx_lambda (lo, i)
!    implicit none
!    integer :: is_idx_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: i
!    is_idx_lambda = 1 + mod((i - lo%llim_world)/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
!  end function is_idx_lambda
!
!  pure function ik_idx_lambda (lo, i)
!    implicit none
!    integer :: ik_idx_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: i
!    ik_idx_lambda = 1 + mod((i - lo%llim_world)/lo%negrid/lo%ntheta0, lo%naky)
!  end function ik_idx_lambda
!
!  pure function it_idx_lambda (lo, i)
!    implicit none
!    integer :: it_idx_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: i
!    it_idx_lambda = 1 + mod((i - lo%llim_world)/lo%negrid, lo%ntheta0)
!  end function it_idx_lambda
!
!  pure function ie_idx_lambda (lo, i)
!    implicit none
!    integer :: ie_idx_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: i
!    ie_idx_lambda = 1 + mod(i - lo%llim_world, lo%negrid)
!  end function ie_idx_lambda
!
!  pure function idx_lambda (lo, ik, it, ie, is)
!    implicit none
!    integer :: idx_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: ik, it, ie, is
!    idx_lambda = ie-1 + lo%negrid*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
!  end function idx_lambda
!
!  pure function proc_id_lambda (lo, i)
!    implicit none
!    integer :: proc_id_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: i
!    proc_id_lambda = i/lo%blocksize
!  end function proc_id_lambda
!
!  pure function idx_local_lambda (lo, ik, it, ie, is)
!    implicit none
!    logical :: idx_local_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: ik, it, ie, is
!
!    idx_local_lambda = idx_local(lo, idx(lo, ik, it, ie, is))
!  end function idx_local_lambda
!
!  pure function ig_local_lambda (lo, ig)
!    implicit none
!    logical :: ig_local_lambda
!    type (lambda_layout_type), intent (in) :: lo
!    integer, intent (in) :: ig
!
!    ig_local_lambda = lo%iproc == proc_id(lo, ig)
!  end function ig_local_lambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_x_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    xxf_lo%iproc = iproc
    xxf_lo%ntgrid = ntgrid
    xxf_lo%nsign = 2
    xxf_lo%naky = naky
    xxf_lo%ntheta0 = ntheta0
    if (nx > ntheta0) then
       xxf_lo%nx = nx
    else
       xxf_lo%nx = (3*ntheta0+1)/2
    end if
    xxf_lo%nadd = xxf_lo%nx - ntheta0
    xxf_lo%nlambda = nlambda
    xxf_lo%negrid = negrid
    xxf_lo%nspec = nspec
    xxf_lo%llim_world = 0
    xxf_lo%ulim_world = naky*(2*ntgrid+1)*2*nlambda*negrid*nspec - 1
    xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
    xxf_lo%llim_proc = xxf_lo%blocksize*iproc
    xxf_lo%ulim_proc &
         = min(xxf_lo%ulim_world, xxf_lo%llim_proc + xxf_lo%blocksize - 1)
    xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)
  end subroutine init_x_transform_layouts

  pure function is_idx_xxf (lo, i)
    implicit none
    integer :: is_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
  end function is_idx_xxf

  pure function ie_idx_xxf (lo, i)
    implicit none
    integer :: ie_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
  end function ie_idx_xxf

  pure function il_idx_xxf (lo, i)
    implicit none
    integer :: il_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
  end function il_idx_xxf

  pure function isign_idx_xxf (lo, i)
    implicit none
    integer :: isign_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    isign_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), lo%nsign)
  end function isign_idx_xxf

  pure function ig_idx_xxf (lo, i)
    implicit none
    integer :: ig_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
  end function ig_idx_xxf

  pure function ik_idx_xxf (lo, i)
    implicit none
    integer :: ik_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
  end function ik_idx_xxf

  pure function idx_xxf (lo, ig, isign, ik, il, ie, is)
    implicit none
    integer :: idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, il, ie, is
    idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
         + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
  end function idx_xxf

  pure function proc_id_xxf (lo, i)
    implicit none
    integer :: proc_id_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_xxf = i/lo%blocksize
  end function proc_id_xxf

  pure function idx_local_xxf (lo, ig, isign, ik, il, ie, is)
    implicit none
    logical :: idx_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, il, ie, is
    idx_local_xxf = idx_local (lo, idx(lo, ig, isign, ik, il, ie, is))
  end function idx_local_xxf

  pure function ig_local_xxf (lo, i)
    implicit none
    logical ig_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_xxf = lo%iproc == proc_id(lo, i)
  end function ig_local_xxf

  pure function xxf_ky_is_zero (lo, i)
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
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    logical, save :: initialized = .false.
    integer :: nnx, nny

    if (initialized) return
    initialized = .true.

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
    yxf_lo%nsign = 2
    yxf_lo%naky = naky
    yxf_lo%ny = nny
    yxf_lo%ntheta0 = ntheta0
    yxf_lo%nx = nnx
    yxf_lo%nlambda = nlambda
    yxf_lo%negrid = negrid
    yxf_lo%nspec = nspec
    yxf_lo%llim_world = 0
    yxf_lo%ulim_world = nnx*(2*ntgrid+1)*2*nlambda*negrid*nspec - 1
    yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
    yxf_lo%llim_proc = yxf_lo%blocksize*iproc
    yxf_lo%ulim_proc &
         = min(yxf_lo%ulim_world, yxf_lo%llim_proc + yxf_lo%blocksize - 1)
    yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)
  end subroutine init_y_transform_layouts

  pure function is_idx_yxf (lo, i)
    implicit none
    integer :: is_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
  end function is_idx_yxf

  pure function ie_idx_yxf (lo, i)
    implicit none
    integer :: ie_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
  end function ie_idx_yxf

  pure function il_idx_yxf (lo, i)
    implicit none
    integer :: il_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
  end function il_idx_yxf

  pure function isign_idx_yxf (lo, i)
    implicit none
    integer :: isign_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    isign_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), lo%nsign)
  end function isign_idx_yxf

  pure function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
  end function ig_idx_yxf

  pure function it_idx_yxf (lo, i)
    implicit none
    integer :: it_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
  end function it_idx_yxf

  pure function idx_yxf (lo, ig, isign, it, il, ie, is)
    implicit none
    integer :: idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, it, il, ie, is
    idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
         + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
  end function idx_yxf

  pure function proc_id_yxf (lo, i)
    implicit none
    integer :: proc_id_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_yxf = i/lo%blocksize
  end function proc_id_yxf

  pure function idx_local_yxf (lo, ig, isign, it, il, ie, is)
    implicit none
    logical :: idx_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, it, il, ie, is
    idx_local_yxf = idx_local (lo, idx(lo, ig, isign, it, il, ie, is))
  end function idx_local_yxf

  pure function ig_local_yxf (lo, i)
    implicit none
    logical ig_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_yxf = lo%iproc == proc_id(lo, i)
  end function ig_local_yxf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Accelerated FFT layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_accel_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    logical, save :: initialized = .false.
    integer :: nnx, nny

    if (initialized) return
    initialized = .true.

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

    accelx_lo%iproc = iproc
    accelx_lo%ntgrid = ntgrid
    accelx_lo%nsign = 2
    accelx_lo%naky = naky
    accelx_lo%ny = nny
    accelx_lo%ntheta0 = ntheta0
    accelx_lo%nx = nnx
    accelx_lo%nxny = nnx*nny
    accelx_lo%nlambda = nlambda
    accelx_lo%negrid = negrid
    accelx_lo%nspec = nspec
    accelx_lo%llim_world = 0
    accelx_lo%ulim_world = nnx*nny*nlambda*negrid*nspec - 1
    accelx_lo%blocksize = accelx_lo%ulim_world/nproc + 1
    accelx_lo%llim_proc = accelx_lo%blocksize*iproc
    accelx_lo%ulim_proc &
         = min(accelx_lo%ulim_world, accelx_lo%llim_proc + accelx_lo%blocksize - 1)
    accelx_lo%ulim_alloc = max(accelx_lo%llim_proc, accelx_lo%ulim_proc)

    accel_lo%iproc = iproc
    accel_lo%ntgrid = ntgrid
    accel_lo%nsign = 2
    accel_lo%naky = naky
    accel_lo%ndky = nny/2+1  ! Note: this is for dealiased k space quantities
    accel_lo%ny = nny
    accel_lo%ntheta0 = ntheta0
    accel_lo%nx = nnx
    accel_lo%nxny = nnx*nny
    accel_lo%nxnky = nnx*(nny/2+1)
    accel_lo%nlambda = nlambda
    accel_lo%negrid = negrid
    accel_lo%nspec = nspec
    accel_lo%nia = nlambda*negrid*nspec/nproc
    accel_lo%llim_world = 0
    accel_lo%ulim_world = nnx*(nny/2+1)*nlambda*negrid*nspec - 1
    accel_lo%blocksize = accel_lo%ulim_world/nproc + 1
    accel_lo%llim_proc = accel_lo%blocksize*iproc
    accel_lo%ulim_proc &
         = min(accel_lo%ulim_world, accel_lo%llim_proc + accel_lo%blocksize - 1)
    accel_lo%ulim_alloc = max(accel_lo%llim_proc, accel_lo%ulim_proc)

  end subroutine init_accel_transform_layouts

  pure function it_idx_accel (lo, i)
    implicit none
    integer :: it_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
  end function it_idx_accel

  pure function ik_idx_accel (lo, i)
    implicit none
    integer :: ik_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
  end function ik_idx_accel

  pure function dealiasing (lo, ia)
    implicit none
    logical dealiasing
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: ia
    
    dealiasing = .true.

    if (it_idx(accel_lo, ia) > lo%ntheta0/2+1 &
  .and. it_idx(accel_lo, ia) <= lo%nx-lo%ntheta0/2) return
    if (ik_idx(accel_lo, ia) > lo%naky) return

    dealiasing = .false.
  end function dealiasing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transformation subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, &
                              il, ilz)
    implicit none
    integer, intent (in) :: ig, isign
    type (g_layout_type), intent (in) :: g_lo
    integer, intent (in) :: iglo
    type (lz_layout_type), intent (in) :: lz_lo
    integer, intent (in) :: ntgrid
    integer, dimension (-ntgrid:), intent (in) :: jend
    integer, intent (out) :: il, ilz
    integer :: je

    je = jend(ig)
    il = il_idx(g_lo,iglo)

    if (je == 0) then
       if (isign == 2) then
          il = 2*g_lo%nlambda+1 - il
       end if
    else

!       if (il > je) then
!          il = 2*je + 1
!       else if (isign == 2) then
!          il = 2*je - il 
!       end if

       if (il == je) then
          if (isign == 1) il = 2*je  ! throw this info away
       else if (il > je) then 
          if (isign == 1) il = il + je
          if (isign == 2) il = 2*g_lo%nlambda + 1 - il + je
       else
          if (isign == 2) il = 2*je - il !+ 1
       end if
          
    end if
    ilz = idx(lz_lo, ig, ik_idx(g_lo,iglo), it_idx(g_lo,iglo), &
         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2lzidx

!  pure subroutine gidx2lamidx (g_lo, iglo, lambda_lo, il, ilam)
!    implicit none
!    type (g_layout_type), intent (in) :: g_lo
!    integer, intent (in) :: iglo
!    type (lambda_layout_type), intent (in) :: lambda_lo
!    integer, intent (out) :: il, ilam
!
!    il = il_idx(g_lo, iglo)
!
!    ilam = idx(lambda_lo, ik_idx(g_lo,iglo), it_idx(g_lo,iglo), &
!         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
!  end subroutine gidx2lamidx

  pure subroutine gidx2gintidx (g_lo, iglo, gint_lo, il, igint)
    implicit none
    type (g_layout_type), intent (in) :: g_lo
    integer, intent (in) :: iglo
    type (gint_layout_type), intent (in) :: gint_lo
    integer, intent (out) :: il, igint

    il = il_idx(g_lo, iglo)

    igint = idx(gint_lo, ik_idx(g_lo,iglo), it_idx(g_lo,iglo), &
         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2gintidx

  pure subroutine gintidx2geidx (gint_lo, igint, ie, geint_lo, igeint)
    implicit none
    type (gint_layout_type), intent (in) :: gint_lo
    integer, intent (in) :: igint
    type (geint_layout_type), intent (in) :: geint_lo
    integer, intent (out) :: ie, igeint

    ie = ie_idx(gint_lo, igint)

    igeint = idx(geint_lo, ik_idx(gint_lo, igint), it_idx(gint_lo, igint), &
         is_idx(gint_lo, igint))

  end subroutine gintidx2geidx

!  pure subroutine lamidx2gintidx (lambda_lo, ilam, gint_lo, igint)
!    implicit none
!    type (lambda_layout_type), intent (in) :: lambda_lo
!    integer, intent (in) :: ilam
!    type (gint_layout_type), intent (in) :: gint_lo
!    integer, intent (out) :: igint
!
!    igint = idx(gint_lo, ik_idx(lambda_lo, ilam), it_idx(lambda_lo, ilam), &
!         ie_idx(lambda_lo, ilam), is_idx(lambda_lo, ilam))
!  end subroutine lamidx2gintidx

! not used 
!  pure subroutine lzidx2gidx (il, ilz, lz_lo, g_lo, jend, ig, isign, iglo)
!    implicit none
!    integer, intent (in) :: il, ilz
!    type (lz_layout_type), intent (in) :: lz_lo
!    type (g_layout_type), intent (in) :: g_lo
!    integer, dimension (:), intent (in) :: jend
!    integer, intent (out) :: ig, isign, iglo
!    integer :: ilfold, je
!
!    ig = ig_idx(lz_lo,ilz)
!    je = jend(ig)
!
!    if (je == 0) then
!       je = lz_lo%ng2
!       if (il <= je) then
!          isign = 1
!          ilfold = il
!       else if (il > je+1) then
!          isign = 2
!          ilfold = 2*je + 1 - il
!       else
!          isign = -999999
!          ilfold = -999999
!       end if
!    else
!       if (il <= je) then
!          isign = 1
!          ilfold = il
!       else if (il <= 2*je) then
!          isign = 2
!          ilfold = 2*je - il + 1
!!          ilfold = 2*je - il 
!       else
!          isign = -999999
!          ilfold = -999999
!       end if
!    end if
!
!    iglo = idx(g_lo, ik_idx(lz_lo,ilz), it_idx(lz_lo,ilz), &
!         ilfold, ie_idx(lz_lo,ilz), is_idx(lz_lo,ilz))
!  end subroutine lzidx2gidx

  pure subroutine gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
    implicit none
    integer, intent (in) :: ig, isign, iglo
    type (g_layout_type), intent (in) :: g_lo
    type (xxf_layout_type), intent (in) :: xxf_lo
    integer, intent (out) :: it, ixxf

    it = it_idx(g_lo,iglo)
    if (it > (xxf_lo%ntheta0+1)/2) then
       it = it - xxf_lo%ntheta0 + xxf_lo%nx
    end if

    ixxf = idx(xxf_lo, ig, isign, ik_idx(g_lo,iglo), il_idx(g_lo,iglo), &
         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2xxfidx

  pure subroutine xxfidx2gidx (it, ixxf, xxf_lo, g_lo, ig, isign, iglo)
    implicit none
    integer, intent (in) :: it, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (g_layout_type), intent (in) :: g_lo
    integer, intent (out) :: ig, isign, iglo
    integer :: it0

    it0 = it
    if (it0 > (xxf_lo%ntheta0+1)/2) then
       it0 = it0 + xxf_lo%ntheta0 - xxf_lo%nx
       if (it0 <= (xxf_lo%ntheta0+1)/2) then
          ig = -999999
          isign = -999999
          iglo = -999999
          return
       end if
    end if

    ig = ig_idx(xxf_lo,ixxf)
    isign = isign_idx(xxf_lo,ixxf)
    iglo = idx(g_lo, ik_idx(xxf_lo,ixxf), it0, il_idx(xxf_lo,ixxf), &
         ie_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
  end subroutine xxfidx2gidx

  pure subroutine xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
    implicit none
    integer, intent (in) :: it, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (yxf_layout_type), intent (in) :: yxf_lo
    integer, intent (out) :: ik, iyxf

    ik = ik_idx(xxf_lo,ixxf)
    iyxf = idx(yxf_lo, ig_idx(xxf_lo,ixxf), isign_idx(xxf_lo,ixxf), &
         it, il_idx(xxf_lo,ixxf), ie_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
  end subroutine xxfidx2yxfidx

  pure subroutine yxfidx2xxfidx (ik, iyxf, yxf_lo, xxf_lo, it, ixxf)
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
    ixxf = idx(xxf_lo, ig_idx(yxf_lo,iyxf), isign_idx(yxf_lo,iyxf), &
         ik0, il_idx(yxf_lo,iyxf), ie_idx(yxf_lo,iyxf), is_idx(yxf_lo,iyxf))
  end subroutine yxfidx2xxfidx

  subroutine pe_layout (char)

    character (1), intent (out) :: char

    char = 'x'

  end subroutine pe_layout

end module gs2_layouts

