module gs2_layouts
  implicit none
  private

  public :: init_dist_fn_layouts
  public :: g_lo, g_layout_type
  public :: gint_lo, gint_layout_type
  public :: geint_lo, geint_layout_type

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_lorentz_layouts
  public :: lz_lo, lz_layout_type
  public :: gidx2lzidx, lzidx2gidx

  public :: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx, if_idx, idx
  public :: proc_id, idx_local

  type :: g_layout_type
     integer :: iproc
     integer :: naky, ntheta0, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, blocksize
  end type g_layout_type

  type :: gint_layout_type
     integer :: iproc
     integer :: naky, ntheta0, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, blocksize
  end type gint_layout_type

  type :: geint_layout_type
     integer :: iproc
     integer :: naky, ntheta0, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, blocksize
  end type geint_layout_type

  type (g_layout_type) :: g_lo
  type (gint_layout_type) :: gint_lo
  type (geint_layout_type) :: geint_lo

  type :: f_layout_type
     integer :: iproc
     integer :: nindex, naky, ntheta0
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, blocksize
  end type f_layout_type

  type (f_layout_type) :: f_lo

  type :: lz_layout_type
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, blocksize
  end type lz_layout_type

  type (lz_layout_type) :: lz_lo

  interface if_idx
     module procedure if_idx_f
  end interface

  interface ig_idx
     module procedure ig_idx_lz
  end interface

  interface ik_idx
     module procedure ik_idx_g
     module procedure ik_idx_gint
     module procedure ik_idx_geint
     module procedure ik_idx_f
     module procedure ik_idx_lz
  end interface

  interface it_idx
     module procedure it_idx_g
     module procedure it_idx_gint
     module procedure it_idx_geint
     module procedure it_idx_f
     module procedure it_idx_lz
  end interface

  interface il_idx
     module procedure il_idx_g
  end interface

  interface ie_idx
     module procedure ie_idx_g
     module procedure ie_idx_gint
     module procedure ie_idx_lz
  end interface

  interface is_idx
     module procedure is_idx_g
     module procedure is_idx_gint
     module procedure is_idx_geint
     module procedure is_idx_lz
  end interface

  interface proc_id
     module procedure proc_id_g
     module procedure proc_id_gint
     module procedure proc_id_geint
     module procedure proc_id_f
     module procedure proc_id_lz
  end interface

  interface idx
     module procedure idx_g
     module procedure idx_gint
     module procedure idx_geint
     module procedure idx_f
     module procedure idx_lz
  end interface

  interface idx_local
     module procedure idx_local_g,     ig_local_g
     module procedure idx_local_gint,  ig_local_gint
     module procedure idx_local_geint, ig_local_geint
     module procedure idx_local_f,     ig_local_f
     module procedure idx_local_lz,    ig_local_lz
  end interface

contains

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
    g_lo%ulim_world = naky*ntheta0*nlambda*negrid*nspec - 1
    g_lo%blocksize = g_lo%ulim_world/nproc + 1
    g_lo%llim_proc = g_lo%blocksize*iproc
    g_lo%ulim_proc = min(g_lo%ulim_world, g_lo%llim_proc + g_lo%blocksize - 1)

    gint_lo%iproc = iproc
    gint_lo%naky = naky
    gint_lo%ntheta0 = ntheta0
    gint_lo%negrid = negrid
    gint_lo%nspec = nspec
    gint_lo%llim_world = 0
    gint_lo%ulim_world = naky*ntheta0*negrid*nspec - 1
    gint_lo%blocksize = gint_lo%ulim_world/nproc + 1
    gint_lo%llim_proc = gint_lo%blocksize*iproc
    gint_lo%ulim_proc &
         = min(gint_lo%ulim_world, gint_lo%llim_proc + gint_lo%blocksize - 1)
    
    geint_lo%iproc = iproc
    geint_lo%naky = naky
    geint_lo%ntheta0 = ntheta0
    geint_lo%nspec = nspec
    geint_lo%llim_world = 0
    geint_lo%ulim_world = naky*ntheta0*nspec - 1
    geint_lo%blocksize = geint_lo%ulim_world/nproc + 1
    geint_lo%llim_proc = geint_lo%blocksize*iproc
    geint_lo%ulim_proc &
         = min(geint_lo%ulim_world, geint_lo%llim_proc + geint_lo%blocksize -1)
  end subroutine init_dist_fn_layouts

  subroutine init_fields_layouts (nindex, naky, ntheta0)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: nindex, naky, ntheta0
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    f_lo%iproc = iproc
    f_lo%nindex = nindex
    f_lo%naky = naky
    f_lo%ntheta0 = ntheta0
    f_lo%llim_world = 0
    f_lo%ulim_world = nindex*naky*ntheta0 - 1
    f_lo%blocksize = f_lo%ulim_world/nproc + 1
    f_lo%llim_proc = f_lo%blocksize*iproc
    f_lo%ulim_proc = min(f_lo%ulim_world, f_lo%llim_proc + f_lo%blocksize - 1)
  end subroutine init_fields_layouts

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
    lz_lo%ulim_world = (2*ntgrid+1)*naky*ntheta0*negrid*nspec - 1
    lz_lo%blocksize = lz_lo%ulim_world/nproc + 1
    lz_lo%llim_proc = lz_lo%blocksize*iproc
    lz_lo%ulim_proc &
         = min(lz_lo%ulim_world, lz_lo%llim_proc + lz_lo%blocksize - 1)
  end subroutine init_lorentz_layouts

  pure function if_idx_f (lo, i)
    implicit none
    integer :: if_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    if_idx_f = 1 + mod(i - lo%llim_world, lo%nindex)
  end function if_idx_f

  pure function ig_idx_lz (lo, i)
    implicit none
    integer :: ig_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
  end function ig_idx_lz

  pure function ik_idx_g (lo, i)
    implicit none
    integer :: ik_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_g = 1 + mod(i - lo%llim_world, lo%naky)
  end function ik_idx_g

  pure function ik_idx_gint (lo, i)
    implicit none
    integer :: ik_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_gint = 1 + mod(i - lo%llim_world, lo%naky)
  end function ik_idx_gint

  pure function ik_idx_geint (lo, i)
    implicit none
    integer :: ik_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_geint = 1 + mod(i - lo%llim_world, lo%naky)
  end function ik_idx_geint

  pure function ik_idx_f (lo, i)
    implicit none
    integer :: ik_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_f = 1 + mod((i - lo%llim_world)/lo%nindex, lo%naky)
  end function ik_idx_f

  pure function ik_idx_lz (lo, i)
    implicit none
    integer :: ik_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%naky)
  end function ik_idx_lz

  pure function it_idx_g (lo, i)
    implicit none
    integer :: it_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_g = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
  end function it_idx_g

  pure function it_idx_gint (lo, i)
    implicit none
    integer :: it_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
  end function it_idx_gint

  pure function it_idx_geint (lo, i)
    implicit none
    integer :: it_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
  end function it_idx_geint

  pure function it_idx_f (lo, i)
    implicit none
    integer :: it_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_f = 1 + mod((i - lo%llim_world)/lo%nindex/lo%naky, lo%ntheta0)
  end function it_idx_f

  pure function it_idx_lz (lo, i)
    implicit none
    integer :: it_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky, &
         lo%ntheta0)
  end function it_idx_lz

  pure function il_idx_g (lo, i)
    implicit none
    integer :: il_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nlambda)
  end function il_idx_g

  pure function ie_idx_g (lo, i)
    implicit none
    integer :: ie_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%nlambda, &
         lo%negrid)
  end function ie_idx_g

  pure function ie_idx_gint (lo, i)
    implicit none
    integer :: ie_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%negrid)
  end function ie_idx_gint

  pure function ie_idx_lz (lo, i)
    implicit none
    integer :: ie_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_lz = 1 + mod((i - lo%llim_world) &
         /(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0, lo%negrid)
  end function ie_idx_lz

  pure function is_idx_g (lo, i)
    implicit none
    integer :: is_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_g = 1 + mod((i - lo%llim_world) &
         /lo%naky/lo%ntheta0/lo%nlambda/lo%negrid, lo%nspec)
  end function is_idx_g

  pure function is_idx_gint (lo, i)
    implicit none
    integer :: is_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%negrid, &
         lo%nspec)
  end function is_idx_gint

  pure function is_idx_geint (lo, i)
    implicit none
    integer :: is_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nspec)
  end function is_idx_geint

  pure function is_idx_lz (lo, i)
    implicit none
    integer :: is_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_lz = 1 + mod((i - lo%llim_world) &
         /(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
  end function is_idx_lz

  pure function proc_id_g (lo, i)
    implicit none
    integer :: proc_id_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_g = i/lo%blocksize
  end function proc_id_g

  pure function proc_id_gint (lo, i)
    implicit none
    integer :: proc_id_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_gint = i/lo%blocksize
  end function proc_id_gint

  pure function proc_id_geint (lo, i)
    implicit none
    integer :: proc_id_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_geint = i/lo%blocksize
  end function proc_id_geint

  pure function proc_id_f (lo, i)
    implicit none
    integer :: proc_id_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_f = i/lo%blocksize
  end function proc_id_f

  pure function proc_id_lz (lo, i)
    implicit none
    integer :: proc_id_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_lz = i/lo%blocksize
  end function proc_id_lz

  pure function idx_g (lo, ik, it, il, ie, is)
    implicit none
    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    idx_g = ik-1 + lo%naky*(it-1 + lo%ntheta0*(il-1 + lo%nlambda*(ie-1 &
         + lo%negrid*(is-1))))
  end function idx_g

  pure function idx_gint (lo, ik, it, ie, is)
    implicit none
    integer :: idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, ie, is
    idx_gint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(ie-1 + lo%negrid*(is-1)))
  end function idx_gint

  pure function idx_geint (lo, ik, it, is)
    implicit none
    integer :: idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, is
    idx_geint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(is-1))
  end function idx_geint

  pure function idx_f (lo, if, ik, it)
    implicit none
    integer :: idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: if, ik, it
    idx_f = if-1 + lo%nindex*(ik-1 + lo%naky*(it-1))
  end function idx_f

  pure function idx_lz (lo, ig, ik, it, ie, is)
    implicit none
    integer :: idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is
    idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 + lo%naky*(it-1 &
         + lo%ntheta0*(ie-1 + lo%negrid*(is-1))))
  end function idx_lz

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

  pure function idx_local_f (lo, if, ik, it)
    implicit none
    logical :: idx_local_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: if, ik, it

    idx_local_f = idx_local(lo, idx(lo, if, ik, it))
  end function idx_local_f

  pure function ig_local_f (lo, ig)
    implicit none
    logical :: ig_local_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_f = lo%iproc == proc_id(lo, ig)
  end function ig_local_f

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

  pure subroutine gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, &
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
       if (il > je) then
          il = 2*je + 1
       else if (isign == 2) then
          il = 2*je - il
       end if
    end if
    ilz = idx(lz_lo, ig, ik_idx(g_lo,iglo), it_idx(g_lo,iglo), &
         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2lzidx

  pure subroutine lzidx2gidx (il, ilz, lz_lo, g_lo, jend, ig, isign, iglo)
    implicit none
    integer, intent (in) :: il, ilz
    type (lz_layout_type), intent (in) :: lz_lo
    type (g_layout_type), intent (in) :: g_lo
    integer, dimension (:), intent (in) :: jend
    integer, intent (out) :: ig, isign, iglo
    integer :: ilfold, je

    ig = ig_idx(lz_lo,ilz)
    je = jend(ig)

    if (je == 0) then
       je = lz_lo%ng2
       if (il <= je) then
          isign = 1
          ilfold = il
       else if (il > je+1) then
          isign = 2
          ilfold = 2*je + 1 - il
       else
          isign = 999999
          ilfold = 999999
       end if
    else
       if (il <= je) then
          isign = 1
          ilfold = il
       else if (il <= 2*je) then
          isign = 2
          ilfold = 2*je - il
       else
          isign = 999999
          ilfold = 999999
       end if
    end if

    iglo = idx(g_lo, ik_idx(lz_lo,ilz), it_idx(lz_lo,ilz), &
         ilfold, ie_idx(lz_lo,ilz), is_idx(lz_lo,ilz))
  end subroutine lzidx2gidx

end module gs2_layouts
