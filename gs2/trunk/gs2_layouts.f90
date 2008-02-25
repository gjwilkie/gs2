module gs2_layouts
  implicit none
  private

  public :: pe_layout, layout, local_field_solve
  public :: is_kx_local
  
  public :: init_dist_fn_layouts, init_gs2_layouts
  public :: g_lo, g_layout_type
  public :: gint_lo, gint_layout_type
  public :: geint_lo, geint_layout_type

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_jfields_layouts
  public :: jf_lo, jf_layout_type
  public :: mj, ij, dj

  public :: init_lorentz_layouts
  public :: lz_lo, lz_layout_type
  public :: gidx2lzidx
  public :: gidx2gintidx, gintidx2geidx

  public :: init_ediffuse_layouts
  public :: e_lo, e_layout_type

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

  logical :: local_field_solve, accel_lxyes, lambda_local
  character (len=5) :: layout

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

  type :: lz_layout_type
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
  end type lz_layout_type

  type (lz_layout_type) :: lz_lo

  type :: e_layout_type
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nlambda, nspec, nsign
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type e_layout_type

  type (e_layout_type) :: e_lo

  type :: xxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ntheta0, nx, nadd, negrid, nlambda, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
  end type xxf_layout_type

  type (xxf_layout_type) :: xxf_lo

  type :: yxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ny, ntheta0, nx, negrid, nlambda, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
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
     module procedure ig_idx_e
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
     module procedure ik_idx_e
     module procedure ik_idx_xxf
     module procedure ik_idx_accel
     module procedure ik_idx_accelx
  end interface

  interface it_idx
     module procedure it_idx_g
     module procedure it_idx_gint
     module procedure it_idx_geint
     module procedure it_idx_jf
     module procedure it_idx_lz
     module procedure it_idx_e
     module procedure it_idx_yxf
     module procedure it_idx_accel
     module procedure it_idx_accelx
  end interface

  interface il_idx
     module procedure il_idx_g
     module procedure il_idx_e
     module procedure il_idx_xxf
     module procedure il_idx_yxf
     module procedure il_idx_accelx
  end interface

  interface ie_idx
     module procedure ie_idx_g
     module procedure ie_idx_gint
     module procedure ie_idx_lz
     module procedure ie_idx_xxf
     module procedure ie_idx_yxf
     module procedure ie_idx_accelx
  end interface

  interface is_idx
     module procedure is_idx_g
     module procedure is_idx_gint
     module procedure is_idx_geint
     module procedure is_idx_lz
     module procedure is_idx_e
     module procedure is_idx_xxf
     module procedure is_idx_yxf
     module procedure is_idx_accelx
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
     module procedure proc_id_e
     module procedure proc_id_xxf
     module procedure proc_id_yxf
     module procedure proc_id_accelx
  end interface

  interface idx
     module procedure idx_g
     module procedure idx_gint
     module procedure idx_geint
     module procedure idx_f
     module procedure idx_jf
     module procedure idx_lz
     module procedure idx_e
     module procedure idx_xxf
     module procedure idx_yxf
     module procedure idx_accelx
  end interface

  interface idx_local
     module procedure idx_local_g,      ig_local_g
     module procedure idx_local_gint,   ig_local_gint
     module procedure idx_local_geint,  ig_local_geint
     module procedure idx_local_f,      ig_local_f
     module procedure idx_local_jf,     ig_local_jf
     module procedure idx_local_lz,     ig_local_lz
     module procedure idx_local_e,      ig_local_e
     module procedure idx_local_xxf,    ig_local_xxf
     module procedure idx_local_yxf,    ig_local_yxf
     module procedure idx_local_accelx, ig_local_accelx
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
    logical :: exist
    namelist /layouts_knobs/ layout, local_field_solve

    local_field_solve = .false.
    layout = 'lxyes'
    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)

  end subroutine read_parameters
    
  subroutine broadcast_results
    use mp, only: broadcast
    implicit none

    call broadcast (layout)
    call broadcast (local_field_solve)

  end subroutine broadcast_results

  subroutine check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)

    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky, ntheta0
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nxfacs, nlfacs, nblock
    integer :: i

    if (.not. layout == 'lxyes') then
       accel_lxyes = .false.
       return
    end if

    accel_lxyes = .true.

    allocate (facs(max(nspec,negrid,naky,ntheta0,nlambda)/2+1,5))
    call factors (nspec,   nsfacs, facs(:,1))
    call factors (negrid,  nefacs, facs(:,2))
    call factors (naky,    nyfacs, facs(:,3))
    call factors (ntheta0, nxfacs, facs(:,4))
    call factors (nlambda, nlfacs, facs(:,5))
    
!    do i=1,nsfacs-1
!       if (nproc == facs(i,1)) then
!          nblock = facs(i,1)
!          goto 100
!       end if
!    end do
!    
!    do i=1,nefacs-1
!       if (nproc == facs(i,2)*nspec) then
!          nblock = facs(i,2)*nspec
!          goto 100
!       end if
!    end do
    
    nblock = negrid*nspec
    do i=1,nyfacs-1
       if (nproc == facs(i,3)*negrid*nspec) goto 100
    end do
    
    do i=1,nxfacs-1
       if (nproc == facs(i,4)*naky*negrid*nspec) goto 100
    end do
    
    do i=1,nlfacs-1
       if (nproc == facs(i,5)*ntheta0*naky*negrid*nspec) goto 100
    end do
    
    accel_lxyes = .false.

100 deallocate (facs)

  end subroutine check_accel

  subroutine check_llocal (ntheta0, naky, nlambda, negrid, nspec)

    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky, ntheta0
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nxfacs, nlfacs
    integer :: i

    lambda_local = .true.

    allocate (facs(max(nspec,negrid,naky,ntheta0,nlambda)/2+1,5))

    select case (layout)

    case ('lxyes')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (naky,    nyfacs, facs(:,3))
       call factors (ntheta0, nxfacs, facs(:,4))
       call factors (nlambda, nlfacs, facs(:,5))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nyfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
       do i=1,nxfacs-1
          if (nproc == facs(i,4)*naky*negrid*nspec) goto 100
       end do
          
       do i=1,nlfacs-1
          if (nproc == facs(i,5)*ntheta0*naky*negrid*nspec) goto 100
       end do
          
    case ('lyxes')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (ntheta0, nxfacs, facs(:,3))
       call factors (naky,    nyfacs, facs(:,4))
       call factors (nlambda, nlfacs, facs(:,5))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nxfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
       do i=1,nyfacs-1
          if (nproc == facs(i,4)*ntheta0*negrid*nspec) goto 100
       end do
          
       do i=1,nlfacs-1
          if (nproc == facs(i,5)*naky*ntheta0*negrid*nspec) goto 100
       end do
          
    case ('lexys')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (naky,    nyfacs, facs(:,2))
       call factors (ntheta0, nxfacs, facs(:,3))
       call factors (negrid,  nefacs, facs(:,4))
       call factors (nlambda, nlfacs, facs(:,5))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nyfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nxfacs-1
          if (nproc == facs(i,3)*naky*nspec) goto 100
       end do
          
       do i=1,nefacs-1
          if (nproc == facs(i,4)*ntheta0*naky*nspec) goto 100
       end do
          
       do i=1,nlfacs-1
          if (nproc == facs(i,5)*negrid*ntheta0*naky*nspec) goto 100
       end do

    end select

    lambda_local = .false.    
                       
100 deallocate (facs)


  end subroutine check_llocal

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
    g_lo%ulim_world = naky*ntheta0*negrid*nlambda*nspec - 1
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
    gint_lo%ulim_world = naky*ntheta0*negrid*nspec - 1
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
    geint_lo%ulim_world = naky*ntheta0*nspec - 1
    geint_lo%blocksize = geint_lo%ulim_world/nproc + 1
    geint_lo%llim_proc = geint_lo%blocksize*iproc
    geint_lo%ulim_proc &
         = min(geint_lo%ulim_world, geint_lo%llim_proc + geint_lo%blocksize -1)
    geint_lo%ulim_alloc = max(geint_lo%llim_proc, geint_lo%ulim_proc)
    
  end subroutine init_dist_fn_layouts

  subroutine is_kx_local(negrid, nspec, nlambda, naky, ntheta0, kx_local)  

    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky, ntheta0
    logical, intent (out) :: kx_local
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nxfacs, nlfacs
    integer :: i

    kx_local = .true.

    allocate (facs(max(nspec,negrid,naky,nlambda)/2+1,3))

    select case (layout)

    case ('yxels')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (nlambda, nlfacs, facs(:,2))
       call factors (negrid,  nefacs, facs(:,3))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nlfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nefacs-1
          if (nproc == facs(i,3)*nlambda*nspec) goto 100
       end do
          
    case ('yxles')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (nlambda, nlfacs, facs(:,3))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nlfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
    case ('lxyes')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (naky,    nyfacs, facs(:,3))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nyfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
    case ('lyxes')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

    case ('lexys')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (naky,    nyfacs, facs(:,2))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nyfacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

    end select

    kx_local = .false.    
                       
100 deallocate (facs)

  end subroutine is_kx_local

  elemental function is_idx_g (lo, i)
    implicit none
    integer :: is_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%negrid/lo%nlambda, lo%nspec)
    case ('yxles')
       is_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%nlambda/lo%negrid, lo%nspec)
    case ('lexys')
       is_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
    case ('lxyes')
       is_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%ntheta0/lo%naky/lo%negrid, lo%nspec)
    case ('lyxes')
       is_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    end select

  end function is_idx_g

  elemental function il_idx_g (lo, i)
    implicit none
    integer :: il_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       il_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%negrid, lo%nlambda)
    case ('yxles')
       il_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nlambda)
    case ('lexys')
       il_idx_g = 1 + mod(i - lo%llim_world, lo%nlambda)
    case ('lxyes')
       il_idx_g = 1 + mod(i - lo%llim_world, lo%nlambda)
    case ('lyxes')
       il_idx_g = 1 + mod(i - lo%llim_world, lo%nlambda)
    end select

  end function il_idx_g

  elemental function ie_idx_g (lo, i)
    implicit none
    integer :: ie_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%negrid)
    case ('yxles')
       ie_idx_g = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%nlambda, lo%negrid)
    case ('lexys')
       ie_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda, lo%negrid)
    case ('lxyes')
       ie_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%ntheta0/lo%naky, lo%negrid)
    case ('lyxes')
       ie_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%naky/lo%ntheta0, lo%negrid)
    end select

  end function ie_idx_g

  elemental function it_idx_g (lo, i)
    implicit none
    integer :: it_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
    select case (layout)
    case ('yxels')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('yxles')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('lexys')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%negrid, lo%ntheta0)
    case ('lxyes')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda, lo%ntheta0)
    case ('lyxes')
       it_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%naky, lo%ntheta0)
    end select

  end function it_idx_g

  elemental function ik_idx_g (lo, i)
    implicit none
    integer :: ik_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ik_idx_g = 1 + mod(i - lo%llim_world, lo%naky)
    case ('yxles')
       ik_idx_g = 1 + mod(i - lo%llim_world, lo%naky)
    case ('lexys')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%negrid/lo%ntheta0, lo%naky)
    case ('lxyes')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda/lo%ntheta0, lo%naky)
    case ('lyxes')
       ik_idx_g = 1 + mod((i - lo%llim_world)/lo%nlambda, lo%naky)
    end select

  end function ik_idx_g

  elemental function proc_id_g (lo, i)
    implicit none
    integer :: proc_id_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_g = i/lo%blocksize
  end function proc_id_g

  elemental function idx_g (lo, ik, it, il, ie, is)
    implicit none
    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    
    select case (layout)
    case ('yxels')
       idx_g = ik-1 + lo%naky*(it-1 + lo%ntheta0*(ie-1 + lo%negrid*(il-1 + lo%nlambda*(is-1))))
    case ('yxles')
       idx_g = ik-1 + lo%naky*(it-1 + lo%ntheta0*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1))))
    case ('lexys')
       idx_g = il-1 + lo%nlambda*(ie-1 + lo%negrid*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1))))
    case ('lxyes')
       idx_g = il-1 + lo%nlambda*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1))))
    case ('lyxes')
       idx_g = il-1 + lo%nlambda*(ik-1 + lo%naky*(it-1 + lo%ntheta0*(ie-1 + lo%negrid*(is-1))))
    end select
  end function idx_g

  elemental function idx_local_g (lo, ik, it, il, ie, is)
    implicit none
    logical :: idx_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is

    idx_local_g = idx_local(lo, idx(lo, ik, it, il, ie, is))
  end function idx_local_g

  elemental function ig_local_g (lo, ig)
    implicit none
    logical :: ig_local_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_g = lo%iproc == proc_id(lo, ig)
  end function ig_local_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Once-integrated distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  elemental function is_idx_gint (lo, i)
    implicit none
    integer :: is_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    case ('yxles')
       is_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    case ('lexys')
       is_idx_gint = 1 + mod((i - lo%llim_world)/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
    case ('lxyes')
       is_idx_gint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky/lo%negrid, lo%nspec)
    case ('lyxes')
       is_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    end select
  end function is_idx_gint

  elemental function ie_idx_gint (lo, i)
    implicit none
    integer :: ie_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%negrid)
    case ('yxles')
       ie_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%negrid)
    case ('lexys')
       ie_idx_gint = 1 + mod(i - lo%llim_world, lo%negrid)
    case ('lxyes')
       ie_idx_gint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%negrid)
    case ('lyxes')
       ie_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%negrid)
    end select

  end function ie_idx_gint

  elemental function it_idx_gint (lo, i)
    implicit none
    integer :: it_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       it_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('yxles')
       it_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('lexys')
       it_idx_gint = 1 + mod((i - lo%llim_world)/lo%negrid, lo%ntheta0)
    case ('lxyes')
       it_idx_gint = 1 + mod((i - lo%llim_world), lo%ntheta0)
    case ('lyxes')
       it_idx_gint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    end select
  end function it_idx_gint

  elemental function ik_idx_gint (lo, i)
    implicit none
    integer :: ik_idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ik_idx_gint = 1 + mod(i - lo%llim_world, lo%naky)
    case ('yxles')
       ik_idx_gint = 1 + mod(i - lo%llim_world, lo%naky)
    case ('lexys')
       ik_idx_gint = 1 + mod((i - lo%llim_world)/lo%negrid/lo%ntheta0, lo%naky)
    case ('lxyes')
       ik_idx_gint = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
    case ('lyxes')
       ik_idx_gint = 1 + mod((i - lo%llim_world), lo%naky)
    end select
  end function ik_idx_gint

  elemental function proc_id_gint (lo, i)
    implicit none
    integer :: proc_id_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_gint = i/lo%blocksize
  end function proc_id_gint

  elemental function idx_gint (lo, ik, it, ie, is)
    implicit none
    integer :: idx_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, ie, is

    select case (layout)
    case ('yxels')
       idx_gint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(ie-1 + lo%negrid*(is-1)))
    case ('yxles')
       idx_gint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(ie-1 + lo%negrid*(is-1)))
    case ('lexys')
       idx_gint = ie-1 + lo%negrid*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
    case ('lxyes')
       idx_gint = it-1 + lo%ntheta0*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1)))
    case ('lyxes')
       idx_gint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(ie-1 + lo%negrid*(is-1)))
    end select
  end function idx_gint

  elemental function idx_local_gint (lo, ik, it, ie, is)
    implicit none
    logical :: idx_local_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, ie, is

    idx_local_gint = idx_local(lo, idx(lo, ik, it, ie, is))
  end function idx_local_gint

  elemental function ig_local_gint (lo, ig)
    implicit none
    logical :: ig_local_gint
    type (gint_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_gint = lo%iproc == proc_id(lo, ig)
  end function ig_local_gint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Twice-integrated distribution function layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  elemental function is_idx_geint (lo, i)
    implicit none
    integer :: is_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nspec)
    case ('yxles')
       is_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nspec)
    case ('lexys')
       is_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%nspec)
    case ('lxyes')
       is_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%nspec)
    case ('lyxes')
       is_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky/lo%ntheta0, lo%nspec)
    end select

  end function is_idx_geint

  elemental function it_idx_geint (lo, i)
    implicit none
    integer :: it_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       it_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('yxles')
       it_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    case ('lexys')
       it_idx_geint = 1 + mod(i - lo%llim_world, lo%ntheta0)
    case ('lxyes')
       it_idx_geint = 1 + mod(i - lo%llim_world, lo%ntheta0)
    case ('lyxes')
       it_idx_geint = 1 + mod((i - lo%llim_world)/lo%naky, lo%ntheta0)
    end select
  end function it_idx_geint

  elemental function ik_idx_geint (lo, i)
    implicit none
    integer :: ik_idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ik_idx_geint = 1 + mod(i - lo%llim_world, lo%naky)
    case ('yxles')
       ik_idx_geint = 1 + mod(i - lo%llim_world, lo%naky)
    case ('lexys')
       ik_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
    case ('lxyes')
       ik_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
    case ('lyxes')
       ik_idx_geint = 1 + mod((i - lo%llim_world), lo%naky)
    end select
  end function ik_idx_geint

  elemental function proc_id_geint (lo, i)
    implicit none
    integer :: proc_id_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_geint = i/lo%blocksize
  end function proc_id_geint

  elemental function idx_geint (lo, ik, it, is)
    implicit none
    integer :: idx_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, is

    select case (layout)
    case ('yxels')
       idx_geint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(is-1))
    case ('yxles')
       idx_geint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(is-1))
    case ('lexys')
       idx_geint = it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1))
    case ('lxyes')
       idx_geint = it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1))
    case ('lyxes')
       idx_geint = ik-1 + lo%naky*(it-1 + lo%ntheta0*(is-1))
    end select
  end function idx_geint

  elemental function idx_local_geint (lo, ik, it, is)
    implicit none
    logical :: idx_local_geint
    type (geint_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, is

    idx_local_geint = idx_local(lo, idx(lo, ik, it, is))
  end function idx_local_geint

  elemental function ig_local_geint (lo, ig)
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
! Energy scattering layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_ediffuse_layouts &
       (ntgrid, naky, ntheta0, nlambda, nspec)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, nspec
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    
    e_lo%iproc = iproc
    e_lo%ntgrid = ntgrid
    e_lo%nsign = 2
    e_lo%ntheta0 = ntheta0
    e_lo%naky = naky
    e_lo%nspec = nspec
    e_lo%nlambda = nlambda
    e_lo%llim_world = 0
    e_lo%ulim_world = (2*ntgrid+1)*naky*ntheta0*nlambda*nspec*2 - 1
    e_lo%blocksize = e_lo%ulim_world/nproc + 1
    e_lo%llim_proc = e_lo%blocksize*iproc
    e_lo%ulim_proc &
         = min(e_lo%ulim_world, e_lo%llim_proc + e_lo%blocksize - 1)
    e_lo%ulim_alloc = max(e_lo%llim_proc, e_lo%ulim_proc)

  end subroutine init_ediffuse_layouts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lorentz layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_lorentz_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2
    integer :: ngroup, nprocset
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

    if (layout == 'lxyes' .or. layout == 'lexys' .or. layout == 'lyxes') &
         call check_llocal (ntheta0, naky, nlambda, negrid, nspec)

    if (lambda_local) then

       lz_lo%groupblocksize = (naky*ntheta0*negrid*nspec-1)/nproc + 1
       
       ngroup = min (nproc, naky*ntheta0*negrid*nspec)
       lz_lo%ngroup = ngroup
       
       nprocset = nproc / lz_lo%ngroup
       lz_lo%nprocset = nprocset

       lz_lo%iset   = mod (iproc, nprocset)
       lz_lo%igroup = mod (iproc/nprocset, ngroup)

       lz_lo%llim_group = 0
       lz_lo%ulim_group = (2*ntgrid+1)*lz_lo%groupblocksize - 1
       lz_lo%gsize      = (2*ntgrid+1)*lz_lo%groupblocksize 
       
       lz_lo%nset = lz_lo%ulim_group/lz_lo%nprocset + 1
       
       lz_lo%llim_proc = lz_lo%igroup*lz_lo%gsize + lz_lo%iset*lz_lo%nset
       lz_lo%ulim_proc = min(lz_lo%ulim_group+lz_lo%igroup*lz_lo%gsize, &
            lz_lo%llim_proc + lz_lo%nset - 1)
       lz_lo%ulim_alloc = max(lz_lo%llim_proc, lz_lo%ulim_proc)
    else

       lz_lo%blocksize = lz_lo%ulim_world/nproc + 1
       lz_lo%llim_proc = lz_lo%blocksize*iproc
       lz_lo%ulim_proc &
            = min(lz_lo%ulim_world, lz_lo%llim_proc + lz_lo%blocksize - 1)
       lz_lo%ulim_alloc = max(lz_lo%llim_proc, lz_lo%ulim_proc)
    end if

  end subroutine init_lorentz_layouts

  elemental function is_idx_lz (lo, i)
    implicit none
    integer :: is_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    case ('yxles')
       is_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    case ('lexys')
       is_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%negrid/lo%ntheta0/lo%naky, lo%nspec)
    case ('lxyes')
       is_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%ntheta0/lo%naky/lo%negrid, lo%nspec)
    case ('lyxes')
       is_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0/lo%negrid, lo%nspec)
    end select

  end function is_idx_lz

  elemental function ie_idx_lz (lo, i)
    implicit none
    integer :: ie_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0, lo%negrid)
    case ('yxles')
       ie_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0, lo%negrid)
    case ('lexys')
       ie_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%negrid)
    case ('lxyes')
       ie_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%ntheta0/lo%naky, lo%negrid)
    case ('lyxes')
       ie_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky/lo%ntheta0, lo%negrid)
    end select
  end function ie_idx_lz

  elemental function it_idx_lz (lo, i)
    implicit none
    integer :: it_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       it_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky, lo%ntheta0)
    case ('yxles')
       it_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky, lo%ntheta0)
    case ('lexys')
       it_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%negrid, lo%ntheta0)
    case ('lxyes')
       it_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%ntheta0)
    case ('lyxes')
       it_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%naky, lo%ntheta0)
    end select
  end function it_idx_lz

  elemental function ik_idx_lz (lo, i)
    implicit none
    integer :: ik_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ik_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%naky)
    case ('yxles')
       ik_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%naky)
    case ('lexys')
       ik_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%negrid/lo%ntheta0, lo%naky)
    case ('lxyes')
       ik_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%ntheta0, lo%naky)
    case ('lyxes')
       ik_idx_lz = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1), lo%naky)
    end select
  end function ik_idx_lz

  elemental function ig_idx_lz (lo, i)
    implicit none
    integer :: ig_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('yxles')
       ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('lexys')
       ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('lxyes')
       ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('lyxes')
       ig_idx_lz = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    end select
  end function ig_idx_lz

  elemental function idx_lz (lo, ig, ik, it, ie, is)
    implicit none
    integer :: idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is

    select case (layout)
    case ('yxels')
       idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 + lo%naky*(it-1 &
            + lo%ntheta0*(ie-1 + lo%negrid*(is-1))))
    case ('yxles')
       idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 + lo%naky*(it-1 &
            + lo%ntheta0*(ie-1 + lo%negrid*(is-1))))
    case ('lexys')
       idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ie-1 + lo%negrid*(it-1 &
            + lo%ntheta0*(ik-1 + lo%naky*(is-1))))
    case ('lxyes')
       idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(it-1 + lo%ntheta0*(ik-1 &
            + lo%naky*(ie-1 + lo%negrid*(is-1))))
    case ('lyxes')
       idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 + lo%naky*(it-1 &
            + lo%ntheta0*(ie-1 + lo%negrid*(is-1))))
    end select

  end function idx_lz

  elemental function proc_id_lz (lo, i)
    implicit none
    integer :: proc_id_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    if (lambda_local) then
       proc_id_lz = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
       proc_id_lz = i/lo%blocksize
    end if

  end function proc_id_lz

  elemental function idx_local_lz (lo, ig, ik, it, ie, is)
    implicit none
    logical :: idx_local_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is

    idx_local_lz = idx_local(lo, idx(lo, ig, ik, it, ie, is))
  end function idx_local_lz

  elemental function ig_local_lz (lo, ig)
    implicit none
    logical :: ig_local_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_lz = lo%iproc == proc_id(lo, ig)
  end function ig_local_lz

  elemental function is_idx_e (lo, i)
    implicit none
    integer :: is_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%naky/lo%ntheta0/lo%nlambda, lo%nspec)
    case ('yxles')
       is_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%naky/lo%ntheta0/lo%nlambda, lo%nspec)
    case ('lexys')
       is_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%ntheta0/lo%naky, lo%nspec)
    case ('lxyes')
       is_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%ntheta0/lo%naky, lo%nspec)
    case ('lyxes')
       is_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%naky/lo%ntheta0, lo%nspec)
    end select
  end function is_idx_e

  elemental function il_idx_e (lo, i)
    implicit none
    integer :: il_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       il_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%naky/lo%ntheta0, lo%nlambda)
    case ('yxles')
       il_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%naky/lo%ntheta0, lo%nlambda)
    case ('lexys')
       il_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lxyes')
       il_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lyxes')
       il_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    end select
  end function il_idx_e

  elemental function it_idx_e (lo, i)
    implicit none
    integer :: it_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       it_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%naky, lo%ntheta0)
    case ('yxles')
       it_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%naky, lo%ntheta0)
    case ('lexys')
       it_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%ntheta0)
    case ('lxyes')
       it_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%ntheta0)
    case ('lyxes')
       it_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%naky, lo%ntheta0)
    end select
  end function it_idx_e

  elemental function ik_idx_e (lo, i)
    implicit none
    integer :: ik_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ik_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign, lo%naky)
    case ('yxles')
       ik_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign, lo%naky)
    case ('lexys')
       ik_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%ntheta0, lo%naky)
    case ('lxyes')
       ik_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%ntheta0, lo%naky)
    case ('lyxes')
       ik_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%naky)
    end select
  end function ik_idx_e

  elemental function ig_idx_e (lo, i)
    implicit none
    integer :: ig_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ig_idx_e = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('yxles')
       ig_idx_e = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('lexys')
       ig_idx_e = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('lxyes')
       ig_idx_e = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    case ('lyxes')
       ig_idx_e = -lo%ntgrid + mod(i - lo%llim_world, 2*lo%ntgrid + 1)
    end select
  end function ig_idx_e

  elemental function idx_e (lo, ig, isign, ik, it, il, is)
    implicit none
    integer :: idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, it, il, is

    select case (layout)
    case ('yxels')
       idx_e = ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 + lo%nsign*(ik-1 &
            + lo%naky*(it-1 + lo%ntheta0*(il-1 + lo%nlambda*(is-1)))))
    case ('yxles')
       idx_e = ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 + lo%nsign*(ik-1 &
            + lo%naky*(it-1 + lo%ntheta0*(il-1 + lo%nlambda*(is-1)))))
    case ('lexys')
       idx_e = ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 + lo%nsign*(il-1 &
            + lo%nlambda*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1)))))
    case ('lxyes')
       idx_e = ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 + lo%nsign*(il-1 &
            + lo%nlambda*(it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1)))))
    case ('lyxes')
       idx_e = ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 + lo%nsign*(il-1 &
            + lo%nlambda*(ik-1 + lo%naky*(it-1 + lo%ntheta0*(is-1)))))
    end select

  end function idx_e

  elemental function proc_id_e (lo, i)
    implicit none
    integer :: proc_id_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_e = i/lo%blocksize

  end function proc_id_e

  elemental function idx_local_e (lo, ig, isign, ik, it, il, is)
    implicit none
    logical :: idx_local_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, it, il, is

    idx_local_e = idx_local(lo, idx(lo, ig, isign, ik, it, il, is))
  end function idx_local_e

  elemental function ig_local_e (lo, ig)
    implicit none
    logical :: ig_local_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_e = lo%iproc == proc_id(lo, ig)
  end function ig_local_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_x_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    logical, save :: initialized = .false.
    integer :: nprocset, ngroup, ip, nblock

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

    call check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)
    if (accel_lxyes) then  

       xxf_lo%groupblocksize = (nblock-1)/nproc + 1

       ngroup = min (nproc, nblock)
       xxf_lo%ngroup = ngroup
       
       nprocset = nproc / xxf_lo%ngroup
       xxf_lo%nprocset = nprocset

       xxf_lo%iset   = mod (iproc, nprocset)
       xxf_lo%igroup = mod (iproc/nprocset, ngroup)

       xxf_lo%llim_group = 0
       xxf_lo%ulim_group = naky*(2*ntgrid+1)*2*nlambda*xxf_lo%groupblocksize - 1
       xxf_lo%gsize      = naky*(2*ntgrid+1)*2*nlambda*xxf_lo%groupblocksize 
       
       xxf_lo%nset = xxf_lo%ulim_group/xxf_lo%nprocset + 1
       
       xxf_lo%llim_proc = xxf_lo%igroup*xxf_lo%gsize + xxf_lo%iset*xxf_lo%nset
       xxf_lo%ulim_proc = min(xxf_lo%ulim_group+xxf_lo%igroup*xxf_lo%gsize, &
            xxf_lo%llim_proc + xxf_lo%nset - 1)
       xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)

    else
       xxf_lo%blocksize = xxf_lo%ulim_world/nproc + 1
       xxf_lo%llim_proc = xxf_lo%blocksize*iproc
       xxf_lo%ulim_proc &
            = min(xxf_lo%ulim_world, xxf_lo%llim_proc + xxf_lo%blocksize - 1)
       xxf_lo%ulim_alloc = max(xxf_lo%llim_proc, xxf_lo%ulim_proc)
    end if

!    call barrier
!    do ip=0,nproc-1
!       if (ip == iproc) then
!          write (*,*) 'iproc= ',ip,' llim= ',xxf_lo%llim_proc,' ulim= ',xxf_lo%ulim_proc, &
!               & ' iset= ',xxf_lo%iset,' igroup= ',xxf_lo%igroup
!       end if
!    call barrier
!    end do

  end subroutine init_x_transform_layouts

  elemental function is_idx_xxf (lo, i)
    implicit none
    integer :: is_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%negrid/lo%nlambda, lo%nspec)
    case ('yxles')
       is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    case ('lexys')
       is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    case ('lxyes')
       is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    case ('lyxes')
       is_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    end select

  end function is_idx_xxf

  elemental function ie_idx_xxf (lo, i)
    implicit none
    integer :: ie_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign, lo%negrid)
    case ('yxles')
       ie_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    case ('lexys')
       ie_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    case ('lxyes')
       ie_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    case ('lyxes')
       ie_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    end select
  end function ie_idx_xxf

  elemental function il_idx_xxf (lo, i)
    implicit none
    integer :: il_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       il_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign/lo%negrid, lo%nlambda)
    case ('yxles')
       il_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lexys')
       il_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lxyes')
       il_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lyxes')
       il_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    end select
  end function il_idx_xxf

  elemental function isign_idx_xxf (lo, i)
    implicit none
    integer :: isign_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       isign_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), lo%nsign)
    case ('yxles')
       isign_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), lo%nsign)
    case ('lexys')
       isign_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), lo%nsign)
    case ('lxyes')
       isign_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), lo%nsign)
    case ('lyxes')
       isign_idx_xxf = 1 + mod((i - lo%llim_world)/lo%naky/(2*lo%ntgrid + 1), lo%nsign)
    end select
  end function isign_idx_xxf

  elemental function ig_idx_xxf (lo, i)
    implicit none
    integer :: ig_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
    case ('yxles')
       ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
    case ('lexys')
       ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
    case ('lxyes')
       ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
    case ('lyxes')
       ig_idx_xxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%naky, 2*lo%ntgrid + 1)
    end select
  end function ig_idx_xxf

  elemental function ik_idx_xxf (lo, i)
    implicit none
    integer :: ik_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
    case ('yxles')
       ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
    case ('lexys')
       ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
    case ('lxyes')
       ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
    case ('lyxes')
       ik_idx_xxf = 1 + mod(i - lo%llim_world, lo%naky)
    end select
  end function ik_idx_xxf

  elemental function idx_xxf (lo, ig, isign, ik, il, ie, is)
    implicit none
    integer :: idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, il, ie, is

    select case (layout)
    case ('yxels')
       idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(ie-1 + lo%negrid*(il-1 + lo%nlambda*(is-1)))))
    case ('yxles')
       idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    case ('lexys')
       idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    case ('lxyes')
       idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    case ('lyxes')
       idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    end select
  end function idx_xxf

  elemental function proc_id_xxf (lo, i)
    implicit none
    integer :: proc_id_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    if (accel_lxyes) then
       proc_id_xxf = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
       proc_id_xxf = i/lo%blocksize
    end if

  end function proc_id_xxf

  elemental function idx_local_xxf (lo, ig, isign, ik, il, ie, is)
    implicit none
    logical :: idx_local_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, il, ie, is
    idx_local_xxf = idx_local (lo, idx(lo, ig, isign, ik, il, ie, is))
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
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    logical, save :: initialized = .false.
    integer :: nnx, nny, ngroup, nprocset, nblock

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

    call check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)

    if (accel_lxyes) then  

       yxf_lo%groupblocksize = (nblock-1)/nproc + 1

       ngroup = min (nproc, nblock)
       yxf_lo%ngroup = ngroup
       
       nprocset = nproc / yxf_lo%ngroup
       yxf_lo%nprocset = nprocset

       yxf_lo%iset   = mod (iproc, nprocset)
       yxf_lo%igroup = mod (iproc/nprocset, ngroup)

       yxf_lo%llim_group = 0
       yxf_lo%ulim_group = nnx*(2*ntgrid+1)*2*nlambda - 1
       yxf_lo%gsize      = nnx*(2*ntgrid+1)*2*nlambda*yxf_lo%groupblocksize
       
       yxf_lo%nset = yxf_lo%ulim_group/yxf_lo%nprocset + 1
       
       yxf_lo%llim_proc = yxf_lo%igroup*yxf_lo%gsize + yxf_lo%iset*yxf_lo%nset
       yxf_lo%ulim_proc = min(yxf_lo%ulim_group+yxf_lo%igroup*yxf_lo%gsize, &
            yxf_lo%llim_proc + yxf_lo%nset - 1)
       yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)

    else
       yxf_lo%blocksize = yxf_lo%ulim_world/nproc + 1
       yxf_lo%llim_proc = yxf_lo%blocksize*iproc
       yxf_lo%ulim_proc &
            = min(yxf_lo%ulim_world, yxf_lo%llim_proc + yxf_lo%blocksize - 1)
       yxf_lo%ulim_alloc = max(yxf_lo%llim_proc, yxf_lo%ulim_proc)
    end if

  end subroutine init_y_transform_layouts

  elemental function is_idx_yxf (lo, i)
    implicit none
    integer :: is_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%negrid/lo%nlambda, lo%nspec)
    case ('yxles')
       is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    case ('lexys')
       is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    case ('lxyes')
       is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    case ('lyxes')
       is_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda/lo%negrid, lo%nspec)
    end select
  end function is_idx_yxf

  elemental function ie_idx_yxf (lo, i)
    implicit none
    integer :: ie_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign, lo%negrid)
    case ('yxles')
       ie_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    case ('lexys')
       ie_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    case ('lxyes')
       ie_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    case ('lyxes')
       ie_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%nlambda, lo%negrid)
    end select
  end function ie_idx_yxf

  elemental function il_idx_yxf (lo, i)
    implicit none
    integer :: il_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       il_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign/lo%negrid, lo%nlambda)
    case ('yxles')
       il_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lexys')
       il_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lxyes')
       il_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    case ('lyxes')
       il_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid + 1)/lo%nsign, lo%nlambda)
    end select
  end function il_idx_yxf

  elemental function isign_idx_yxf (lo, i)
    implicit none
    integer :: isign_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       isign_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), lo%nsign)
    case ('yxles')
       isign_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), lo%nsign)
    case ('lexys')
       isign_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), lo%nsign)
    case ('lxyes')
       isign_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), lo%nsign)
    case ('lyxes')
       isign_idx_yxf = 1 + mod((i - lo%llim_world)/lo%nx/(2*lo%ntgrid+1), lo%nsign)
    end select
  end function isign_idx_yxf

  elemental function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
    case ('yxles')
       ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
    case ('lexys')
       ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
    case ('lxyes')
       ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
    case ('lyxes')
       ig_idx_yxf = -lo%ntgrid + mod((i - lo%llim_world)/lo%nx, 2*lo%ntgrid + 1)
    end select
  end function ig_idx_yxf

  elemental function it_idx_yxf (lo, i)
    implicit none
    integer :: it_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
    case ('yxles')
       it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
    case ('lexys')
       it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
    case ('lxyes')
       it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
    case ('lyxes')
       it_idx_yxf = 1 + mod(i - lo%llim_world, lo%nx)
    end select
  end function it_idx_yxf

  elemental function idx_accelx (lo, ik, it, il, ie, is)
    implicit none
    integer :: idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is

    select case (layout)
    case ('yxels')
       idx_accelx = ik-1 + lo%ny*(it-1 + lo%nx*(ie-1 + &
            lo%negrid*(il-1 + lo%nlambda*(is-1))))
    case ('yxles')
       idx_accelx = ik-1 + lo%ny*(it-1 + lo%nx*(il-1 + &
            lo%nlambda*(ie-1 + lo%negrid*(is-1))))
    end select
  end function idx_accelx

  elemental function idx_yxf (lo, ig, isign, it, il, ie, is)
    implicit none
    integer :: idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, it, il, ie, is

    select case (layout)
    case ('yxels')
       idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(ie-1 + lo%negrid*(il-1 + lo%nlambda*(is-1)))))
    case ('yxles')
       idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    case ('lexys')
       idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    case ('lxyes')
       idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    case ('lyxes')
       idx_yxf = it-1 + lo%nx*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    end select
  end function idx_yxf

  elemental function proc_id_accelx (lo, i)
    implicit none
    integer :: proc_id_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_accelx = i/lo%blocksize

  end function proc_id_accelx

  elemental function proc_id_yxf (lo, i)
    implicit none
    integer :: proc_id_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    if (accel_lxyes) then
       proc_id_yxf = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
       proc_id_yxf = i/lo%blocksize
    end if
  end function proc_id_yxf

  elemental function idx_local_accelx (lo, ik, it, il, ie, is)
    implicit none
    logical :: idx_local_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
    idx_local_accelx = idx_local (lo, idx(lo, ik, it, il, ie, is))
  end function idx_local_accelx

  elemental function ig_local_accelx (lo, i)
    implicit none
    logical ig_local_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_local_accelx = lo%iproc == proc_id(lo, i)
  end function ig_local_accelx

  elemental function idx_local_yxf (lo, ig, isign, it, il, ie, is)
    implicit none
    logical :: idx_local_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, it, il, ie, is
    idx_local_yxf = idx_local (lo, idx(lo, ig, isign, it, il, ie, is))
  end function idx_local_yxf

  elemental function ig_local_yxf (lo, i)
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
    accel_lo%nia = negrid*nlambda*nspec/nproc
    accel_lo%llim_world = 0
    accel_lo%ulim_world = nnx*(nny/2+1)*nlambda*negrid*nspec - 1
    accel_lo%blocksize = accel_lo%ulim_world/nproc + 1
    accel_lo%llim_proc = accel_lo%blocksize*iproc
    accel_lo%ulim_proc &
         = min(accel_lo%ulim_world, accel_lo%llim_proc + accel_lo%blocksize - 1)
    accel_lo%ulim_alloc = max(accel_lo%llim_proc, accel_lo%ulim_proc)

  end subroutine init_accel_transform_layouts

  elemental function is_idx_accelx (lo, i)
    implicit none
    integer :: is_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%negrid/lo%nlambda, lo%nspec)
  end function is_idx_accelx

  elemental function ie_idx_accelx (lo, i)
    implicit none
    integer :: ie_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       ie_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny, lo%negrid)
    case ('yxles')
       ie_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%nlambda, lo%negrid)
    end select
  end function ie_idx_accelx

  elemental function il_idx_accelx (lo, i)
    implicit none
    integer :: il_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    select case (layout)
    case ('yxels')
       il_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%negrid, lo%nlambda)
    case ('yxles')
       il_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny, lo%nlambda)
    end select
  end function il_idx_accelx

  elemental function it_idx_accelx (lo, i)
    implicit none
    integer :: it_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
  end function it_idx_accelx

  elemental function ik_idx_accelx (lo, i)
    implicit none
    integer :: ik_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
  end function ik_idx_accelx

  elemental function it_idx_accel (lo, i)
    implicit none
    integer :: it_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
  end function it_idx_accel

  elemental function ik_idx_accel (lo, i)
    implicit none
    integer :: ik_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
  end function ik_idx_accel

  elemental function dealiasing (lo, ia)
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

  elemental subroutine gidx2gintidx (g_lo, iglo, gint_lo, il, igint)
    implicit none
    type (g_layout_type), intent (in) :: g_lo
    integer, intent (in) :: iglo
    type (gint_layout_type), intent (in) :: gint_lo
    integer, intent (out) :: il, igint

    il = il_idx(g_lo, iglo)

    igint = idx(gint_lo, ik_idx(g_lo,iglo), it_idx(g_lo,iglo), &
         ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
  end subroutine gidx2gintidx

  elemental subroutine gintidx2geidx (gint_lo, igint, ie, geint_lo, igeint)
    implicit none
    type (gint_layout_type), intent (in) :: gint_lo
    integer, intent (in) :: igint
    type (geint_layout_type), intent (in) :: geint_lo
    integer, intent (out) :: ie, igeint

    ie = ie_idx(gint_lo, igint)

    igeint = idx(geint_lo, ik_idx(gint_lo, igint), it_idx(gint_lo, igint), &
         is_idx(gint_lo, igint))

  end subroutine gintidx2geidx

  elemental subroutine gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
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

  elemental subroutine xxfidx2gidx (it, ixxf, xxf_lo, g_lo, ig, isign, iglo)
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

  elemental subroutine xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
    implicit none
    integer, intent (in) :: it, ixxf
    type (xxf_layout_type), intent (in) :: xxf_lo
    type (yxf_layout_type), intent (in) :: yxf_lo
    integer, intent (out) :: ik, iyxf

    ik = ik_idx(xxf_lo,ixxf)
    iyxf = idx(yxf_lo, ig_idx(xxf_lo,ixxf), isign_idx(xxf_lo,ixxf), &
         it, il_idx(xxf_lo,ixxf), ie_idx(xxf_lo,ixxf), is_idx(xxf_lo,ixxf))
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
    ixxf = idx(xxf_lo, ig_idx(yxf_lo,iyxf), isign_idx(yxf_lo,iyxf), &
         ik0, il_idx(yxf_lo,iyxf), ie_idx(yxf_lo,iyxf), is_idx(yxf_lo,iyxf))
  end subroutine yxfidx2xxfidx

  subroutine pe_layout (char)

    character (1), intent (out) :: char

    select case (layout)
    case ('yxels')
       char = 'v'
    case ('yxles')
       char = 'v'
    case ('lexys')
       char = 'x'
    case ('lxyes')
       char = 'm'    ! mixed
    case ('lyxes')
       char = 'b'    ! big, since this layout makes sense mainly for big runs?
    end select

  end subroutine pe_layout

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

end module gs2_layouts



! save:

!       xxf_lo%groupblocksize = (negrid*nspec-1)/nproc + 1
!
!       xxf_lo%ngroup = min (nproc, negrid*nspec)
!       xxf_lo%nprocset = nproc / xxf_lo%ngroup
!
!       xxf_lo%llim_group = xxf_lo%groupblocksize*(iproc/nprocset)
!       xxf_lo%ulim_group = xxf_lo%llim_group + xxf_lo%groupblocksize - 1
!       
!       xxf_lo%nset = (naky*(2*ntgrid+1)*2*nlambda - 1)/nprocset + 1
!       
!       xxf_lo%llim_set = xxf_lo%nset*mod(iproc, nprocset)
!       xxf_lo%ulim_set = min( &
!            naky*(2*ntgrid+1)*2*nlambda - 1, xxf_lo%llim_set + xxf_lo%nset - 1)
!
!       llim_proc:
!       ulim_proc:
!       ulim_alloc:
!
!
!       ixxf = iset + xxf_lo%nset*igroup + (iproc/xxf_lo%ngroup)*xxf_lo%nset*xxf_lo%ngroup

! igroup: llim_group:ulim_group
! iset: llim_set:ulim_set
!
