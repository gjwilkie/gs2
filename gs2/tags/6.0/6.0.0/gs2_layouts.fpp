! Modifications for unbalanced decomposition functionality:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2012
! on behalf of EPSRC for the HECToR project
# include "define.inc"

module gs2_layouts

! TT: What are gint_layout_type and geint_layout_type?

! TT>
  use layouts_type, only: g_layout_type, lz_layout_type, e_layout_type
  use layouts_type, only: le_layout_type, p_layout_type
! <TT

  implicit none

  private

  public :: factors
  public :: pe_layout, layout, local_field_solve
  public :: is_kx_local
  
  public :: init_dist_fn_layouts, init_gs2_layouts
  public :: wnml_gs2_layouts
  public :: init_parity_layouts ! MAB
! TT>
  public :: g_lo
! <TT
  public :: gint_lo, gint_layout_type
  public :: geint_lo, geint_layout_type
  public :: p_lo

! HJL>
  public :: finish_layouts
! <HJL
  public :: finish_fields_layouts, finish_jfields_layouts

  public :: init_fields_layouts
  public :: f_lo, f_layout_type

  public :: init_jfields_layouts
  public :: jf_lo, jf_layout_type
  public :: mj, ij, dj

  public :: init_lambda_layouts
! TT>
  public :: lz_lo
! <TT
  public :: gidx2lzidx
  public :: gidx2gintidx, gintidx2geidx

  public :: init_energy_layouts
! TT>
  public :: e_lo
! <TT

! MAB>
! ported le_lo from agk
  public :: init_le_layouts
  public :: le_lo
! <MAB

  public :: init_x_transform_layouts, init_y_transform_layouts
  public :: calculate_unbalanced_x, calculate_unbalanced_y, calculate_idle_processes
  public :: xxf_lo, xxf_layout_type, yxf_lo, yxf_layout_type
  public :: gidx2xxfidx, xxfidx2yxfidx, yxfidx2xxfidx, xxfidx2gidx
  public :: xxf_ky_is_zero

  public :: init_accel_transform_layouts
  public :: accel_lo, accel_layout_type, dealiasing
  public :: accelx_lo, accelx_layout_type

  public :: ig_idx, ik_idx, it_idx, il_idx, ie_idx, is_idx, if_idx, isign_idx
  public :: im_idx, in_idx, ij_idx, ifield_idx
  public :: idx, proc_id, idx_local

  public :: opt_local_copy
  public :: opt_redist_init !Do we use optimised redist init routines?
  public :: intmom_sub !Do we use sub-communicators in velocity space?
  public :: intspec_sub !Do we use sub-communicators in integrate_species?

  public :: fft_use_wisdom, fft_wisdom_file, fft_measure_plan


  logical :: initialized_x_transform = .false.
  logical :: initialized_y_transform = .false.
  logical :: initialized_layouts = .false.
  logical :: initialized_dist_fn_layouts = .false.
  logical :: initialized_fields_layouts = .false.
  logical :: initialized_jfields_layouts = .false.
  logical :: initialized_le_layouts = .false.
  logical :: initialized_energy_layouts = .false.
  logical :: initialized_lambda_layouts = .false.
  logical :: initialized_parity_layouts = .false.
  logical :: initialized_accel_transform_layouts = .false.

  logical :: opt_local_copy
  logical :: opt_redist_init
  logical :: intmom_sub
  logical :: intspec_sub
  logical :: local_field_solve, accel_lxyes, lambda_local, unbalanced_xxf, unbalanced_yxf
  real :: max_unbalanced_xxf, max_unbalanced_yxf
  character (len=5) :: layout
  character (len=1000) :: fft_wisdom_file
  logical :: fft_use_wisdom, fft_measure_plan
  logical :: exist

! TT>
!!$  type :: g_layout_type
!!$     integer :: iproc
!!$     integer :: naky, ntheta0, nlambda, negrid, nspec
!!$     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!!$  end type g_layout_type
! <TT

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
! TT>
  type (le_layout_type) :: le_lo  ! new type
! <TT
  type (p_layout_type) :: p_lo

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
     integer :: ig_min,ig_max,ifield_min,ifield_max,it_min,it_max,ik_min,ik_max
  end type jf_layout_type

  type (jf_layout_type) :: jf_lo

  integer, dimension(:), allocatable :: ij, mj
  integer, dimension(:,:), allocatable :: dj

! TT>
!!$  type :: lz_layout_type
!!$     integer :: iproc
!!$     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2
!!$     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
!!$     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
!!$  end type lz_layout_type
! <TT

  type (lz_layout_type) :: lz_lo

! TT>
!!$  type :: e_layout_type
!!$     integer :: iproc
!!$     integer :: ntgrid, naky, ntheta0, nlambda, nspec, nsign
!!$     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!!$  end type e_layout_type
! <TT

  type (e_layout_type) :: e_lo

  type :: xxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ntheta0, nx, nadd, negrid, nlambda, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
     !<DD>
     integer :: ig_ord, isgn_ord, ik_ord, il_ord, ie_ord, is_ord
     integer :: ig_comp, isgn_comp, ik_comp, il_comp, ie_comp, is_comp
     integer,dimension(6) :: compound_count
  end type xxf_layout_type

  type (xxf_layout_type) :: xxf_lo

  type :: yxf_layout_type
     integer :: iproc
     integer :: ntgrid, nsign, naky, ny, ntheta0, nx, negrid, nlambda, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
     integer :: small_block_size, block_multiple, large_block_size, num_small, num_large
     integer :: small_block_balance_factor, large_block_balance_factor
     !<DD>
     integer :: ig_ord, isgn_ord, it_ord, il_ord, ie_ord, is_ord
     integer :: ig_comp, isgn_comp, it_comp, il_comp, ie_comp, is_comp
     integer,dimension(6) :: compound_count
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
     module procedure ifield_idx_jf !<DD> Added
  end interface

  interface if_idx
     module procedure if_idx_f
     module procedure if_idx_jf
  end interface

  interface ig_idx
     module procedure ig_idx_lz
     module procedure ig_idx_e
! TT>
     module procedure ig_idx_le
! <TT
     module procedure ig_idx_xxf
     module procedure ig_idx_yxf
     module procedure ig_idx_f
     module procedure ig_idx_jf !<DD> Added
  end interface

  interface ik_idx
     module procedure ik_idx_g
     module procedure ik_idx_gint
     module procedure ik_idx_geint
     module procedure ik_idx_jf
     module procedure ik_idx_lz
     module procedure ik_idx_e
! TT>
     module procedure ik_idx_le
! <TT
     module procedure ik_idx_xxf
     module procedure ik_idx_accel
     module procedure ik_idx_accelx
     module procedure ik_idx_parity
     module procedure ik_idx_f !<DD> Added
  end interface

  interface it_idx
     module procedure it_idx_g
     module procedure it_idx_gint
     module procedure it_idx_geint
     module procedure it_idx_jf
     module procedure it_idx_lz
     module procedure it_idx_e
! TT>
     module procedure it_idx_le
! <TT
     module procedure it_idx_yxf
     module procedure it_idx_accel
     module procedure it_idx_accelx
     module procedure it_idx_f !<DD> Added
  end interface

  interface il_idx
     module procedure il_idx_g
     module procedure il_idx_e
     module procedure il_idx_xxf
     module procedure il_idx_yxf
     module procedure il_idx_accelx
     module procedure il_idx_parity
  end interface

  interface ie_idx
     module procedure ie_idx_g
     module procedure ie_idx_gint
     module procedure ie_idx_lz
     module procedure ie_idx_xxf
     module procedure ie_idx_yxf
     module procedure ie_idx_accelx
     module procedure ie_idx_parity
  end interface

  interface is_idx
     module procedure is_idx_g
     module procedure is_idx_gint
     module procedure is_idx_geint
     module procedure is_idx_lz
     module procedure is_idx_e
! TT>
     module procedure is_idx_le
! <TT
     module procedure is_idx_xxf
     module procedure is_idx_yxf
     module procedure is_idx_accelx
     module procedure is_idx_parity
  end interface

  interface isign_idx
     module procedure isign_idx_xxf
     module procedure isign_idx_yxf
     module procedure isign_idx_e
  end interface

  interface proc_id
     module procedure proc_id_g
     module procedure proc_id_gint
     module procedure proc_id_geint
     module procedure proc_id_f
     module procedure proc_id_jf
     module procedure proc_id_lz
     module procedure proc_id_e
! TT>
     module procedure proc_id_le
! <TT
     module procedure proc_id_xxf
     module procedure proc_id_yxf
     module procedure proc_id_accelx
     module procedure proc_id_parity
  end interface

  interface idx
     module procedure idx_g
     module procedure idx_gint
     module procedure idx_geint
     module procedure idx_f
     module procedure idx_jf
     module procedure idx_lz
     module procedure idx_e
! TT>
     module procedure idx_le
! <TT
     module procedure idx_xxf
     module procedure idx_yxf
     module procedure idx_accelx
     module procedure idx_parity
  end interface

  interface idx_local
     module procedure idx_local_g,      ig_local_g
     module procedure idx_local_gint,   ig_local_gint
     module procedure idx_local_geint,  ig_local_geint
     module procedure idx_local_f,      ig_local_f
     module procedure idx_local_jf,     ig_local_jf
     module procedure idx_local_lz,     ig_local_lz
     module procedure idx_local_e,      ig_local_e
! TT>
     module procedure idx_local_le,      ig_local_le
! <TT
     module procedure idx_local_xxf,    ig_local_xxf
     module procedure idx_local_yxf,    ig_local_yxf
     module procedure idx_local_accelx, ig_local_accelx
     module procedure idx_local_parity, ig_local_parity
  end interface

contains
  subroutine wnml_gs2_layouts(unit)
    implicit none
    integer, intent(in) :: unit
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

    if (initialized_layouts) return
    initialized_layouts = .true.
    
    if (proc0) call read_parameters
    call broadcast_results

  end subroutine init_gs2_layouts

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist, error_unit
    use redistribute, only: opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap
    use mp, only: mp_abort
    implicit none
    integer :: in_file
    namelist /layouts_knobs/ layout, local_field_solve, unbalanced_xxf, &
         max_unbalanced_xxf, unbalanced_yxf, max_unbalanced_yxf, &
         opt_local_copy, opt_redist_nbk, opt_redist_init, intmom_sub, &
         intspec_sub, opt_redist_persist, opt_redist_persist_overlap, fft_measure_plan, &
         fft_use_wisdom, fft_wisdom_file

    local_field_solve = .false.
    unbalanced_xxf = .false.
    unbalanced_yxf = .false.
    opt_local_copy = .false. 
    max_unbalanced_xxf = 0.0
    max_unbalanced_yxf = 0.0
    layout = 'lxyes'
    opt_redist_nbk = .false. !<DD>True=>Use nonblocking redistributes
    opt_redist_persist = .false. !<DD>True=>Use persistent communications in redistributes
    opt_redist_persist_overlap = .false. !<DD>True=>Start comms before doing local copy
    opt_redist_init= .false. !<DD>True=>Use optimised routines to init redist objects
    fft_measure_plan = .true. !<DD>False=>use heuristics rather than measure in fft plan create
    fft_wisdom_file = 'default'  ! Location of fftw wisdom... if left as default, this is set to 
        ! run_name//'.fftw_wisdom', unless overriden by the environment variable
        ! GK_FFTW3_WISDOM. If set to anything other default, overrides  GK_FFTW3_WISDOM
    fft_use_wisdom = .true. ! Try to load and save wisdom about fftw plans to fft_wisdom_file
    intmom_sub=.false.
    intspec_sub=.false.
    in_file=input_unit_exist("layouts_knobs", exist)
    if (exist) read (unit=input_unit("layouts_knobs"), nml=layouts_knobs)
    if (layout.ne.'yxels' .and. layout.ne.'yxles' .and. layout.ne.'lexys'&
    .and. layout.ne.'lxyes' .and. layout.ne.'lyxes' .and. layout.ne.'xyles') &
    then
       write(6,*) "gs2_layouts: read_parameters finds illegal layout=",layout," =>stop"
       call mp_abort("gs2_layouts: read_parameters finds illegal layout=")
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

    !Disable settings if dependent settings not set
    opt_redist_persist=opt_redist_persist.and.opt_redist_nbk
    opt_redist_persist_overlap=opt_redist_persist_overlap.and.opt_redist_persist
    call get_wisdom_file(fft_wisdom_file)
  end subroutine read_parameters
    

  !> Set the fft_wisdom_file based on the value set in the namelist (if any)
  !! and the value in the GK_FFTW3_WISDOM environment variable
  subroutine get_wisdom_file(wisdom_file)
    use file_utils, only: run_name
    character(len=*), intent(inout) :: wisdom_file
    character(len=1000) :: env_wisdom_file
    call get_environment_variable("GK_FFTW3_WISDOM", env_wisdom_file)
    !read (env_wisdom_file,'(I10)') verbosity

    if (wisdom_file .eq. 'default') then
      if (trim(env_wisdom_file) .eq. '') then
        wisdom_file = trim(run_name)//'.fftw_wisdom'
      else
        wisdom_file = env_wisdom_file
      end if
    end if
    !write (*,*) 'wisdom_file', wisdom_file
  end subroutine get_wisdom_file

  subroutine broadcast_results
    use mp, only: broadcast
    use redistribute, only: opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap
    implicit none

    call broadcast (opt_redist_init)
    call broadcast (opt_redist_nbk)
    call broadcast (opt_redist_persist)
    call broadcast (opt_redist_persist_overlap)
    call broadcast (intmom_sub)
    call broadcast (intspec_sub)
    call broadcast (layout)
    call broadcast (local_field_solve)
    call broadcast (unbalanced_xxf)
    call broadcast (unbalanced_yxf)
    call broadcast (max_unbalanced_xxf)
    call broadcast (max_unbalanced_yxf)
    call broadcast (opt_local_copy)
    call broadcast (fft_measure_plan)
    call broadcast (fft_use_wisdom)
    call broadcast (fft_wisdom_file)
  end subroutine broadcast_results

  subroutine check_accel (ntheta0, naky, nlambda, negrid, nspec, nblock)
    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky, ntheta0
    integer, intent (out) :: nblock
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nxfacs, nlfacs
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
       (naky, ntheta0, nlambda, negrid, nspec)
    use mp, only: iproc, nproc, proc0
    use mp, only: split, mp_comm, nproc_comm, rank_comm
    use mp, only: sum_allreduce_sub
! TT>
    use file_utils, only: error_unit
! <TT
    implicit none
    integer, intent (in) :: naky, ntheta0, nlambda, negrid, nspec
!<DD>
    integer :: iglo,ik,it,il,ie,is,col,mycol
    integer :: ik_min,ik_max,it_min,it_max,il_min,il_max,ie_min,ie_max,is_min,is_max,ip
    integer :: nproc_subcomm, rank_subcomm
    integer :: nproc_tmp, idim, tmp
    real :: tmp_r
    integer,dimension(5) :: nproc_dim
    integer,dimension(:),allocatable :: les_kxky_range
    logical,dimension(5) :: dim_divides,dim_local
    logical,dimension(:),allocatable :: it_local, ik_local, il_local, ie_local, is_local
! TT>
# ifdef USE_C_INDEX
    integer :: ierr
    interface
       function init_indices_glo_c (layout)
         integer :: init_indices_glo_c
         character(*) :: layout
       end function init_indices_glo_c
    end interface
# endif
! <TT

    if (initialized_dist_fn_layouts) return
    initialized_dist_fn_layouts = .true.
   
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

!<DD>Calculate constants used in the index lookup routines
!(rather than calculating them at each call as done previously)
!This allows the removal of the select case statements in the index routines
    select case (layout)
    case ('yxels')
       g_lo%ie_ord=3
       g_lo%il_ord=4
       g_lo%it_ord=2
       g_lo%ik_ord=1
       g_lo%is_ord=5
       g_lo%compound_count(1)=1
       g_lo%compound_count(2)=naky
       g_lo%compound_count(3)=naky*ntheta0
       g_lo%compound_count(4)=naky*ntheta0*negrid
       g_lo%compound_count(5)=naky*ntheta0*negrid*nlambda
       g_lo%dim_size(1)=naky
       g_lo%dim_size(2)=ntheta0
       g_lo%dim_size(3)=negrid
       g_lo%dim_size(4)=nlambda
       g_lo%dim_size(5)=nspec
    case ('yxles')
       g_lo%ie_ord=4
       g_lo%il_ord=3
       g_lo%it_ord=2
       g_lo%ik_ord=1
       g_lo%is_ord=5
       g_lo%compound_count(1)=1
       g_lo%compound_count(2)=naky
       g_lo%compound_count(3)=naky*ntheta0
       g_lo%compound_count(4)=naky*ntheta0*nlambda
       g_lo%compound_count(5)=naky*ntheta0*nlambda*negrid
       g_lo%dim_size(1)=naky
       g_lo%dim_size(2)=ntheta0
       g_lo%dim_size(3)=nlambda
       g_lo%dim_size(4)=negrid
       g_lo%dim_size(5)=nspec
    case ('lexys')
       g_lo%ie_ord=2
       g_lo%il_ord=1
       g_lo%it_ord=3
       g_lo%ik_ord=4
       g_lo%is_ord=5
       g_lo%compound_count(1)=1
       g_lo%compound_count(2)=nlambda
       g_lo%compound_count(3)=nlambda*negrid
       g_lo%compound_count(4)=nlambda*negrid*ntheta0
       g_lo%compound_count(5)=nlambda*negrid*ntheta0*naky
       g_lo%dim_size(1)=nlambda
       g_lo%dim_size(2)=negrid
       g_lo%dim_size(3)=ntheta0
       g_lo%dim_size(4)=naky
       g_lo%dim_size(5)=nspec
    case ('lxyes')
       g_lo%ie_ord=4
       g_lo%il_ord=1
       g_lo%it_ord=2
       g_lo%ik_ord=3
       g_lo%is_ord=5
       g_lo%compound_count(1)=1
       g_lo%compound_count(2)=nlambda
       g_lo%compound_count(3)=nlambda*ntheta0
       g_lo%compound_count(4)=nlambda*ntheta0*naky
       g_lo%compound_count(5)=nlambda*ntheta0*naky*negrid
       g_lo%dim_size(1)=nlambda
       g_lo%dim_size(2)=ntheta0
       g_lo%dim_size(3)=naky
       g_lo%dim_size(4)=negrid
       g_lo%dim_size(5)=nspec
    case ('lyxes')
       g_lo%ie_ord=4
       g_lo%il_ord=1
       g_lo%it_ord=3
       g_lo%ik_ord=2
       g_lo%is_ord=5
       g_lo%compound_count(1)=1
       g_lo%compound_count(2)=nlambda
       g_lo%compound_count(3)=nlambda*naky
       g_lo%compound_count(4)=nlambda*naky*ntheta0
       g_lo%compound_count(5)=nlambda*naky*ntheta0*negrid
       g_lo%dim_size(1)=nlambda
       g_lo%dim_size(2)=naky
       g_lo%dim_size(3)=ntheta0
       g_lo%dim_size(4)=negrid
       g_lo%dim_size(5)=nspec
    case ('xyles')
       g_lo%ie_ord=4
       g_lo%il_ord=3
       g_lo%it_ord=1
       g_lo%ik_ord=2
       g_lo%is_ord=5
       g_lo%compound_count(1)=1
       g_lo%compound_count(2)=ntheta0
       g_lo%compound_count(3)=ntheta0*naky
       g_lo%compound_count(4)=ntheta0*naky*nlambda
       g_lo%compound_count(5)=ntheta0*naky*nlambda*negrid
       g_lo%dim_size(1)=ntheta0
       g_lo%dim_size(2)=naky
       g_lo%dim_size(3)=nlambda
       g_lo%dim_size(4)=negrid
       g_lo%dim_size(5)=nspec
    end select
    g_lo%ik_comp=g_lo%compound_count(g_lo%ik_ord)
    g_lo%it_comp=g_lo%compound_count(g_lo%it_ord)
    g_lo%il_comp=g_lo%compound_count(g_lo%il_ord)
    g_lo%ie_comp=g_lo%compound_count(g_lo%ie_ord)
    g_lo%is_comp=g_lo%compound_count(g_lo%is_ord)

    !<DD>Work out min and max of the five compound indices and the locality
    !These are currently only used to slice the arrays used in the velocity space integrals
    !but could well be useful in a number of places.
    !<DD>NOTE: These are typically only useful if the dimensions are split nicely (i.e. we're using
    !a sweet spot) otherwise our local ik section may look like, [naky-1,naky,1,2], so whilst the 
    !min and max values (1,naky) are technically correct, we cannot assume we cover the range 1-naky
    ik_max=0
    ik_min=naky+1
    it_max=0
    it_min=ntheta0+1
    il_max=0
    il_min=nlambda+1
    ie_max=0
    ie_min=negrid+1
    is_max=0
    is_min=nspec+1
    allocate(ik_local(naky),it_local(ntheta0),il_local(nlambda),ie_local(negrid),is_local(nspec))
    ik_local=.false.
    it_local=.false.
    il_local=.false.
    ie_local=.false.
    is_local=.false.
    do iglo=g_lo%llim_proc,g_lo%ulim_proc
       ik=ik_idx(g_lo,iglo)
       it=it_idx(g_lo,iglo)
       il=il_idx(g_lo,iglo)
       ie=ie_idx(g_lo,iglo)
       is=is_idx(g_lo,iglo)
       ik_max=MAX(ik,ik_max)
       it_max=MAX(it,it_max)
       il_max=MAX(il,il_max)
       ie_max=MAX(ie,ie_max)
       is_max=MAX(is,is_max)
       ik_min=MIN(ik,ik_min)
       it_min=MIN(it,it_min)
       il_min=MIN(il,il_min)
       ie_min=MIN(ie,ie_min)
       is_min=MIN(is,is_min)
       ik_local(ik)=.true.
       it_local(it)=.true.
       il_local(il)=.true.
       ie_local(ie)=.true.
       is_local(is)=.true.
    enddo

    !Store dimension locality
    g_lo%y_local=all(ik_local)
    g_lo%x_local=all(it_local)
    g_lo%l_local=all(il_local)
    g_lo%e_local=all(ie_local)
    g_lo%s_local=all(is_local)
    
    !//Deallocate locality data. Note this may actually be useful in some
    !places, so we may want to think about attaching it to the g_lo object
    deallocate(ik_local,it_local,il_local,ie_local,is_local)

    !<DD>Now store min/max in appropriate dimension vars
    !This data isn't currently used but is useful to allow us to
    !replace loops over iglo with nested loops over the 5 generic dimensions
    !which are certain to be in optimised memory order.
    !Here is an example of this.
    !
    !ORIGINAL ROUTINE
    !subroutine Test
    ! use gs2_layouts, only: g_lo, ik_idx,it_idx
    ! use theta_grids, only: ntgrid
    ! use kt_grids, only: naky, ntheta0
    ! implicit none
    ! integer :: iglo, it,ik
    ! real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lim_alloc)::g
    ! real, dimension(-ntgrid:ntgrid,naky,ntheta0) :: total
    ! do iglo=g_lo%llim_proc,g_lim_alloc
    !     it=it_idx(g_lo,iglo)
    !     kt=ik_idx(g_lo,iglo)
    !     total(:,ik,it)=total(:,ik,it)+g(:,1,iglo)+g(:,2,iglo)
    ! enddo
    !end subroutine
    !
    !NEW ROUTINE
    !subroutine Test
    ! use gs2_layouts, only: g_lo, ik_idx,it_idx
    ! use theta_grids, only: ntgrid
    ! use kt_grids, only: naky, ntheta0
    ! implicit none
    ! integer :: iglo, d1,d2,d3,d4,d5
    ! integer, pointer :: it,ik
    ! integer, dimension(5),target :: dind
    ! real, dimension(-ntgrid:ntgrid,1:2,g_lo%llim_proc:g_lim_alloc)::g
    ! real, dimension(-ntgrid:ntgrid,naky,ntheta0) :: total
    ! EQUIVALENCE(dind(1),d1),(dind(2),d2),(dind(3),d3),(dind(4),d4),(dind(5),d5) !<--This just allows us to access d1..d5 using array notation
    ! !Setup ik and it to point to appropriate d1..d5
    ! ik=>dind(g_lo%ik_ord)
    ! it=>dind(g_lo%it_ord)
    ! iglo=g_lo%llim_proc
    ! do d5=g_lo%d5_min,g_lo%d5_max
    !    do d4=.... ; do d3=.... ; do d2=.... ; do d1=....
    !          total(:,ik,it)=total(:,ik,it)+g(:,1,iglo)+g(:,2,iglo)
    !          iglo=iglo+1
    !    enddo;enddo;enddo;enddo
    ! enddo
    !end subroutine
    select case (g_lo%ik_ord)
    case (1)
       g_lo%d1_min=ik_min
       g_lo%d1_max=ik_max
    case (2)
       g_lo%d2_min=ik_min
       g_lo%d2_max=ik_max
    case (3)
       g_lo%d3_min=ik_min
       g_lo%d3_max=ik_max
    case (4)
       g_lo%d4_min=ik_min
       g_lo%d4_max=ik_max
    case (5)
       g_lo%d5_min=ik_min
       g_lo%d5_max=ik_max
    end select

    select case (g_lo%it_ord)
    case (1)
       g_lo%d1_min=it_min
       g_lo%d1_max=it_max
    case (2)
       g_lo%d2_min=it_min
       g_lo%d2_max=it_max
    case (3)
       g_lo%d3_min=it_min
       g_lo%d3_max=it_max
    case (4)
       g_lo%d4_min=it_min
       g_lo%d4_max=it_max
    case (5)
       g_lo%d5_min=it_min
       g_lo%d5_max=it_max
    end select
       
    select case (g_lo%il_ord)
    case (1)
       g_lo%d1_min=il_min
       g_lo%d1_max=il_max
    case (2)
       g_lo%d2_min=il_min
       g_lo%d2_max=il_max
    case (3)
       g_lo%d3_min=il_min
       g_lo%d3_max=il_max
    case (4)
       g_lo%d4_min=il_min
       g_lo%d4_max=il_max
    case (5)
       g_lo%d5_min=il_min
       g_lo%d5_max=il_max
    end select

    select case (g_lo%ie_ord)
    case (1)
       g_lo%d1_min=ie_min
       g_lo%d1_max=ie_max
    case (2)
       g_lo%d2_min=ie_min
       g_lo%d2_max=ie_max
    case (3)
       g_lo%d3_min=ie_min
       g_lo%d3_max=ie_max
    case (4)
       g_lo%d4_min=ie_min
       g_lo%d4_max=ie_max
    case (5)
       g_lo%d5_min=ie_min
       g_lo%d5_max=ie_max
    end select

    select case (g_lo%is_ord)
    case (1)
       g_lo%d1_min=is_min
       g_lo%d1_max=is_max
    case (2)
       g_lo%d2_min=is_min
       g_lo%d2_max=is_max
    case (3)
       g_lo%d3_min=is_min
       g_lo%d3_max=is_max
    case (4)
       g_lo%d4_min=is_min
       g_lo%d4_max=is_max
    case (5)
       g_lo%d5_min=is_min
       g_lo%d5_max=is_max
    end select

    g_lo%ik_min=ik_min
    g_lo%ik_max=ik_max
    g_lo%it_min=it_min
    g_lo%it_max=it_max
    g_lo%il_min=il_min
    g_lo%il_max=il_max
    g_lo%ie_min=ie_min
    g_lo%ie_max=ie_max
    g_lo%is_min=is_min
    g_lo%is_max=is_max
    
    !Work out if x comes before y in layout
    do iglo=1,LEN_TRIM(layout)
       if(layout(iglo:iglo).eq.'x') then
          g_lo%x_before_y=.true.
          exit
       elseif (layout(iglo:iglo).eq.'y') then
          g_lo%x_before_y=.false.
          exit
       endif
    enddo

    !Now work out how the compound dimension splits amongst processors
    nproc_tmp=nproc
    !/Initialise
    nproc_dim=-1
    dim_divides=.false.
    dim_local=.false.
    do idim=5,1,-1

       !If we are down to 1 proc then everything else is local
       !so say dim_divides and exit loop
       !If nproc_tmp is 0 then this indicates the last dimension
       !had more entries than procs available
       if(nproc_tmp.le.1) then
          dim_divides(idim:1:-1)=.true.
          dim_local(idim:1:-1)=.true.
          exit
       endif
       
       !Does the current dimension divide evenly
       !amongst the remaining processors
       tmp_r=(1.0*nproc_tmp)/g_lo%dim_size(idim)

       !Check if the current dimension is larger than the number of processors
       !if so then we can still divide nicely but can't rely on division being integer
       !Invert division to check!
       if(int(tmp_r).eq.0) then
          tmp_r=g_lo%dim_size(idim)/(1.0*nproc_tmp)
       endif

       !Store logical to indicate if things divide perfectly
       dim_divides(idim)=((tmp_r-int(tmp_r)).eq.0)

       !If this dimension doesn't divide perfectly then we exit
       !loop and leave other dims as "un-initialised"
       if(.not.dim_divides(idim)) then
          nproc_dim(idim)=int(tmp_r-int(tmp_r))
          exit
       endif

       !If we do split nicely then store how many procs there are
       !to deal with the subsequent dims
       !If this is zero then it indicates a dimension that is not
       !fully distributed.
       nproc_dim(idim)=nproc_tmp/g_lo%dim_size(idim)

       !Update the available processor count
       nproc_tmp=nproc_dim(idim)
    enddo

    !Now print some stuff for debugging
!     if (proc0) then
!        write(6,'("Idim | Divides | Nproc_per_block | Size | Local")')
!        do idim=5,1,-1
!           write(6,'(I0,6(" "),L,9(" "),I0,16(" "),I0,8(" "),L)') idim,dim_divides(idim),nproc_dim(idim),g_lo%dim_size(idim),dim_local(idim)
!        enddo
!     endif

    !Now use the dimension splitting to work out if the various subcommunicators are
    !allowed.
    if(.not.(dim_divides(g_lo%is_ord).and.dim_divides(g_lo%it_ord).and.dim_divides(g_lo%ik_ord))) then
       if(intmom_sub.and.proc0) write(error_unit(),'("Disabling intmom_sub -- is,it,ik split nicely ? ",3(L1," "))') dim_divides(g_lo%is_ord),dim_divides(g_lo%it_ord),dim_divides(g_lo%ik_ord)
       intmom_sub=.false.
    endif
    if(.not.(dim_divides(g_lo%it_ord).and.dim_divides(g_lo%ik_ord))) then
       if(intspec_sub.and.proc0) write(error_unit(),'("Disabling intspec_sub -- it,ik split nicely ? ",2(L1," "))') dim_divides(g_lo%it_ord),dim_divides(g_lo%ik_ord)
       intspec_sub=.false.
    endif
    !Note that because we currently need to gather amongst LES blocks at the end of
    !integrate_species we need to make sure that LES blocks are also sensible
    !This could be removed if we don't want to gather in integrate_species
    if(.not.(dim_divides(g_lo%il_ord).and.dim_divides(g_lo%ie_ord).and.dim_divides(g_lo%is_ord))) then
       if(intspec_sub.and.proc0) write(error_unit(),'("Disabling intspec_sub -- il,ie,is split nicely ? ",3(L1," "))') dim_divides(g_lo%il_ord),dim_divides(g_lo%ie_ord),dim_divides(g_lo%is_ord)
       intspec_sub=.false.
    endif

    !Update various flags based on locality
    !/Don't use subcommunicator and  gather for integrate_species if x and y are entirely local
    if((g_lo%x_local.and.g_lo%y_local)) then
       if(proc0.and.intspec_sub) write(error_unit(),'("Disabling intspec_sub as x and y are local")')
       intspec_sub=.false.
    endif

    !We could(should) use a logical allreduction here to ensure that ALL processors have the same flag settings
    !(if not then we'll get a program hang in the following splits).
    !The logic above **should** ensure everyone has the right answer anyway.

    !<DD>Now work out which xy block we live in for velocity space comms
    !Set sub communicators to mp_comm if we want to keep collective comms global.
    !NOTE: We might want to consider using a communicator type so that we can attach
    !some useful information to the communicators such as number of procs in comm,
    !this procs rank in comm, if this comm is equivalent to mp_comm (i.e. if split
    !doesn't actually split the global comm) etc.
    !This could help reduce the amount of mp(i) calls placed throughout the code
    !as things like the number of procs and the local rank are useful when using the sub
    !comms but don't currently have a logical place to be stored and hence tend to be
    !(re)calculated when needed.
    if(intspec_sub) then
       !This is for unique xy blocks
       col=1
       mycol=0
       do it=1,ntheta0
          do ik=1,naky
             !This is inefficient but who cares
             if(it_min.eq.it.and.ik_min.eq.ik) then
                mycol=col
                exit
             endif
             col=col+1
          enddo
       enddo

       !<DD>Now split the global communicator
       call split(mycol,g_lo%xyblock_comm)

       !This is for unique les blocks
       col=1
       mycol=0
       do is=1,nspec
          do il=1,nlambda
             do ie=1,negrid
                !This is inefficient but who cares
                if(il_min.eq.il.and.ie_min.eq.ie.and.is_min.eq.is) then
                   mycol=col
                   exit
                endif
                col=col+1
             enddo
          enddo
       enddo

       !Now split
       call split(mycol,g_lo%lesblock_comm)

    else
       g_lo%xyblock_comm=mp_comm
       g_lo%lesblock_comm=mp_comm
    endif

    if(intmom_sub)then
       !<DD>NOTE: A simple test that all blocks are equal would be:
       !tmp=it_max*ik_max
       !call max_allreduce(tmp,g_lo%xyblock_comm)
       !if(tmp.ne.(it_max*ik_max)) then
       !   print'(a)',"ERROR: In subcomm g_lo%xyblock_comm"
       !   print'(a,I0,a)',"Proc with global rank ",iproc," doesn't have the same it*ik extent as other procs in its subcomm"
       !   call mp_abort("Sub-comm creation error -- Consider changing the number of procs/number of xy gridpoints or setting intmom_sub=.false.")
       !endif
       !Similar tests can be done for all subcomms created here
       !</DD>

       !<DD>Now work out which xys block we live in for velocity space comms
       !This is for unique xys blocks
       col=1
       mycol=0
       do it=1,ntheta0
          do ik=1,naky
             do is=1,nspec
                !This is inefficient but who cares
                if(it_min.eq.it.and.ik_min.eq.ik.and.is_min.eq.is) then
                   mycol=col
                   exit
                endif
                col=col+1
             enddo
          enddo
       enddo

       !<DD>Now split the global communicator
       call split(mycol,g_lo%xysblock_comm)
    else
       g_lo%xysblock_comm=mp_comm
    endif

    !Now allocate and fill various arrays
    if(intspec_sub.and.(.not.allocated(g_lo%les_kxky_range))) then
       !->Kx-Ky range for procs in lesblock_comm sub comm
       !Get number of procs in sub-comm
       call nproc_comm(g_lo%lesblock_comm,nproc_subcomm)
       !Get rank of proc in sub-comm
       call rank_comm(g_lo%lesblock_comm,rank_subcomm)
       !Allocate array | This will store the start and stop indices 
       !for flattened kxky indices in velocity space independent (e.g. field)
       !variables
       allocate(g_lo%les_kxky_range(2,0:nproc_subcomm-1))
       allocate(les_kxky_range(0:nproc_subcomm*2-1))
       !Initialise as we do sum_allreduce later
       !can remove if we use a gather
       les_kxky_range=0
       
       !Now calculate our personal kxky range
       it_min=it_idx(g_lo,g_lo%llim_proc)
       it_max=it_idx(g_lo,g_lo%ulim_alloc)
       ik_min=ik_idx(g_lo,g_lo%llim_proc)
       ik_max=ik_idx(g_lo,g_lo%ulim_alloc)
       
       if(g_lo%x_before_y) then
          les_kxky_range(rank_subcomm*2)=it_min-1+ntheta0*(ik_min-1)
          tmp=les_kxky_range(rank_subcomm*2)
          do iglo=g_lo%llim_proc,g_lo%ulim_alloc
             tmp=MAX(tmp,it_idx(g_lo,iglo)-1+ntheta0*(ik_idx(g_lo,iglo)-1))
          enddo
          les_kxky_range(rank_subcomm*2+1)=tmp
       else
          les_kxky_range(rank_subcomm*2)=ik_min-1+naky*(it_min-1)
          tmp=les_kxky_range(rank_subcomm*2)
          do iglo=g_lo%llim_proc,g_lo%ulim_alloc
             tmp=MAX(tmp,ik_idx(g_lo,iglo)-1+naky*(it_idx(g_lo,iglo)-1))
          enddo
          les_kxky_range(rank_subcomm*2+1)=tmp
       endif
       !Now gather values and store in g_lo, sum_allreduce would also work
       !    call mpi_allgather(les_kxky_range(rank_subcomm*2:rank_subcomm*2+1),2,MPI_INTEGER,&
       !         les_kxky_range,2,MPI_INTEGER,g_lo%lesblock_comm)
       !Should really use gather but this is easier
       call sum_allreduce_sub(les_kxky_range,g_lo%lesblock_comm)
       do ip=0,nproc_subcomm-1
          g_lo%les_kxky_range(:,ip)=les_kxky_range(ip*2:ip*2+1)
       enddo
       !Free memory
       deallocate(les_kxky_range)
    endif
!</DD>

!Note: gint_lo isn't used anywhere!
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
    
!Note: geint_lo isn't used anywhere!
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

! TT>
# ifdef USE_C_INDEX
    ierr = init_indices_glo_c (layout)
    if (ierr /= 0) &
         & write (error_unit(),*) 'ERROR: layout not found: ', trim(layout)
# endif
! <TT

  end subroutine init_dist_fn_layouts

  subroutine is_kx_local(negrid, nspec, nlambda, naky, kx_local)  

    use mp, only: nproc
    implicit none
    integer, intent (in) :: negrid, nspec, nlambda, naky
    logical, intent (out) :: kx_local
    integer, dimension(:,:), allocatable :: facs
    integer :: nsfacs, nefacs, nyfacs, nlfacs
    integer :: i

    kx_local = .true.

    allocate (facs(max(nspec,negrid,naky,nlambda)/2+1,4))

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

    case ('xyles')

       call factors (nspec,   nsfacs, facs(:,1))
       call factors (negrid,  nefacs, facs(:,2))
       call factors (nlambda, nlfacs, facs(:,3))
       call factors (naky, nyfacs, facs(:,4))

       do i=1,nsfacs-1
          if (nproc == facs(i,1)) goto 100
       end do
       
       do i=1,nefacs-1
          if (nproc == facs(i,2)*nspec) goto 100
       end do

       do i=1,nlfacs-1
          if (nproc == facs(i,3)*negrid*nspec) goto 100
       end do
          
       do i=1,nyfacs-1
          if (nproc == facs(i,4)*nlambda*negrid*nspec) goto 100
       end do
          
    end select

    kx_local = .false.    
                       
100 deallocate (facs)

  end subroutine is_kx_local

! TT>
# ifdef USE_C_INDEX
  function is_idx_g (lo, i)
# else
  elemental function is_idx_g (lo, i)
# endif
! <TT
    implicit none
    integer :: is_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
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
    is_idx_g=1+mod((i-lo%llim_world)/lo%is_comp,lo%nspec)
# endif
! <TT

  end function is_idx_g

! TT>
# ifdef USE_C_INDEX
  function il_idx_g (lo, i)
# else
  elemental function il_idx_g (lo, i)
# endif
! <TT
    implicit none
    integer :: il_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function il_idx_g_c (lo,num)
         use layouts_type, only: g_layout_type
         integer :: il_idx_g_c
         type (g_layout_type) :: lo
         integer :: num
       end function il_idx_g_c
    end interface
    il_idx_g = il_idx_g_c (lo,i)
# else
    il_idx_g=1+mod((i-lo%llim_world)/lo%il_comp,lo%nlambda)
# endif
! <TT

  end function il_idx_g

! TT>
# ifdef USE_C_INDEX
  function ie_idx_g (lo, i)
# else
  elemental function ie_idx_g (lo, i)
# endif
! <TT
    implicit none
    integer :: ie_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function ie_idx_g_c (lo,num)
         use layouts_type, only: g_layout_type
         integer :: ie_idx_g_c
         type (g_layout_type) :: lo
         integer :: num
       end function ie_idx_g_c
    end interface
    ie_idx_g = ie_idx_g_c (lo,i)
# else
    ie_idx_g=1+mod((i-lo%llim_world)/lo%ie_comp,lo%negrid)
# endif
! <TT

  end function ie_idx_g

! TT>
# ifdef USE_C_INDEX
  function it_idx_g (lo, i)
# else
  elemental function it_idx_g (lo, i)
# endif
! <TT
    implicit none
    integer :: it_idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    
! TT>
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
    it_idx_g=1+mod((i-lo%llim_world)/lo%it_comp,lo%ntheta0)
# endif
! <TT

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
! TT>
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
    ik_idx_g=1+mod((i-lo%llim_world)/lo%ik_comp,lo%naky)
# endif
! <TT

  end function ik_idx_g

  elemental function proc_id_g (lo, i)
    implicit none
    integer :: proc_id_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_g = i/lo%blocksize

  end function proc_id_g

! TT>
# ifdef USE_C_INDEX
  function idx_g (lo, ik, it, il, ie, is)
# else
  elemental function idx_g (lo, ik, it, il, ie, is)
# endif
! <TT
    implicit none
    integer :: idx_g
    type (g_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, it, il, ie, is
! TT>
# ifdef USE_C_INDEX
    interface
       function idx_g_c (lo,ik,it,il,ie,is)
         use layouts_type, only: g_layout_type
         integer :: idx_g_c
         type (g_layout_type) :: lo
         integer :: ik,it,il,ie,is
       end function idx_g_c
    end interface
    idx_g = idx_g_c (lo, ik, it, il, ie, is)
# else
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
    case ('xyles')
       idx_g = it-1 + lo%ntheta0*(ik-1 + lo%naky*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1))))
    end select
# endif
! <TT

  end function idx_g

! TT>
# ifdef USE_C_INDEX
  function idx_local_g (lo, ik, it, il, ie, is)
# else
  elemental function idx_local_g (lo, ik, it, il, ie, is)
# endif
! <TT
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
    case ('xyles')
       is_idx_gint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky/lo%negrid, lo%nspec)
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
    case ('xyles')
       ie_idx_gint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%negrid)
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
    case ('xyles')
       it_idx_gint = 1 + mod((i - lo%llim_world), lo%ntheta0)
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
    case ('xyles')
       ik_idx_gint = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
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
    case ('xyles')
       idx_gint = it-1 + lo%ntheta0*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1)))
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
    case ('xyles')
       is_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0/lo%naky, lo%nspec)
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
    case ('xyles')
       it_idx_geint = 1 + mod(i - lo%llim_world, lo%ntheta0)
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
    case ('xyles')
       ik_idx_geint = 1 + mod((i - lo%llim_world)/lo%ntheta0, lo%naky)
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
    case ('xyles')
       idx_geint = it-1 + lo%ntheta0*(ik-1 + lo%naky*(is-1))
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
    integer :: i, utmp, btmp

    if (initialized_fields_layouts) return
    initialized_fields_layouts = .true.
    
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

  subroutine finish_fields_layouts
    implicit none

    integer :: i,f_size

    if(allocated(f_lo))then
       f_size = size(f_lo)
       do i=1,f_size
          if(associated(f_lo(i)%ik)) deallocate (f_lo(i)%ik)
          if(associated(f_lo(i)%it)) deallocate (f_lo(i)%it)
       enddo
       deallocate(f_lo)    
    endif
    initialized_fields_layouts = .false.
  end subroutine finish_fields_layouts

! HJL < Finish routine so that gs2 can be restarted with a different
! number of processors by trinity
  subroutine finish_layouts
    implicit none
    call finish_fields_layouts
    call finish_jfields_layouts

    initialized_layouts = .false.
    initialized_dist_fn_layouts = .false.
    initialized_le_layouts = .false.
    initialized_energy_layouts = .false.
    initialized_lambda_layouts = .false.
    initialized_parity_layouts = .false.
    initialized_accel_transform_layouts = .false.
    initialized_x_transform = .false.
    initialized_y_transform = .false.

  end subroutine finish_layouts
! > HJL

!<DD>Added
  function ik_idx_f (lo, i)
    implicit none
    integer :: ik_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_f = lo%ik(im_idx(lo,i),in_idx(lo,i))
  end function ik_idx_f

!<DD>Added
  function it_idx_f (lo, i)
    implicit none
    integer :: it_idx_f
    type (f_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_f = lo%it(im_idx(lo,i),in_idx(lo,i))
  end function it_idx_f

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
    integer :: jlo
    if (initialized_jfields_layouts) return
    initialized_jfields_layouts = .true.
    
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

    !Now work out the range of dimensions stored on this proc !<DD>
    jf_lo%ig_max=-jf_lo%ntgrid-1
    jf_lo%ig_min=jf_lo%ntgrid+1
    jf_lo%ifield_max=0
    jf_lo%ifield_min=nfield+1
    jf_lo%ik_max=0
    jf_lo%ik_min=naky+1
    jf_lo%it_max=0
    jf_lo%it_min=0
    do jlo=jf_lo%llim_proc,jf_lo%ulim_alloc
       jf_lo%ig_max=MAX(jf_lo%ig_max,ig_idx(jf_lo,jlo))
       jf_lo%ig_min=MIN(jf_lo%ig_min,ig_idx(jf_lo,jlo))
       jf_lo%ifield_max=MAX(jf_lo%ifield_max,ifield_idx(jf_lo,jlo))
       jf_lo%ifield_min=MIN(jf_lo%ifield_min,ifield_idx(jf_lo,jlo))
       jf_lo%ik_max=MAX(jf_lo%ik_max,ik_idx(jf_lo,jlo))
       jf_lo%ik_min=MIN(jf_lo%ik_min,ik_idx(jf_lo,jlo))
       jf_lo%it_max=MAX(jf_lo%it_max,it_idx(jf_lo,jlo))
       jf_lo%it_min=MIN(jf_lo%it_min,it_idx(jf_lo,jlo))
    enddo

    allocate (ij(jf_lo%llim_proc:jf_lo%ulim_alloc))
    allocate (mj(jf_lo%llim_proc:jf_lo%ulim_alloc))
    allocate (dj(i_class,jf_lo%llim_proc:jf_lo%ulim_alloc))
    ij = 1  ; mj = 1;  dj = 0
    
  end subroutine init_jfields_layouts

  subroutine finish_jfields_layouts
    implicit none
    if(allocated(ij)) deallocate(ij)
    if(allocated(mj)) deallocate(mj)
    if(allocated(dj)) deallocate(dj)
    initialized_jfields_layouts = .false.
  end subroutine finish_jfields_layouts

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

  !Here if is the compound theta*field index
  function if_idx_jf (lo, i)
    implicit none
    integer :: if_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    if_idx_jf = 1 + mod(i - lo%llim_world, lo%nindex)
  end function if_idx_jf

  !<DD>Added
  function ig_idx_jf (lo, i)
    implicit none
    integer :: ig_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_jf = -lo%ntgrid + mod(i - lo%llim_world, lo%ntgrid*2+1)
  end function ig_idx_jf

  !<DD>Added
  function ifield_idx_jf (lo, i)
    implicit none
    integer :: ifield_idx_jf
    type (jf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ifield_idx_jf = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid+1), lo%nfield)
  end function ifield_idx_jf

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

! TT>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lambda-Energy layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_le_layouts (ntgrid, naky, ntheta0, nspec)
    use mp, only: iproc, nproc
    use file_utils, only: error_unit
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nspec
    integer :: ik_min, it_min, ig_min, is_min
    integer :: ik_max, it_max, ig_max, is_max
    integer :: ik, it, ig, is, ile

# ifdef USE_C_INDEX
    integer :: ierr
    interface
       function init_indices_lelo_c (layout)
         integer :: init_indices_lelo_c
         character(*) :: layout
       end function init_indices_lelo_c
    end interface
# endif

    if (initialized_le_layouts) return
    initialized_le_layouts = .true.
    
    le_lo%iproc = iproc
    le_lo%ntgrid = ntgrid
    le_lo%ntgridtotal=(2*ntgrid+1)
    le_lo%ntheta0 = ntheta0
    le_lo%naky = naky
    le_lo%nspec = nspec
    le_lo%llim_world = 0
    le_lo%ulim_world = (2*ntgrid+1) * naky * ntheta0 * nspec - 1
    le_lo%blocksize = le_lo%ulim_world / nproc + 1
    le_lo%llim_proc = le_lo%blocksize * iproc
    le_lo%ulim_proc &
         = min(le_lo%ulim_world, le_lo%llim_proc + le_lo%blocksize - 1)
    le_lo%ulim_alloc = max(le_lo%llim_proc, le_lo%ulim_proc)

!<DD>See g_lo init
    select case (layout)
    case ('yxels')
       le_lo%ig_ord=1
       le_lo%it_ord=3
       le_lo%ik_ord=2
       le_lo%is_ord=4
       le_lo%compound_count(1)=1
       le_lo%compound_count(2)=le_lo%ntgridtotal
       le_lo%compound_count(3)=le_lo%compound_count(2)*naky
       le_lo%compound_count(4)=le_lo%compound_count(3)*ntheta0
    case ('yxles')
       le_lo%ig_ord=1
       le_lo%it_ord=3
       le_lo%ik_ord=2
       le_lo%is_ord=4
       le_lo%compound_count(1)=1
       le_lo%compound_count(2)=le_lo%ntgridtotal
       le_lo%compound_count(3)=le_lo%compound_count(2)*naky
       le_lo%compound_count(4)=le_lo%compound_count(3)*ntheta0
    case ('lexys')
       le_lo%ig_ord=1
       le_lo%it_ord=2
       le_lo%ik_ord=3
       le_lo%is_ord=4
       le_lo%compound_count(1)=1
       le_lo%compound_count(2)=le_lo%ntgridtotal
       le_lo%compound_count(3)=le_lo%compound_count(2)*ntheta0
       le_lo%compound_count(4)=le_lo%compound_count(3)*naky
    case ('lxyes')
       le_lo%ig_ord=1
       le_lo%it_ord=2
       le_lo%ik_ord=3
       le_lo%is_ord=4
       le_lo%compound_count(1)=1
       le_lo%compound_count(2)=le_lo%ntgridtotal
       le_lo%compound_count(3)=le_lo%compound_count(2)*ntheta0
       le_lo%compound_count(4)=le_lo%compound_count(3)*naky
    case ('lyxes')
       le_lo%ig_ord=1
       le_lo%it_ord=3
       le_lo%ik_ord=2
       le_lo%is_ord=4
       le_lo%compound_count(1)=1
       le_lo%compound_count(2)=le_lo%ntgridtotal
       le_lo%compound_count(3)=le_lo%compound_count(2)*naky
       le_lo%compound_count(4)=le_lo%compound_count(3)*ntheta0
    case ('xyles')
       le_lo%ig_ord=1
       le_lo%it_ord=2
       le_lo%ik_ord=3
       le_lo%is_ord=4
       le_lo%compound_count(1)=1
       le_lo%compound_count(2)=le_lo%ntgridtotal
       le_lo%compound_count(3)=le_lo%compound_count(2)*ntheta0
       le_lo%compound_count(4)=le_lo%compound_count(3)*naky
    end select
    le_lo%ig_comp=le_lo%compound_count(le_lo%ig_ord)
    le_lo%ik_comp=le_lo%compound_count(le_lo%ik_ord)
    le_lo%it_comp=le_lo%compound_count(le_lo%it_ord)
    le_lo%is_comp=le_lo%compound_count(le_lo%is_ord)
    !<DD>Work out min and max of the five compound indices
    !These are currently only used to slice the arrays used in the velocity space integrals
    !but could well be useful in a number of places.
    ik_max=0
    ik_min=naky+1
    it_max=0
    it_min=ntheta0+1
    ig_max=-ntgrid-1
    ig_min=ntgrid+1
    is_max=0
    is_min=nspec+1
    do ile=le_lo%llim_proc,le_lo%ulim_proc
       ik=ik_idx(le_lo,ile)
       it=it_idx(le_lo,ile)
       ig=ig_idx(le_lo,ile)
       is=is_idx(le_lo,ile)
       ik_max=MAX(ik,ik_max)
       it_max=MAX(it,it_max)
       ig_max=MAX(ig,ig_max)
       is_max=MAX(is,is_max)
       ik_min=MIN(ik,ik_min)
       it_min=MIN(it,it_min)
       ig_min=MIN(ig,ig_min)
       is_min=MIN(is,is_min)
    enddo

    !Store dimension locality
    le_lo%x_local=(it_min.eq.1).and.(it_max.eq.ntheta0)
    le_lo%y_local=(ik_min.eq.1).and.(ik_max.eq.naky)
    le_lo%t_local=(ig_min.eq.-ntgrid).and.(ig_max.eq.ntgrid)
    le_lo%s_local=(is_min.eq.1).and.(is_max.eq.nspec)
!</DD>

# ifdef USE_C_INDEX
    ierr = init_indices_lelo_c (layout)
    if (ierr /= 0) &
         & write (error_unit(),*) 'ERROR: layout not found: ', trim(layout)
# endif
  end subroutine init_le_layouts

  elemental function is_idx_le (lo, i)
    implicit none
    integer :: is_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_le = 1 + mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_le

# ifdef USE_C_INDEX
  function it_idx_le (lo, i)
# else
  elemental function it_idx_le (lo, i)
# endif
    implicit none
    integer :: it_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
# ifdef USE_C_INDEX
    interface
       function it_idx_le_c (lo,num)
         use layouts_type, only: le_layout_type
         integer :: it_idx_le_c
         type (le_layout_type) :: lo
         integer :: num
       end function it_idx_le_c
    end interface
    it_idx_le = it_idx_le_c (lo,i)
# else
    it_idx_le = 1 + mod((i - lo%llim_world)/lo%it_comp, lo%ntheta0)
# endif
  end function it_idx_le

# ifdef USE_C_INDEX
  function ik_idx_le (lo, i)
# else
  elemental function ik_idx_le (lo, i)
# endif
    implicit none
    integer :: ik_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
# ifdef USE_C_INDEX
    interface
       function ik_idx_le_c (lo,num)
         use layouts_type, only: le_layout_type
         integer :: ik_idx_le_c
         type (le_layout_type) :: lo
         integer :: num
       end function ik_idx_le_c
    end interface
    ik_idx_le = ik_idx_le_c (lo,i)
# else
    ik_idx_le = 1 + mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
# endif
  end function ik_idx_le

  elemental function ig_idx_le (lo, i)
    implicit none
    integer :: ig_idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_le = -lo%ntgrid + mod((i - lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_le

# ifdef USE_C_INDEX
  function idx_le (lo, ig, ik, it, is)
# else
  elemental function idx_le (lo, ig, ik, it, is)
# endif
    implicit none
    integer :: idx_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, is
# ifdef USE_C_INDEX
    interface
       function idx_le_c (lo, ig, ik, it, is)
         use layouts_type, only: le_layout_type
         integer :: idx_le_c
         type (le_layout_type) :: lo
         integer :: ig, ik, it, is
       end function idx_le_c
    end interface
    idx_le = idx_le_c (lo, ig, ik, it, is)
# else
    select case (layout)
    case ('yxels')
       idx_le = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 &
            + lo%naky*(it-1 + lo%ntheta0*(is-1)))
    case ('yxles')
       idx_le = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 &
            + lo%naky*(it-1 + lo%ntheta0*(is-1)))
    case ('lexys')
       idx_le = ig+lo%ntgrid + (2*lo%ntgrid+1)*(it-1 &
            + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
    case ('lxyes')
       idx_le = ig+lo%ntgrid + (2*lo%ntgrid+1)*(it-1 &
            + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
    case ('lyxes')
       idx_le = ig+lo%ntgrid + (2*lo%ntgrid+1)*(ik-1 &
            + lo%naky*(it-1 + lo%ntheta0*(is-1)))
    case ('xyles')
       idx_le = ig+lo%ntgrid + (2*lo%ntgrid+1)*(it-1 &
            + lo%ntheta0*(ik-1 + lo%naky*(is-1)))
    end select
# endif
  end function idx_le

  elemental function proc_id_le (lo, i)
    implicit none
    integer :: proc_id_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_le = i / lo%blocksize

  end function proc_id_le

# ifdef USE_C_INDEX
  function idx_local_le (lo, ig, ik, it, is)
# else
  elemental function idx_local_le (lo, ig, ik, it, is)
# endif
    implicit none
    logical :: idx_local_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, is

    idx_local_le = idx_local(lo, idx(lo, ig, ik, it, is))
  end function idx_local_le

  elemental function ig_local_le (lo, ig)
    implicit none
    logical :: ig_local_le
    type (le_layout_type), intent (in) :: lo
    integer, intent (in) :: ig

    ig_local_le = lo%iproc == proc_id(lo, ig)
  end function ig_local_le
! <TT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Energy layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_energy_layouts &
       (ntgrid, naky, ntheta0, nlambda, nspec)
    use mp, only: iproc, nproc
! TT>
    use file_utils, only: error_unit
! <TT
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, nspec
! TT>
# ifdef USE_C_INDEX
    integer :: ierr
    interface
       function init_indices_elo_c (layout)
         integer :: init_indices_elo_c
         character(*) :: layout
       end function init_indices_elo_c
    end interface
# endif
! <TT

    if (initialized_energy_layouts) return
    initialized_energy_layouts = .true.
    
    e_lo%iproc = iproc
    e_lo%ntgrid = ntgrid
    e_lo%nsign = 2
    e_lo%ntheta0 = ntheta0
    e_lo%naky = naky
    e_lo%nspec = nspec
    e_lo%nlambda = nlambda
    e_lo%ntgridtotal=2*ntgrid+1
    e_lo%llim_world = 0
    e_lo%ulim_world = (2*ntgrid+1)*naky*ntheta0*nlambda*nspec*2 - 1
    e_lo%blocksize = e_lo%ulim_world/nproc + 1
    e_lo%llim_proc = e_lo%blocksize*iproc
    e_lo%ulim_proc &
         = min(e_lo%ulim_world, e_lo%llim_proc + e_lo%blocksize - 1)
    e_lo%ulim_alloc = max(e_lo%llim_proc, e_lo%ulim_proc)

!<DD>See g_lo init
    select case (layout)
    case ('yxels')
       e_lo%ig_ord=1
       e_lo%isgn_ord=2
       e_lo%il_ord=5
       e_lo%it_ord=4
       e_lo%ik_ord=3
       e_lo%is_ord=6
       e_lo%compound_count(1)=1
       e_lo%compound_count(2)=e_lo%ntgridtotal
       e_lo%compound_count(3)=e_lo%compound_count(2)*2
       e_lo%compound_count(4)=e_lo%compound_count(3)*naky
       e_lo%compound_count(5)=e_lo%compound_count(4)*ntheta0
       e_lo%compound_count(6)=e_lo%compound_count(5)*nlambda
    case ('yxles')
       e_lo%ig_ord=1
       e_lo%isgn_ord=2
       e_lo%il_ord=5
       e_lo%it_ord=4
       e_lo%ik_ord=3
       e_lo%is_ord=6
       e_lo%compound_count(1)=1
       e_lo%compound_count(2)=e_lo%ntgridtotal
       e_lo%compound_count(3)=e_lo%compound_count(2)*2
       e_lo%compound_count(4)=e_lo%compound_count(3)*naky
       e_lo%compound_count(5)=e_lo%compound_count(4)*ntheta0
       e_lo%compound_count(6)=e_lo%compound_count(5)*nlambda
    case ('lexys')
       e_lo%ig_ord=1
       e_lo%isgn_ord=2
       e_lo%il_ord=3
       e_lo%it_ord=4
       e_lo%ik_ord=5
       e_lo%is_ord=6
       e_lo%compound_count(1)=1
       e_lo%compound_count(2)=e_lo%ntgridtotal
       e_lo%compound_count(3)=e_lo%compound_count(2)*2
       e_lo%compound_count(4)=e_lo%compound_count(3)*nlambda
       e_lo%compound_count(5)=e_lo%compound_count(4)*ntheta0
       e_lo%compound_count(6)=e_lo%compound_count(5)*naky
    case ('lxyes')
       e_lo%ig_ord=1
       e_lo%isgn_ord=2
       e_lo%il_ord=3
       e_lo%it_ord=4
       e_lo%ik_ord=5
       e_lo%is_ord=6
       e_lo%compound_count(1)=1
       e_lo%compound_count(2)=e_lo%ntgridtotal
       e_lo%compound_count(3)=e_lo%compound_count(2)*2
       e_lo%compound_count(4)=e_lo%compound_count(3)*nlambda
       e_lo%compound_count(5)=e_lo%compound_count(4)*ntheta0
       e_lo%compound_count(6)=e_lo%compound_count(5)*naky
    case ('lyxes')
       e_lo%ig_ord=1
       e_lo%isgn_ord=2
       e_lo%il_ord=3
       e_lo%it_ord=5
       e_lo%ik_ord=4
       e_lo%is_ord=6
       e_lo%compound_count(1)=1
       e_lo%compound_count(2)=e_lo%ntgridtotal
       e_lo%compound_count(3)=e_lo%compound_count(2)*2
       e_lo%compound_count(4)=e_lo%compound_count(3)*nlambda
       e_lo%compound_count(5)=e_lo%compound_count(4)*naky
       e_lo%compound_count(6)=e_lo%compound_count(5)*ntheta0
    case ('xyles')
       e_lo%ig_ord=1
       e_lo%isgn_ord=2
       e_lo%il_ord=5
       e_lo%it_ord=3
       e_lo%ik_ord=4
       e_lo%is_ord=6
       e_lo%compound_count(1)=1
       e_lo%compound_count(2)=e_lo%ntgridtotal
       e_lo%compound_count(3)=e_lo%compound_count(2)*2
       e_lo%compound_count(4)=e_lo%compound_count(3)*ntheta0
       e_lo%compound_count(5)=e_lo%compound_count(4)*naky
       e_lo%compound_count(6)=e_lo%compound_count(5)*nlambda
    end select
    e_lo%ig_comp=e_lo%compound_count(e_lo%ig_ord)
    e_lo%isgn_comp=e_lo%compound_count(e_lo%isgn_ord)
    e_lo%ik_comp=e_lo%compound_count(e_lo%ik_ord)
    e_lo%it_comp=e_lo%compound_count(e_lo%it_ord)
    e_lo%il_comp=e_lo%compound_count(e_lo%il_ord)
    e_lo%is_comp=e_lo%compound_count(e_lo%is_ord)
!</DD>

! TT>
# ifdef USE_C_INDEX
    ierr = init_indices_elo_c (layout)
    if (ierr /= 0) &
         & write (error_unit(),*) 'ERROR: layout not found: ', trim(layout)
# endif
! <TT

  end subroutine init_energy_layouts

! TT>
# ifdef USE_C_INDEX
  function is_idx_e (lo, i)
# else
  elemental function is_idx_e (lo, i)
# endif
! <TT
    implicit none
    integer :: is_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function is_idx_e_c (lo,num)
         use layouts_type, only: e_layout_type
         integer :: is_idx_e_c
         type (e_layout_type) :: lo
         integer :: num
       end function is_idx_e_c
    end interface
    is_idx_e = is_idx_e_c (lo,i)
# else
    is_idx_e = 1 + mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
# endif
! <TT

  end function is_idx_e

! TT>
# ifdef USE_C_INDEX
  function il_idx_e (lo, i)
# else
  elemental function il_idx_e (lo, i)
# endif
! <TT
    implicit none
    integer :: il_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function il_idx_e_c (lo,num)
         use layouts_type, only: e_layout_type
         integer :: il_idx_e_c
         type (e_layout_type) :: lo
         integer :: num
       end function il_idx_e_c
    end interface
    il_idx_e = il_idx_e_c (lo,i)
# else
    il_idx_e = 1 + mod((i - lo%llim_world)/lo%il_comp, lo%nlambda)
# endif
! <TT

  end function il_idx_e


  elemental function isign_idx_e (lo, i)
    implicit none
    integer :: isign_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    isign_idx_e = 1 + mod((i - lo%llim_world)/(2*lo%ntgrid + 1),lo%nsign)
  end function isign_idx_e



! TT>
# ifdef USE_C_INDEX
  function it_idx_e (lo, i)
# else
  elemental function it_idx_e (lo, i)
# endif
! <TT
    implicit none
    integer :: it_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function it_idx_e_c (lo,num)
         use layouts_type, only: e_layout_type
         integer :: it_idx_e_c
         type (e_layout_type) :: lo
         integer :: num
       end function it_idx_e_c
    end interface
    it_idx_e = it_idx_e_c (lo,i)
# else
    it_idx_e = 1 + mod((i - lo%llim_world)/lo%it_comp, lo%ntheta0)
# endif
! <TT

  end function it_idx_e

! TT>
# ifdef USE_C_INDEX
  function ik_idx_e (lo, i)
# else
  elemental function ik_idx_e (lo, i)
# endif
! <TT
    implicit none
    integer :: ik_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function ik_idx_e_c (lo,num)
         use layouts_type, only: e_layout_type
         integer :: ik_idx_e_c
         type (e_layout_type) :: lo
         integer :: num
       end function ik_idx_e_c
    end interface
    ik_idx_e = ik_idx_e_c (lo,i)
# else
    ik_idx_e = 1 + mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
# endif
! <TT

  end function ik_idx_e

  elemental function ig_idx_e (lo, i)
    implicit none
    integer :: ig_idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_e = -lo%ntgrid + mod((i - lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_e

! TT>
# ifdef USE_C_INDEX
  function idx_e (lo, ig, isign, ik, it, il, is)
# else
  elemental function idx_e (lo, ig, isign, ik, it, il, is)
# endif
! <TT
    implicit none
    integer :: idx_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, isign, ik, it, il, is
! TT>
# ifdef USE_C_INDEX
    interface
       function idx_e_c (lo,ig,isign,ik,it,il,is)
         use layouts_type, only: e_layout_type
         integer :: idx_e_c
         type (e_layout_type) :: lo
         integer :: ig,isign,ik,it,il,is
       end function idx_e_c
    end interface
    idx_e = idx_e_c (lo, ig, isign, ik, it, il, is)
# else
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
    case ('xyles')
       idx_e = ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 + lo%nsign*(it-1 &
            + lo%ntheta0*(ik-1 + lo%naky*(il-1 + lo%nlambda*(is-1)))))
    end select
# endif
! <TT

  end function idx_e

  elemental function proc_id_e (lo, i)
    implicit none
    integer :: proc_id_e
    type (e_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    proc_id_e = i/lo%blocksize

  end function proc_id_e

! TT>
# ifdef USE_C_INDEX
  function idx_local_e (lo, ig, isign, ik, it, il, is)
# else
  elemental function idx_local_e (lo, ig, isign, ik, it, il, is)
# endif
! <TT
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
! Lambda layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_lambda_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2)
    use mp, only: iproc, nproc
! TT>
    use file_utils, only: error_unit
! <TT
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, ng2
    integer :: ngroup, nprocset
! TT>
# ifdef USE_C_INDEX
    integer :: ierr
    interface
       function init_indices_lzlo_c (layout)
         integer :: init_indices_lzlo_c
         character(*) :: layout
       end function init_indices_lzlo_c
    end interface
# endif
! <TT

!    if(proc0) write(6,*) 'Initializing lambda layouts', initialized_lambda_layouts
    if (initialized_lambda_layouts) return
    initialized_lambda_layouts = .true.

    lz_lo%iproc = iproc
    lz_lo%ntgrid = ntgrid
    lz_lo%naky = naky
    lz_lo%ntheta0 = ntheta0
    lz_lo%negrid = negrid
    lz_lo%nspec = nspec
    lz_lo%ntgridtotal=2*ntgrid+1
    lz_lo%ng2 = ng2
    lz_lo%llim_world = 0
    lz_lo%ulim_world = (2*ntgrid+1)*naky*ntheta0*negrid*nspec - 1

    lambda_local = .false.

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

!<DD>See g_lo init
    select case (layout)
    case ('yxels')
       lz_lo%ig_ord=1
       lz_lo%ie_ord=4
       lz_lo%it_ord=3
       lz_lo%ik_ord=2
       lz_lo%is_ord=5
       lz_lo%compound_count(1)=1
       lz_lo%compound_count(2)=lz_lo%ntgridtotal
       lz_lo%compound_count(3)=lz_lo%compound_count(2)*naky
       lz_lo%compound_count(4)=lz_lo%compound_count(3)*ntheta0
       lz_lo%compound_count(5)=lz_lo%compound_count(4)*negrid
    case ('yxles')
       lz_lo%ig_ord=1
       lz_lo%ie_ord=4
       lz_lo%it_ord=3
       lz_lo%ik_ord=2
       lz_lo%is_ord=5
       lz_lo%compound_count(1)=1
       lz_lo%compound_count(2)=lz_lo%ntgridtotal
       lz_lo%compound_count(3)=lz_lo%compound_count(2)*naky
       lz_lo%compound_count(4)=lz_lo%compound_count(3)*ntheta0
       lz_lo%compound_count(5)=lz_lo%compound_count(4)*negrid
    case ('lexys')
       lz_lo%ig_ord=1
       lz_lo%ie_ord=2
       lz_lo%it_ord=3
       lz_lo%ik_ord=4
       lz_lo%is_ord=5
       lz_lo%compound_count(1)=1
       lz_lo%compound_count(2)=lz_lo%ntgridtotal
       lz_lo%compound_count(3)=lz_lo%compound_count(2)*negrid
       lz_lo%compound_count(4)=lz_lo%compound_count(3)*ntheta0
       lz_lo%compound_count(5)=lz_lo%compound_count(4)*naky
    case ('lxyes')
       lz_lo%ig_ord=1
       lz_lo%ie_ord=4
       lz_lo%it_ord=2
       lz_lo%ik_ord=3
       lz_lo%is_ord=5
       lz_lo%compound_count(1)=1
       lz_lo%compound_count(2)=lz_lo%ntgridtotal
       lz_lo%compound_count(3)=lz_lo%compound_count(2)*ntheta0
       lz_lo%compound_count(4)=lz_lo%compound_count(3)*naky
       lz_lo%compound_count(5)=lz_lo%compound_count(4)*negrid
    case ('lyxes')
       lz_lo%ig_ord=1
       lz_lo%ie_ord=4
       lz_lo%it_ord=3
       lz_lo%ik_ord=2
       lz_lo%is_ord=5
       lz_lo%compound_count(1)=1
       lz_lo%compound_count(2)=lz_lo%ntgridtotal
       lz_lo%compound_count(3)=lz_lo%compound_count(2)*naky
       lz_lo%compound_count(4)=lz_lo%compound_count(3)*ntheta0
       lz_lo%compound_count(5)=lz_lo%compound_count(4)*negrid
    case ('xyles')
       lz_lo%ig_ord=1
       lz_lo%ie_ord=4
       lz_lo%it_ord=2
       lz_lo%ik_ord=3
       lz_lo%is_ord=5
       lz_lo%compound_count(1)=1
       lz_lo%compound_count(2)=lz_lo%ntgridtotal
       lz_lo%compound_count(3)=lz_lo%compound_count(2)*ntheta0
       lz_lo%compound_count(4)=lz_lo%compound_count(3)*naky
       lz_lo%compound_count(5)=lz_lo%compound_count(4)*negrid
    end select
    lz_lo%ig_comp=lz_lo%compound_count(lz_lo%ig_ord)
    lz_lo%ik_comp=lz_lo%compound_count(lz_lo%ik_ord)
    lz_lo%it_comp=lz_lo%compound_count(lz_lo%it_ord)
    lz_lo%ie_comp=lz_lo%compound_count(lz_lo%ie_ord)
    lz_lo%is_comp=lz_lo%compound_count(lz_lo%is_ord)
!</DD>

! TT>
# ifdef USE_C_INDEX
    ierr = init_indices_lzlo_c (layout)
    if (ierr /= 0) &
         & write (error_unit(),*) 'ERROR: layout not found: ', trim(layout)
# endif
! <TT

!    if(proc0) write(6,*) 'Initialized lambda layouts', initialized_lambda_layouts


  end subroutine init_lambda_layouts

! TT>
# ifdef USE_C_INDEX
  function is_idx_lz (lo, i)
# else
  elemental function is_idx_lz (lo, i)
# endif
! <TT
    implicit none
    integer :: is_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function is_idx_lz_c (lo,num)
         use layouts_type, only: lz_layout_type
         integer :: is_idx_lz_c
         type (lz_layout_type) :: lo
         integer :: num
       end function is_idx_lz_c
    end interface
    is_idx_lz = is_idx_lz_c (lo, i)
# else
    is_idx_lz = 1 + mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
# endif
! <TT

  end function is_idx_lz

! TT>
# ifdef USE_C_INDEX
  function ie_idx_lz (lo, i)
# else
  elemental function ie_idx_lz (lo, i)
# endif
! <TT
    implicit none
    integer :: ie_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i

! TT>
# ifdef USE_C_INDEX
    interface
       function ie_idx_lz_c (lo,num)
         use layouts_type, only: lz_layout_type
         integer :: ie_idx_lz_c
         type (lz_layout_type) :: lo
         integer :: num
       end function ie_idx_lz_c
    end interface
    ie_idx_lz = ie_idx_lz_c (lo, i)
# else
    ie_idx_lz = 1 + mod((i - lo%llim_world)/lo%ie_comp, lo%negrid)
# endif
! <TT

  end function ie_idx_lz

! TT>
# ifdef USE_C_INDEX
  function it_idx_lz (lo, i)
# else
  elemental function it_idx_lz (lo, i)
# endif
! <TT
    implicit none
    integer :: it_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function it_idx_lz_c (lo,num)
         use layouts_type, only: lz_layout_type
         integer :: it_idx_lz_c
         type (lz_layout_type) :: lo
         integer :: num
       end function it_idx_lz_c
    end interface
    it_idx_lz = it_idx_lz_c (lo, i)
# else
    it_idx_lz = 1 + mod((i - lo%llim_world)/lo%it_comp, lo%ntheta0)
# endif
! <TT

  end function it_idx_lz

! TT>
# ifdef USE_C_INDEX
  function ik_idx_lz (lo, i)
# else
  elemental function ik_idx_lz (lo, i)
# endif
! <TT
    implicit none
    integer :: ik_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
! TT>
# ifdef USE_C_INDEX
    interface
       function ik_idx_lz_c (lo,num)
         use layouts_type, only: lz_layout_type
         integer :: ik_idx_lz_c
         type (lz_layout_type) :: lo
         integer :: num
       end function ik_idx_lz_c
    end interface
    ik_idx_lz = ik_idx_lz_c (lo, i)
# else
    ik_idx_lz = 1 + mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
# endif
! <TT

  end function ik_idx_lz

  elemental function ig_idx_lz (lo, i)
    implicit none
    integer :: ig_idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_lz = -lo%ntgrid + mod((i - lo%llim_world)/lo%ig_comp, lo%ntgridtotal)
  end function ig_idx_lz

! TT>
# ifdef USE_C_INDEX
  function idx_lz (lo, ig, ik, it, ie, is)
# else
  elemental function idx_lz (lo, ig, ik, it, ie, is)
# endif
! <TT
    implicit none
    integer :: idx_lz
    type (lz_layout_type), intent (in) :: lo
    integer, intent (in) :: ig, ik, it, ie, is
! TT>
# ifdef USE_C_INDEX
    interface
       function idx_lz_c (lo,ig,ik,it,ie,is)
         use layouts_type, only: lz_layout_type
         integer :: idx_lz_c
         type (lz_layout_type) :: lo
         integer :: ig,ik,it,ie,is
       end function idx_lz_c
    end interface
    idx_lz = idx_lz_c (lo, ig, ik, it, ie, is)
# else
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
    case ('xyles')
       idx_lz = ig+lo%ntgrid + (2*lo%ntgrid+1)*(it-1 + lo%ntheta0*(ik-1 &
            + lo%naky*(ie-1 + lo%negrid*(is-1))))
    end select
# endif
! <TT

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

! TT>
# ifdef USE_C_INDEX
  function idx_local_lz (lo, ig, ik, it, ie, is)
# else
  elemental function idx_local_lz (lo, ig, ik, it, ie, is)
# endif
! <TT
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! layouts for parity diagnostic 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_parity_layouts &
       (naky, nlambda, negrid, nspec)
    use mp, only: iproc, nproc
    implicit none
    integer, intent (in) :: naky, nlambda, negrid, nspec

    if (initialized_parity_layouts) return
    initialized_parity_layouts = .true.
    
    p_lo%iproc = iproc
    p_lo%naky = naky
    p_lo%nlambda = nlambda
    p_lo%negrid = negrid
    p_lo%nspec = nspec
    p_lo%llim_world = 0
    p_lo%ulim_world = naky*negrid*nlambda*nspec - 1
    p_lo%blocksize = p_lo%ulim_world/nproc + 1
    p_lo%llim_proc = p_lo%blocksize*iproc
    p_lo%ulim_proc = min(p_lo%ulim_world, p_lo%llim_proc + p_lo%blocksize - 1)
    p_lo%ulim_alloc = max(p_lo%llim_proc, p_lo%ulim_proc)
!<DD>See g_lo init
    select case (layout)
    case ('yxels')
       p_lo%ie_ord=2
       p_lo%il_ord=3
       p_lo%ik_ord=1
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%naky
       p_lo%compound_count(3)=p_lo%compound_count(2)*negrid
       p_lo%compound_count(4)=p_lo%compound_count(3)*nlambda
    case ('yxles')
       p_lo%ie_ord=3
       p_lo%il_ord=2
       p_lo%ik_ord=1
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%naky
       p_lo%compound_count(3)=p_lo%compound_count(2)*nlambda
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    case ('lexys')
       p_lo%ie_ord=2
       p_lo%il_ord=1
       p_lo%ik_ord=3
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%nlambda
       p_lo%compound_count(3)=p_lo%compound_count(2)*negrid
       p_lo%compound_count(4)=p_lo%compound_count(3)*naky
    case ('lxyes')
       p_lo%ie_ord=3
       p_lo%il_ord=1
       p_lo%ik_ord=2
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%nlambda
       p_lo%compound_count(3)=p_lo%compound_count(2)*naky
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    case ('lyxes')
       p_lo%ie_ord=3
       p_lo%il_ord=1
       p_lo%ik_ord=2
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%nlambda
       p_lo%compound_count(3)=p_lo%compound_count(2)*naky
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    case ('xyles')
       p_lo%ie_ord=3
       p_lo%il_ord=2
       p_lo%ik_ord=1
       p_lo%is_ord=4
       p_lo%compound_count(1)=1
       p_lo%compound_count(2)=p_lo%naky
       p_lo%compound_count(3)=p_lo%compound_count(2)*nlambda
       p_lo%compound_count(4)=p_lo%compound_count(3)*negrid
    end select
    p_lo%ie_comp=p_lo%compound_count(p_lo%ie_ord)
    p_lo%ik_comp=p_lo%compound_count(p_lo%ik_ord)
    p_lo%il_comp=p_lo%compound_count(p_lo%il_ord)
    p_lo%is_comp=p_lo%compound_count(p_lo%is_ord)
!</DD>

  end subroutine init_parity_layouts

  elemental function idx_parity (lo, ik, il, ie, is)
    implicit none
    integer :: idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, il, ie, is
    select case (layout)
    case ('yxels')
       idx_parity = ik-1 + lo%naky*(ie-1 + lo%negrid*(il-1 + lo%nlambda*(is-1)))
    case ('yxles')
       idx_parity = ik-1 + lo%naky*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))
    case ('lexys')
       idx_parity = il-1 + lo%nlambda*(ie-1 + lo%negrid*(ik-1 + lo%naky*(is-1)))
    case ('lxyes')
       idx_parity = il-1 + lo%nlambda*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1)))
    case ('lyxes')
       idx_parity = il-1 + lo%nlambda*(ik-1 + lo%naky*(ie-1 + lo%negrid*(is-1)))
    case ('xyles')
       idx_parity = ik-1 + lo%naky*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))
    end select
  end function idx_parity

  elemental function idx_local_parity (lo, ik, il, ie, is)
    implicit none
    logical :: idx_local_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: ik, il, ie, is

    idx_local_parity = idx_local(lo, idx(lo, ik, il, ie, is))
  end function idx_local_parity

  elemental function ig_local_parity (lo, ip)
    implicit none
    logical :: ig_local_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: ip

    ig_local_parity = lo%iproc == proc_id(lo, ip)
  end function ig_local_parity

  elemental function proc_id_parity (lo, i)
    implicit none
    integer :: proc_id_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    proc_id_parity = i/lo%blocksize
  end function proc_id_parity

  elemental function ik_idx_parity (lo, i)
    implicit none
    integer :: ik_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_parity = 1+mod((i - lo%llim_world)/lo%ik_comp, lo%naky)
  end function ik_idx_parity

  elemental function is_idx_parity (lo, i)
    
    implicit none

    integer :: is_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    is_idx_parity = 1+mod((i - lo%llim_world)/lo%is_comp, lo%nspec)
  end function is_idx_parity

  elemental function il_idx_parity (lo, i)

    implicit none

    integer :: il_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_parity = 1+mod((i - lo%llim_world)/lo%il_comp, lo%nlambda)
  end function il_idx_parity

  elemental function ie_idx_parity (lo, i)

    implicit none

    integer :: ie_idx_parity
    type (p_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_parity = 1+mod((i - lo%llim_world)/lo%ie_comp, lo%negrid)
  end function ie_idx_parity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X-space layouts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_x_transform_layouts &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use mp, only: iproc, nproc, proc0, barrier
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    integer :: nprocset, ngroup, nblock, ntgridtotal, nsign
    real :: unbalanced_amount

    if (initialized_x_transform) return
    initialized_x_transform = .true.

    xxf_lo%iproc = iproc
    xxf_lo%ntgrid = ntgrid
    ntgridtotal=(2*ntgrid+1)
    xxf_lo%ntgridtotal = ntgridtotal
    nsign=2
    xxf_lo%nsign = nsign
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

!<DD>See g_lo init
    select case (layout)
    case ('yxels')
       xxf_lo%ig_ord=2
       xxf_lo%isgn_ord=3
       xxf_lo%il_ord=5
       xxf_lo%ie_ord=4
       xxf_lo%ik_ord=1
       xxf_lo%is_ord=6
       xxf_lo%compound_count(1)=1
       xxf_lo%compound_count(2)=xxf_lo%naky
       xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*ntgridtotal
       xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*nsign
       xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*negrid
       xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*nlambda
    case ('yxles')
       xxf_lo%ig_ord=2
       xxf_lo%isgn_ord=3
       xxf_lo%il_ord=4
       xxf_lo%ie_ord=5
       xxf_lo%ik_ord=1
       xxf_lo%is_ord=6
       xxf_lo%compound_count(1)=1
       xxf_lo%compound_count(2)=xxf_lo%naky
       xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*ntgridtotal
       xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*nsign
       xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*nlambda
       xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*negrid
    case ('lexys')
       xxf_lo%ig_ord=2
       xxf_lo%isgn_ord=3
       xxf_lo%il_ord=4
       xxf_lo%ie_ord=5
       xxf_lo%ik_ord=1
       xxf_lo%is_ord=6
       xxf_lo%compound_count(1)=1
       xxf_lo%compound_count(2)=xxf_lo%naky
       xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*ntgridtotal
       xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*nsign
       xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*nlambda
       xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*negrid
    case ('lxyes')
       xxf_lo%ig_ord=2
       xxf_lo%isgn_ord=3
       xxf_lo%il_ord=4
       xxf_lo%ie_ord=5
       xxf_lo%ik_ord=1
       xxf_lo%is_ord=6
       xxf_lo%compound_count(1)=1
       xxf_lo%compound_count(2)=xxf_lo%naky
       xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*ntgridtotal
       xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*nsign
       xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*nlambda
       xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*negrid
    case ('lyxes')
       xxf_lo%ig_ord=2
       xxf_lo%isgn_ord=3
       xxf_lo%il_ord=4
       xxf_lo%ie_ord=5
       xxf_lo%ik_ord=1
       xxf_lo%is_ord=6
       xxf_lo%compound_count(1)=1
       xxf_lo%compound_count(2)=xxf_lo%naky
       xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*ntgridtotal
       xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*nsign
       xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*nlambda
       xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*negrid
    case ('xyles')
       xxf_lo%ig_ord=2
       xxf_lo%isgn_ord=3
       xxf_lo%il_ord=4
       xxf_lo%ie_ord=5
       xxf_lo%ik_ord=1
       xxf_lo%is_ord=6
       xxf_lo%compound_count(1)=1
       xxf_lo%compound_count(2)=xxf_lo%naky
       xxf_lo%compound_count(3)=xxf_lo%compound_count(2)*ntgridtotal
       xxf_lo%compound_count(4)=xxf_lo%compound_count(3)*nsign
       xxf_lo%compound_count(5)=xxf_lo%compound_count(4)*nlambda
       xxf_lo%compound_count(6)=xxf_lo%compound_count(5)*negrid
    end select
    xxf_lo%ig_comp=xxf_lo%compound_count(xxf_lo%ig_ord)
    xxf_lo%isgn_comp=xxf_lo%compound_count(xxf_lo%isgn_ord)
    xxf_lo%ik_comp=xxf_lo%compound_count(xxf_lo%ik_ord)
    xxf_lo%ie_comp=xxf_lo%compound_count(xxf_lo%ie_ord)
    xxf_lo%il_comp=xxf_lo%compound_count(xxf_lo%il_ord)
    xxf_lo%is_comp=xxf_lo%compound_count(xxf_lo%is_ord)
!</DD>

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
       ! naky*(2*ntgrid+1)*nsign*nlambda*negrid*nspec
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
             write(*, fmt="('Using unbalanced decomposition for xxf. ',&
              & 'Unbalanced fraction',F6.2)") unbalanced_amount
          end if

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
       
       ! For the xxf decomposition the next factor to be used
       ! depends on the layout chosen.  If it is the yxels
       ! decompsoition then the next factor is nlambda, otherwise
       ! it is negrid.
       ! Multiple the previous factor by the next factor to be
       ! be used.  If the previous factor divided exactly by the
       ! processes available then k will be 1 so we will only be
       ! considering this factor.
       select case(layout)
       case ('yxels')
          k = xxf_lo%nlambda * k
       case default
          k = xxf_lo%negrid * k
       end select
       
       ! The next code follows the pattern described above
       ! moving through all the other factor, i.e.:
       ! for yxels: negrid*nsign*ntgridtotal*naky
       ! and for the other layouts: nlambda*nsign*ntgridtotal*naky
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
          
          select case(layout)
          case ('yxels')
             k = xxf_lo%negrid * k
          case default
             k = xxf_lo%nlambda * k
          end select
          
          call calculate_block_breakdown(k, i, m, j, level_proc_num)
          
          if(j .eq. 0) then
             
             if(m .eq. 0) then
                level_proc_num = level_proc_num/k
                k = 1
             end if
             
             k = xxf_lo%nsign * k
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
                      call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                           1, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, &
                           xxf_lo%ulim_alloc)
                      
                   end if
                   
                else
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)     
                      ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                      ! naky for this level of the decomposition.
                      call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                           xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple,  xxf_lo%llim_proc, xxf_lo%ulim_proc, &
                           xxf_lo%ulim_alloc)
                      
                   end if
                   
                end if
                
             else
                
                if(i .ne. 0) then
                   
                   call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)                         
                   ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                   ! naky,ntgridtotal for this level of the decomposition.
                   call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                        xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, &
                        xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
                   
                end if
                
             end if
             
             
          else
             
             if(i .ne. 0) then
                
                call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nsign for this level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                     xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, xxf_lo%block_multiple, &
                     xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
                
             end if
             
          end if
          
       else
          
          if(i .ne. 0) then
             
             call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
             
             select case(layout)
             case('yxels')
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nsign,negrid for this layout and level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                     xxf_lo%negrid*xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, &
                     xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
             case default
                ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
                ! naky,ntgridtotal,nsign,nlambda for these layouts and this level of the decomposition.
                call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
                     xxf_lo%nlambda*xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, xxf_lo%large_block_size, &
                     xxf_lo%block_multiple, xxf_lo%llim_proc, xxf_lo%ulim_proc, xxf_lo%ulim_alloc)
             end select
             
          end if
          
       end if
       
    else
       
       if(i .ne. 0) then
          
          call calculate_unbalanced_decomposition(k, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, xxf_lo%num_small, xxf_lo%num_large, level_proc_num)
          ! Calculate the block sizes using the factor of the decomposition that has not yet been split up, namely
          ! naky,ntgridtotal,nsign,nlambda,negrid for this level of the decomposition.
          call calculate_block_size(iproc, xxf_lo%num_small, xxf_lo%num_large, xxf_lo%small_block_balance_factor, xxf_lo%large_block_balance_factor, &
               xxf_lo%negrid*xxf_lo%nlambda*xxf_lo%nsign*xxf_lo%ntgridtotal*xxf_lo%naky, xxf_lo%blocksize, xxf_lo%small_block_size, &
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


  subroutine calculate_block_size(iproc, numsmall, numlarge, smalldecomp, largedecomp, sizeblock, blocksize, smallblocksize, largeblocksize, block_multiple, llim, ulim, ulim_alloc)
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

    integer, intent(in) :: iproc, numsmall, numlarge, smalldecomp, largedecomp, sizeblock
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
!AJ       llim = llim + smallblocksize * (modproc - 1)
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
    is_idx_xxf = 1+mod((i-lo%llim_world)/lo%is_comp,lo%nspec)
  end function is_idx_xxf

  elemental function ie_idx_xxf (lo, i)
    implicit none
    integer :: ie_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_xxf=1+mod((i-lo%llim_world)/lo%ie_comp,lo%negrid)
  end function ie_idx_xxf

  elemental function il_idx_xxf (lo, i)
    implicit none
    integer :: il_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_xxf = 1+mod((i-lo%llim_world)/lo%il_comp,lo%nlambda)
  end function il_idx_xxf

  elemental function isign_idx_xxf (lo, i)
    implicit none
    integer :: isign_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    isign_idx_xxf = 1+mod((i-lo%llim_world)/lo%isgn_comp,lo%nsign)
  end function isign_idx_xxf

  elemental function ig_idx_xxf (lo, i)
    implicit none
    integer :: ig_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_xxf = -lo%ntgrid+mod((i-lo%llim_world)/lo%ig_comp,lo%ntgridtotal)
  end function ig_idx_xxf

  elemental function ik_idx_xxf (lo, i)
    implicit none
    integer :: ik_idx_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_xxf = 1+mod((i-lo%llim_world)/lo%ik_comp,lo%naky)
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
    case ('xyles')
       idx_xxf = ik-1 + lo%naky*(ig+lo%ntgrid + (2*lo%ntgrid+1)*(isign-1 &
            + lo%nsign*(il-1 + lo%nlambda*(ie-1 + lo%negrid*(is-1)))))
    end select
  end function idx_xxf


  elemental function proc_id_xxf (lo, i)
    implicit none
    integer :: proc_id_xxf
    type (xxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    integer :: block_offset, offset_block_number, j, k, tempi
    
    if (accel_lxyes) then
       proc_id_xxf = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
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
    use mp, only: iproc, nproc, proc0
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    integer :: nnx, nny, ngroup, nprocset, nblock, ntgridtotal, nsign
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
    ntgridtotal=(2*ntgrid+1)
    yxf_lo%ntgridtotal = ntgridtotal
    nsign=2
    yxf_lo%nsign = nsign
    yxf_lo%naky = naky
    yxf_lo%ny = nny
    yxf_lo%ntheta0 = ntheta0
    yxf_lo%nx = nnx
    yxf_lo%nlambda = nlambda
    yxf_lo%negrid = negrid
    yxf_lo%nspec = nspec
    yxf_lo%llim_world = 0
    yxf_lo%ulim_world = nnx*(2*ntgrid+1)*2*nlambda*negrid*nspec - 1

!<DD>Calculate compound index divisors
    select case (layout)
    case ('yxels')
       yxf_lo%ig_ord=2
       yxf_lo%isgn_ord=3
       yxf_lo%il_ord=5
       yxf_lo%ie_ord=4
       yxf_lo%it_ord=1
       yxf_lo%is_ord=6
       yxf_lo%compound_count(1)=1
       yxf_lo%compound_count(2)=yxf_lo%nx
       yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*ntgridtotal
       yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*nsign
       yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*negrid
       yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*nlambda
    case ('yxles')
       yxf_lo%ig_ord=2
       yxf_lo%isgn_ord=3
       yxf_lo%il_ord=4
       yxf_lo%ie_ord=5
       yxf_lo%it_ord=1
       yxf_lo%is_ord=6
       yxf_lo%compound_count(1)=1
       yxf_lo%compound_count(2)=yxf_lo%nx
       yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*ntgridtotal
       yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*nsign
       yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nlambda
       yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*negrid
    case ('lexys')
       yxf_lo%ig_ord=2
       yxf_lo%isgn_ord=3
       yxf_lo%il_ord=4
       yxf_lo%ie_ord=5
       yxf_lo%it_ord=1
       yxf_lo%is_ord=6
       yxf_lo%compound_count(1)=1
       yxf_lo%compound_count(2)=yxf_lo%nx
       yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*ntgridtotal
       yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*nsign
       yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nlambda
       yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*negrid
    case ('lxyes')
       yxf_lo%ig_ord=2
       yxf_lo%isgn_ord=3
       yxf_lo%il_ord=4
       yxf_lo%ie_ord=5
       yxf_lo%it_ord=1
       yxf_lo%is_ord=6
       yxf_lo%compound_count(1)=1
       yxf_lo%compound_count(2)=yxf_lo%nx
       yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*ntgridtotal
       yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*nsign
       yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nlambda
       yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*negrid
    case ('lyxes')
       yxf_lo%ig_ord=2
       yxf_lo%isgn_ord=3
       yxf_lo%il_ord=4
       yxf_lo%ie_ord=5
       yxf_lo%it_ord=1
       yxf_lo%is_ord=6
       yxf_lo%compound_count(1)=1
       yxf_lo%compound_count(2)=yxf_lo%nx
       yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*ntgridtotal
       yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*nsign
       yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nlambda
       yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*negrid
    case ('xyles')
       yxf_lo%ig_ord=2
       yxf_lo%isgn_ord=3
       yxf_lo%il_ord=4
       yxf_lo%ie_ord=5
       yxf_lo%it_ord=1
       yxf_lo%is_ord=6
       yxf_lo%compound_count(1)=1
       yxf_lo%compound_count(2)=yxf_lo%nx
       yxf_lo%compound_count(3)=yxf_lo%compound_count(2)*ntgridtotal
       yxf_lo%compound_count(4)=yxf_lo%compound_count(3)*nsign
       yxf_lo%compound_count(5)=yxf_lo%compound_count(4)*nlambda
       yxf_lo%compound_count(6)=yxf_lo%compound_count(5)*negrid
    end select
    yxf_lo%ig_comp=yxf_lo%compound_count(yxf_lo%ig_ord)
    yxf_lo%isgn_comp=yxf_lo%compound_count(yxf_lo%isgn_ord)
    yxf_lo%it_comp=yxf_lo%compound_count(yxf_lo%it_ord)
    yxf_lo%ie_comp=yxf_lo%compound_count(yxf_lo%ie_ord)
    yxf_lo%il_comp=yxf_lo%compound_count(yxf_lo%il_ord)
    yxf_lo%is_comp=yxf_lo%compound_count(yxf_lo%is_ord)
!</DD>

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
       ! nxx*(2*ntgrid+1)*nsign*nlambda*negrid*nspec
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
             write(*, fmt="('Using unbalanced decomposition for yxf. ',&
             & 'Unbalanced fraction',F6.2)") unbalanced_amount
          end if

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
       
       select case(layout)
       case ('yxels')
          k = yxf_lo%nlambda * k
       case default
          k = yxf_lo%negrid * k
       end select
       call calculate_block_breakdown(k, i, m, j, level_proc_num)
       
       if(j .eq. 0) then
          
          if(m .eq. 0) then
             level_proc_num = level_proc_num/k
             k = 1
          end if
          
          select case(layout)
          case ('yxels')
             k = yxf_lo%negrid * k
          case default
             k = yxf_lo%nlambda * k
          end select
          
          call calculate_block_breakdown(k, i, m, j, level_proc_num)
          
          if(j .eq. 0) then
             
             if(m .eq. 0) then
                level_proc_num = level_proc_num/k
                k = 1
             end if
             
             k = yxf_lo%nsign * k
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
                      call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                           1, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                      
                   end if
                   
                else
                   
                   if(i .ne. 0) then
                      
                      call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)     
                      call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                           yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                      
                   end if
                   
                end if
                
             else
                
                if(i .ne. 0) then
                   
                   call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)                         
                   call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                        yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, &
                        yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                   
                end if
                
             end if
             
             
          else
             
             if(i .ne. 0) then
                
                call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                     yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, yxf_lo%block_multiple, &
                     yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
                
             end if
             
          end if
          
       else
          
          if(i .ne. 0) then
             
             call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
             
             select case(layout)
             case('yxels')
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                     yxf_lo%negrid*yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, &
                     yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
             case default
                call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
                     yxf_lo%nlambda*yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, yxf_lo%large_block_size, &
                     yxf_lo%block_multiple, yxf_lo%llim_proc, yxf_lo%ulim_proc, yxf_lo%ulim_alloc)
             end select
             
          end if
          
       end if
       
    else
       
       if(i .ne. 0) then
          
          call calculate_unbalanced_decomposition(k, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, yxf_lo%num_small, yxf_lo%num_large, level_proc_num)
          call calculate_block_size(iproc, yxf_lo%num_small, yxf_lo%num_large, yxf_lo%small_block_balance_factor, yxf_lo%large_block_balance_factor, &
               yxf_lo%negrid*yxf_lo%nlambda*yxf_lo%nsign*yxf_lo%ntgridtotal*yxf_lo%nx, yxf_lo%blocksize, yxf_lo%small_block_size, &
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
    is_idx_yxf = 1+mod((i-lo%llim_world)/lo%is_comp,lo%nspec)
  end function is_idx_yxf

  elemental function ie_idx_yxf (lo, i)
    implicit none
    integer :: ie_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ie_idx_yxf = 1+mod((i-lo%llim_world)/lo%ie_comp,lo%negrid)
  end function ie_idx_yxf

  elemental function il_idx_yxf (lo, i)
    implicit none
    integer :: il_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    il_idx_yxf = 1+mod((i-lo%llim_world)/lo%il_comp,lo%nlambda)
  end function il_idx_yxf

  elemental function isign_idx_yxf (lo, i)
    implicit none
    integer :: isign_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    isign_idx_yxf = 1+mod((i-lo%llim_world)/lo%isgn_comp,lo%nsign)
 end function isign_idx_yxf

  elemental function ig_idx_yxf (lo, i)
    implicit none
    integer :: ig_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ig_idx_yxf = -lo%ntgrid+mod((i-lo%llim_world)/lo%ig_comp,lo%ntgridtotal)
  end function ig_idx_yxf

  elemental function it_idx_yxf (lo, i)
    implicit none
    integer :: it_idx_yxf
    type (yxf_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_yxf = 1+mod((i-lo%llim_world)/lo%it_comp,lo%nx)
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
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       idx_accelx = it-1 + lo%nx*(ik-1 + lo%ny*(il-1 + &
!            lo%nlambda*(ie-1 + lo%negrid*(is-1))))
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
    case ('xyles')
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
    integer :: block_offset, offset_block_number, j, k, tempi

    if (accel_lxyes) then
       proc_id_yxf = (i/lo%gsize)*lo%nprocset + mod(i, lo%gsize)/lo%nset
    else
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
    integer :: nnx, nny

    if (initialized_accel_transform_layouts) return
    initialized_accel_transform_layouts = .true.

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
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       ie_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny/lo%nlambda, lo%negrid)
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
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!   case ('xyles')
!       il_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nxny, lo%nlambda)
    end select
  end function il_idx_accelx

  elemental function it_idx_accelx (lo, i)
    implicit none
    integer :: it_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
!    select case (layout)
!    case ('yxels')
!       it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
!    case ('yxles')
!       it_idx_accelx = 1 + mod((i - lo%llim_world)/lo%ny, lo%nx)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       it_idx_accelx = 1 + mod((i - lo%llim_world), lo%nx)
!    end select
  end function it_idx_accelx

  elemental function ik_idx_accelx (lo, i)
    implicit none
    integer :: ik_idx_accelx
    type (accelx_layout_type), intent (in) :: lo
    integer, intent (in) :: i

    ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
!    select case (layout)
!    case ('yxels')
!       ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
!    case ('yxles')
!       ik_idx_accelx = 1 + mod((i - lo%llim_world), lo%ny)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       ik_idx_accelx = 1 + mod((i - lo%llim_world)/lo%nx, lo%ny)
!    end select
  end function ik_idx_accelx

  elemental function it_idx_accel (lo, i)
    implicit none
    integer :: it_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
!    select case (layout)
!    case ('yxels')
!       it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
!    case ('yxles')
!       it_idx_accel = 1 + mod((i - lo%llim_world)/lo%ndky, lo%nx)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       it_idx_accel = 1 + mod((i - lo%llim_world), lo%nx)
!    end select
  end function it_idx_accel

  elemental function ik_idx_accel (lo, i)
    implicit none
    integer :: ik_idx_accel
    type (accel_layout_type), intent (in) :: lo
    integer, intent (in) :: i
    ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
!    select case (layout)
!    case ('yxels')
!       ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
!    case ('yxles')
!       ik_idx_accel = 1 + mod(i - lo%llim_world, lo%ndky)
! CMR: 'xyles' possible if 2D accel FFTs handled x faster than y 
!    case ('xyles')
!       ik_idx_accel = 1 + mod((i - lo%llim_world)/lo%nx, lo%ndky)
!    end select
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


! TT>
# ifdef USE_C_INDEX
  subroutine gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, il, ilz)
# else
  pure subroutine gidx2lzidx (ig, isign, g_lo, iglo, lz_lo, ntgrid, jend, &
                              il, ilz)
# endif
! <TT
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

! TT>
# ifdef USE_C_INDEX
  subroutine gidx2gintidx (g_lo, iglo, gint_lo, il, igint)
# else
  elemental subroutine gidx2gintidx (g_lo, iglo, gint_lo, il, igint)
# endif
! <TT
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

! TT>
# ifdef USE_C_INDEX
  subroutine gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
# else
  elemental subroutine gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
# endif
! <TT
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

! TT>
# ifdef USE_C_INDEX
  subroutine xxfidx2gidx (it, ixxf, xxf_lo, g_lo, ig, isign, iglo)
# else
  elemental subroutine xxfidx2gidx (it, ixxf, xxf_lo, g_lo, ig, isign, iglo)
# endif
! <TT
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
    case ('xyles')
!CMR: set char="?"
! char only acts if set to "v", to force use of accel layout.
! accel layouts NOT implemented for "xyles", which assume throughout that 
! y is fastest index, including in calls to 2D FFTs
       char = '?'    
    end select

  end subroutine pe_layout

  subroutine factors (n, j, div)
    !CMR, 3/10/13:
    !    Find all the factors of n and return in div(j)
    integer, intent (in) :: n
    integer, intent (out) :: j
    integer, dimension (:), intent (out) :: div
    integer :: i, imax

    ! find: i = lowest factor of n
    ! and therefore imax=n/1 is the HIGHEST factor of n
    do i=2,n
       if (mod(n,i)==0) exit
    end do
    imax = n/i
    j=1
    ! loop over all possible factors of n, and return in div(j)
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
