module layouts_type
  ! This can be made by just replacing nakx by ntheta0
  !   from AstroGK's layouts_type.f90

  implicit none

  type :: g_layout_type
     sequence
     integer :: iproc
     integer :: naky, ntheta0, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ik_ord, it_ord, il_ord, ie_ord, is_ord
     integer :: ik_comp, it_comp, il_comp, ie_comp, is_comp
     integer,dimension(5) :: compound_count, dim_size
!    For replacing loops over iglo with 5 separate loops
!    This is not currently implemented, but involves replacing "do iglo=...." with
!    EQUIVALENCE(dind(1),d1),(dind(2),d2).... !Note this is just a short hand for storing currently d values in an array   
!    ik=>dind(g_lo%ik_ord) !This is a pointer assignment, so that whatever layout we use ik points to the correct d value
!    do d5=d5_min,d5max ; do d4=d4_min,d4_max ...
     integer :: d1_min,d1_max,d2_min,d2_max,d3_min,d3_max
     integer :: d4_min,d4_max,d5_min,d5_max
!    Used in deciding colour for sub-communicator creation.
!    Also used to select slices of array to take part in sub-comm global reductions 
!    (in order to reduce message size by avoiding regions of arrays outside this procs domain)
     integer :: ik_min, ik_max, it_min, it_max, il_min, il_max, ie_min, ie_max, is_min, is_max
     integer :: xyblock_comm !Sub comms for xy blocks (i.e. l-e-s integrals)
     integer :: xysblock_comm !Sub comms for xys blocks (i.e. l-e integrals)
     integer :: lesblock_comm !Sub comms for les blocks (i.e. x-y integrals)
     integer,dimension(:,:),allocatable :: les_kxky_range !Used in integrate species, holds start and stop indices for flattened kxky dimension
     logical :: x_before_y !Information about layout order
     logical :: x_local, y_local, l_local, e_local, s_local !Used to record if various dimensions are entirely local
! shm data
     integer :: ppn, nnd, bpn, ndbklsize, llim_node, ulim_node
     integer, allocatable :: proc_map(:) ! llimit on all ranks
  end type g_layout_type

  type :: lz_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
!<DD>
!    For branchless index lookup routines
     integer :: ig_ord, isgn_ord, ik_ord, it_ord, ie_ord, is_ord
     integer :: ig_comp, isgn_comp, ik_comp, it_comp, ie_comp, is_comp
     integer,dimension(6) :: compound_count

  end type lz_layout_type

  type :: e_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nlambda, nspec, nsign, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ig_ord, isgn_ord, ik_ord, it_ord, il_ord, is_ord
     integer :: ig_comp, isgn_comp, ik_comp, it_comp, il_comp, is_comp
     integer,dimension(6) :: compound_count
  end type e_layout_type

  type :: le_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nspec, ntgridtotal
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ig_ord, ik_ord, it_ord, is_ord
     integer :: ig_comp, ik_comp, it_comp, is_comp
     integer,dimension(4) :: compound_count
     integer :: ik_min, it_min, ig_min, is_min 
     integer :: ik_max, it_max, ig_max, is_max 
     logical :: x_local, y_local, t_local, s_local !Used to record if various dimensions are entirely local
  end type le_layout_type

  type :: p_layout_type
     sequence
     integer :: iproc
     integer :: naky, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
!<DD>
!    For branchless index lookup routines
     integer :: ik_ord, il_ord, ie_ord, is_ord
     integer :: ik_comp, il_comp, ie_comp, is_comp
     integer,dimension(4) :: compound_count
  end type p_layout_type

end module layouts_type
