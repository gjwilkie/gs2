module layouts_type
  ! This can be made by just replacing nakx by ntheta0
  !   from AstroGK's layouts_type.f90

  implicit none

  type :: g_layout_type
     sequence
     integer :: iproc
     integer :: naky, ntheta0, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type g_layout_type

  type :: lz_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, negrid, nspec, ng2
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize, gsize
     integer :: llim_group, ulim_group, igroup, ngroup, nprocset, iset, nset, groupblocksize
  end type lz_layout_type

  type :: e_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nlambda, nspec, nsign
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type e_layout_type

  type :: le_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type le_layout_type

  type :: p_layout_type
     sequence
     integer :: iproc
     integer :: naky, nlambda, negrid, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type p_layout_type

end module layouts_type
