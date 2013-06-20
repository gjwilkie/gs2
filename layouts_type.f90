module layouts_type
  ! This can be made by just replacing nakx by ntheta0
  !   from AstroGK's layouts_type.f90

  implicit none

  type :: g_layout_type
     sequence
     integer :: iproc
     integer :: ntgrid, naky, ntheta0, nvgrid, nmu, nspec
     integer :: llim_world, ulim_world, llim_proc, ulim_proc, ulim_alloc, blocksize
  end type g_layout_type

end module layouts_type
