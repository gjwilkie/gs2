module gs2_save

  implicit none

  public :: gs2_restore, gs2_save_for_restart
  public :: init_save, init_dt, init_tstart

  interface gs2_restore
     module procedure gs2_restore_many, gs2_restore_one
  end interface
  
  private
  character(300), save :: restart_file

contains

  subroutine gs2_save_for_restart (g, t0, delt0, istatus, exit_in)

    use theta_grid, only: ntgrid
! Must include g_layout_type here to avoid obscure bomb while compiling
! gs2_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
    use gs2_layouts, only: g_lo, g_layout_type
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, intent (in) :: t0, delt0
    integer, intent (out) :: istatus
    logical, intent (in), optional :: exit_in

  end subroutine gs2_save_for_restart


  subroutine gs2_restore_many (g, scale, istatus, many)

    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
    logical, intent (in) :: many

  end subroutine gs2_restore_many

  subroutine gs2_restore_one (g, scale, istatus)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo

    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus

  end subroutine gs2_restore_one

  subroutine init_save (file)
    character(300), intent (in) :: file
    
    restart_file = file

  end subroutine init_save

  subroutine init_dt (delt0, istatus)

    implicit none
    real, intent (in out) :: delt0
    integer, intent (out) :: istatus

  end subroutine init_dt

  subroutine init_tstart (tstart, istatus)

    implicit none
    real, intent (in out) :: tstart
    integer, intent (out) :: istatus
    
  end subroutine init_tstart

end module gs2_save
