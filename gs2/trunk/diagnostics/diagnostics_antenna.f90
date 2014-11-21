!> A module to calculate the properties of the
!! antenna used to drive turbulence 
module diagnostics_antenna
  implicit none

  private

  public :: init_diagnostics_antenna, finish_diagnostics_antenna, write_jext
  public :: write_lorentzian

  real, dimension(:,:,:), allocatable ::  j_ext_hist

contains
  subroutine init_diagnostics_antenna(gnostics)
    use diagnostics_config, only: diagnostics_type
    use kt_grids, only: ntheta0, naky
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    if (gnostics%write_jext) allocate (j_ext_hist(ntheta0, naky,0:gnostics%navg-1)) 
  end subroutine init_diagnostics_antenna

  subroutine finish_diagnostics_antenna(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    if (gnostics%write_jext.and.allocated(j_ext_hist)) deallocate(j_ext_hist)
  end subroutine finish_diagnostics_antenna

  subroutine write_jext(gnostics)
    use kt_grids, only: ntheta0, naky
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_config, only: diagnostics_type
    use mp, only: proc0
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real:: t
    integer:: ik, it

    !GGH J_external
    real, dimension(:,:), allocatable ::  j_ext

    allocate (j_ext(ntheta0, naky)); j_ext=0.

    !Get z-centered j_ext at current time
    call calc_jext(gnostics,j_ext)

    !Write to netcdf
    call create_and_write_variable(gnostics, gnostics%rtype, "antenna_j_ext", "XYt", &
         "Time averaged external current in the antenna, real(kperp^2 A_antenna), &
         &  as a function of kx and ky", "radians", j_ext)

    !Write to ascii
    if (proc0 .and. gnostics%write_ascii) then 
       t = gnostics%user_time
       do ik=1,naky
          do it = 1, ntheta0
             if (j_ext(it,ik) .ne. 0.) then
                write (unit=gnostics%ascii_files%jext, fmt="(es12.4,i4,i4,es12.4)")  &
                     t,it,ik,j_ext(it,ik)
             endif
          enddo
       enddo
    end if

    deallocate(j_ext)
  end subroutine write_jext

  subroutine write_lorentzian(gnostics)
    use antenna, only: antenna_w
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: tmp2
    tmp2 = real(antenna_w())
    call create_and_write_variable(gnostics, gnostics%rtype, "antenna_w", "rt", &
         "antenna_w?? ", "TBC", tmp2)
  end subroutine write_lorentzian

  !> A subroutine to calculate the time-averaged antenna current
  !! j_ext = kperp^2 A_antenna
  subroutine calc_jext (gnostics, j_ext)
    use mp, only: proc0
    use dist_fn, only: get_jext
    use antenna, only: antenna_apar
    use volume_averages, only: average_theta
    use theta_grid, only: ntgrid
    use kt_grids, only: kperp2, ntheta0, naky
    use diagnostics_config, only: diagnostics_type
    implicit none
    !Passed
    type(diagnostics_type), intent(in) :: gnostics
    integer:: istep
    real, dimension(:,:), intent(out) ::  j_ext
    complex, dimension(:,:,:), allocatable :: j_extz
    !Local 
    integer :: i

    istep = gnostics%istep

    !call get_jext(j_ext)    
    allocate (j_extz(-ntgrid:ntgrid, ntheta0, naky)) ; j_extz = 0.
    call antenna_apar (kperp2, j_extz)       
    call average_theta(real(j_extz), j_ext, gnostics%distributed)
    deallocate (j_extz)

    !Do averages with a history variable
    if (proc0) then
       !Save variable to history
       if (gnostics%navg > 1) then
          if (istep > 1) j_ext_hist(:,:,mod(istep,gnostics%navg)) = j_ext(:,:)
          
          !Use average of history
          if (istep >= gnostics%navg) then
             j_ext=0.
             do i=0,gnostics%navg-1
                j_ext(:,:) = j_ext(:,:) + j_ext_hist(:,:,i) / real(gnostics%navg)
             end do
          end if
       end if
    end if
  end subroutine calc_jext
end module diagnostics_antenna
