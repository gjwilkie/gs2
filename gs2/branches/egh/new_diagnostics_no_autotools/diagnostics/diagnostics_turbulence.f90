!> A module for calculating properties of the turbulence
!! such as cross phases and correlations
module diagnostics_turbulence
  use diagnostics_config
  implicit none
  contains
    subroutine write_cross_phase(gnostics)
      use mp, only: proc0
      implicit none
      type(diagnostics_type), intent(in) :: gnostics
      real :: phase_tot, phase_theta
      real :: t

      t = gnostics%user_time

      call get_cross_phase (phase_tot, phase_theta)
      if (proc0) write (unit=phase_unit, fmt="('t= ',e17.10,' phase_tot= ',e11.4,' phase_theta= ',e11.4)") &
           & t, phase_tot, phase_theta
    end subroutine write_cross_phase

    !> This is a highly simplified synthetic diagnostic which 
    !! calculates the cross phase between the electron density and the 
    !! perpendicular electron temperature for comparisons with DIII-D.  
    !! Returns the value of the cross-phase at the outboard midplane and 
    !! integrated over all v and x. We can generalize this routine to 
    !! other fields at some point, but for now this is just a skeleton for 
    !! a more realistic synthetic diagnostic. 
    subroutine get_cross_phase (phase_tot, phase_theta)
      use species, only: nspec, spec, electron_species
      use kt_grids, only: ntheta0, naky
      use theta_grid, only: ntgrid
      use dist_fn, only: getemoms
      use fields_arrays, only: phinew, bparnew
      implicit none
      real, intent (out) :: phase_tot, phase_theta
      complex, dimension (:,:,:,:), allocatable :: ntot, tperp
      complex, dimension (ntheta0, naky) :: nTp_by_mode
      complex :: nTp
  !    real, dimension (ntheta0, naky) :: n2_by_mode, T2_by_mode
  !    real :: n2, T2

      integer ::  is

      allocate ( ntot(-ntgrid:ntgrid,ntheta0,naky,nspec))
      allocate (tperp(-ntgrid:ntgrid,ntheta0,naky,nspec))

      call getemoms (phinew, bparnew, ntot, tperp)

      do is = 1,nspec
         if (spec(is)%type == electron_species) then
            ! get cross_phase at outboard midplane
            call get_vol_int (ntot(0,:,:,is), tperp(0,:,:,is), nTp, nTp_by_mode)
  !          call get_vol_average (ntot(0,:,:,is), ntot(0,:,:,is), n2, n2_by_mode)
  !          call get_vol_average (tperp(0,:,:,is), tperp(0,:,:,is), T2, T2_by_mode)
            phase_theta = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
            ! get integrated cross_phase 
            call get_vol_int (ntot(:,:,:,is), tperp(:,:,:,is), nTp, nTp_by_mode)
  !          call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), n2, n2_by_mode)
  !          call get_vol_average (tperp(:,:,:,is), tperp(:,:,:,is), T2, T2_by_mode)
            phase_tot = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
         end if
      end do

      deallocate (ntot, tperp)

    end subroutine get_cross_phase
end module diagnostics_turbulence
