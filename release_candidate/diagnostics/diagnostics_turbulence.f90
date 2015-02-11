!> A module for calculating properties of the turbulence
!! such as cross phases and correlations
module diagnostics_turbulence
  implicit none
  private

  public :: write_correlation_extend, write_cross_phase
  public :: write_correlation

  integer :: ntg_out, ntg_extend, nth0_extend
contains

  subroutine write_correlation(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky
    !use mp, only: proc0
    use gs2_io, only: nc_loop_corr
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    complex, dimension(:,:), allocatable :: phi_corr_2pi
    allocate (phi_corr_2pi(-ntgrid:ntgrid,naky))
    call correlation (phi_corr_2pi)
    call create_and_write_variable(gnostics, gnostics%rtype, "phi_corr_2pi", &
         dim_string([gnostics%dims%ri,gnostics%dims%theta,gnostics%dims%ky,gnostics%dims%time]), &
         "2 point correlation function calculated from the electric potential &
         & as a function of ky and parallel sep", "1", phi_corr_2pi)
    !if (proc0) call nc_loop_corr (nout,phi_corr_2pi)
    deallocate (phi_corr_2pi)
  end subroutine write_correlation

  subroutine correlation (cfnc_2pi)
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use fields_arrays, only: phinew
    implicit none
    complex, dimension (-ntgrid:,:), intent (out) :: cfnc_2pi
    integer :: ik, it, ig
    real :: fac
    
    cfnc_2pi = 0.0
    
    do ik = 1, naky
       if (ik==1) then
          fac = 1.0
       else
          fac = 0.5
       end if
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             cfnc_2pi(ig,ik) = cfnc_2pi(ig,ik) + phinew(0,it,ik)*conjg(phinew(ig,it,ik))*fac
          end do
       end do
    end do
  end subroutine correlation
  subroutine write_correlation_extend(gnostics)
    use theta_grid, only: ntgrid
    use kt_grids, only: jtwist_out, ntheta0, naky
    !use mp, only: proc0
    use gs2_io, only: nc_loop_corr_extend
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: t,t_old
    integer :: istep
    complex, dimension (:,:,:), allocatable, save:: phicorr_sum
    real, dimension (:,:,:), allocatable, save :: phiextend_sum
    complex, dimension (:,:,:), allocatable :: phi_corr
    real, dimension (:,:,:), allocatable :: phi2_extend
    real, save :: tcorr0 = 0.0
    
    t= gnostics%user_time
    t_old = gnostics%user_time_old
    istep = gnostics%istep
    
    if (.not. allocated(phicorr_sum)) then
       ntg_extend = (2*ntgrid+1)*((ntheta0-1)/jtwist_out+1)
       nth0_extend = min(ntheta0,jtwist_out*(naky-1))
       allocate (phicorr_sum(ntg_extend,ntheta0,naky)) ; phicorr_sum = 0.0
       allocate (phiextend_sum(ntg_extend,ntheta0,naky)) ; phiextend_sum = 0.0
       tcorr0 = t_old
    end if
    
    allocate (phi_corr(ntg_extend,ntheta0,naky))
    allocate (phi2_extend(ntg_extend,ntheta0,naky))
    call correlation_extend (phi_corr, phi2_extend)
    phicorr_sum = phicorr_sum + phi_corr*(t-t_old)
    phiextend_sum = phiextend_sum + phi2_extend*(t-t_old)
    call create_and_write_variable(gnostics, gnostics%rtype, "phi_corr", &
         dim_string([gnostics%dims%ri,gnostics%dims%theta_ext,gnostics%dims%kx,&
         gnostics%dims%ky,gnostics%dims%time]), &
         "2 point correlation function calculated from the electric potential &
         & as a function of ky and parallel sep... time averaged and &
         & calculated along the extended domain. ", "1", phicorr_sum/(t-tcorr0))
    !if (proc0 .and. mod(istep,nwrite_mult*nwrite)==0) then
       !call nc_loop_corr_extend (nout_big, phicorr_sum/(t-tcorr0), phiextend_sum/(t-tcorr0))
       !nout_big = nout_big + 1
    !end if
    deallocate (phi_corr, phi2_extend)
  end subroutine write_correlation_extend

  subroutine correlation_extend (cfnc,phi2extend)
!    use constants, only: pi
    use fields_arrays, only: phinew
    use theta_grid, only: ntgrid, jacob, delthet
!   use theta_grid, only: theta
    use kt_grids, only: ntheta0, naky, jtwist_out
    implicit none
    real, dimension (:,:,:), intent (out) :: phi2extend
    complex, dimension (:,:,:), intent (out) :: cfnc
    integer :: ig, it, ik, im, igmod
    integer :: itshift, nconnect, offset
!    real :: fac
    real, dimension (:), allocatable :: dl_over_b
    complex, dimension (:,:,:), allocatable :: phiextend
    complex, dimension (:,:), allocatable :: phir
!    real, dimension (:), allocatable :: phisum

    allocate (dl_over_b(2*ntgrid+1))
    dl_over_b = delthet*jacob

!    allocate (theta_extend(ntg_extend)) ; theta_extend = 0.0
!    do ig = 1, (ntheta0-1)/jtwist_out+1
!       theta_extend((ig-1)*(2*ntgrid+1)+1:ig*(2*ntgrid+1)) = pi+theta+(ig-1)*2.*pi
!    end do
!    theta_extend = theta_extend - theta_extend(size(theta_extend))/2.

    allocate (phiextend(ntg_extend,ntheta0,naky)) ; phiextend = 0.0
    allocate (phir(-ntgrid:ntgrid,ntheta0))
!    allocate (phisum(naky))

!    phisum = 0.0 ; cfnc_2pi = 0.0
!     do ik = 1, naky
!        if (ik==1) then
!           fac = 1.0
!        else
!           fac = 0.5
!        end if
!        ! integrate over x
!        do it = 1, ntheta0
!           do ig = -ntgrid, ntgrid
!              cfnc_2pi(ig,ik) = cfnc_2pi(ig,ik) + phinew(0,it,ik)*conjg(phinew(ig,it,ik))*fac
!           end do
! !          phisum(ik) = phisum(ik) + phinew(0,it,ik)*conjg(0,it,ik)*fac
!        end do
! !       cfnc_2pi(:,ik) = cfnc_2pi(:,ik)/phisum(ik)
!     end do

    cfnc = 0.0 ; phiextend = 0.0
    offset = ((ntheta0-1)/jtwist_out)*(2*ntgrid+1)/2
    do ig = -ntgrid, ntgrid
       call reorder_kx (phinew(ig,:,1), phir(ig,:))
    end do
    phiextend(offset+1:offset+(2*ntgrid+1),:,1) = phir
    do it = 1, ntheta0
       do im = 1, 2*ntgrid+1
          do ig = im+offset, offset+(2*ntgrid+1)
             igmod = mod(ig-offset-1,2*ntgrid+1)+1
             cfnc(im,it,1) = cfnc(im,it,1) &
                  + phiextend(ig,it,1)*conjg(phiextend(ig-im+1,it,1)) &
                  * dl_over_b(igmod)
          end do
       end do
       cfnc(:,it,1) = cfnc(:,it,1) / cfnc(1,it,1)
    end do
    
    do ik = 2, naky
       do ig = -ntgrid, ntgrid
          call reorder_kx (phinew(ig,:,ik), phir(ig,:))
       end do
       ! shift in kx due to parallel boundary condition
       ! also the number of independent theta0s
       itshift = jtwist_out*(ik-1)
       do it = 1, min(itshift,ntheta0)
          ! number of connections between kx's
          nconnect = (ntheta0-it)/itshift
          ! shift of theta index to account for fact that not all ky's
          ! have same length in extended theta
          offset = (2*ntgrid+1)*((ntheta0-1)/jtwist_out-nconnect)/2
          do ig = 1, nconnect+1
             phiextend(offset+(ig-1)*(2*ntgrid+1)+1:offset+ig*(2*ntgrid+1),it,ik) &
                  = phir(:,ntheta0-it-(ig-1)*itshift+1)
          end do
          do im = 1, (2*ntgrid+1)*(nconnect+1)
             do ig = im+offset, offset+(2*ntgrid+1)*(nconnect+1)
                igmod = mod(ig-offset-1,2*ntgrid+1)+1
                cfnc(im,it,ik) = cfnc(im,it,ik) &
                     + phiextend(ig,it,ik)*conjg(phiextend(ig-im+1,it,ik)) &
                     * dl_over_b(igmod)
             end do
          end do
          cfnc(:,it,ik) = cfnc(:,it,ik) / cfnc(1,it,ik)
       end do
    end do
    
    phi2extend = real(phiextend*conjg(phiextend))

!    deallocate (dl_over_b, phir, phiextend, phisum)
    deallocate (dl_over_b, phir, phiextend)
    
  end subroutine correlation_extend

  subroutine reorder_kx (unord, ord)
    use kt_grids, only: ntheta0
    implicit none
    complex, dimension (:), intent (in) :: unord
    complex, dimension (:), intent (out) :: ord
    
    ord(:ntheta0/2) = unord(ntheta0/2+2:)
    ord(ntheta0/2+1:) = unord(:ntheta0/2+1)
    
  end subroutine reorder_kx

  subroutine write_cross_phase(gnostics)
    use mp, only: proc0
    use diagnostics_config, only: diagnostics_type
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real :: phase_tot, phase_theta
    real :: t
    
    t = gnostics%user_time
    
    call get_cross_phase (gnostics, phase_tot, phase_theta)
    
    call create_and_write_variable(gnostics, gnostics%rtype, "electron_cross_phase_theta", &
         dim_string([gnostics%dims%time]), &
         "Cross phase between electron density and perpendicular electron temperature, &
         &  given at the outboard midplane and averaged across x and y", "radians", phase_theta)
    call create_and_write_variable(gnostics, gnostics%rtype, "electron_cross_phase_tot", &
         dim_string([gnostics%dims%time]), &
         "Cross phase between electron density and perpendicular electron temperature, &
         & averaged across all space", "radians", phase_tot)
    
    if (proc0 .and. gnostics%write_ascii) write (unit=gnostics%ascii_files%phase, &
         fmt="('t= ',e17.10,' phase_tot= ',e11.4,' phase_theta= ',e11.4)") &
         & t, phase_tot, phase_theta
  end subroutine write_cross_phase

  !> This is a highly simplified synthetic diagnostic which 
  !! calculates the cross phase between the electron density and the 
  !! perpendicular electron temperature for comparisons with DIII-D.  
  !! Returns the value of the cross-phase at the outboard midplane and 
  !! integrated over all v and x. We can generalize this routine to 
  !! other fields at some point, but for now this is just a skeleton for 
  !! a more realistic synthetic diagnostic. 
  subroutine get_cross_phase (gnostics, phase_tot, phase_theta)
    use species, only: nspec, spec, electron_species
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    use dist_fn, only: getemoms
    use fields_arrays, only: phinew, bparnew
    use volume_averages, only: average_all
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(in) :: gnostics
    real, intent (out) :: phase_tot, phase_theta
    complex, dimension (:,:,:,:), allocatable :: ntot, tperp
    !complex, dimension (ntheta0, naky) :: nTp_by_mode
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
          !call get_vol_int (ntot(0,:,:,is), tperp(0,:,:,is), nTp, nTp_by_mode)
          call average_all(ntot(0,:,:,is), tperp(0,:,:,is), nTp, gnostics%distributed)
!!!          call get_vol_average (ntot(0,:,:,is), ntot(0,:,:,is), n2, n2_by_mode)
          !          call get_vol_average (tperp(0,:,:,is), tperp(0,:,:,is), T2, T2_by_mode)
          phase_theta = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
          ! get integrated cross_phase 
          call average_all(ntot(:,:,:,is), tperp(:,:,:,is), nTp, gnostics%distributed)
          !call get_vol_int (ntot(:,:,:,is), tperp(:,:,:,is), nTp, nTp_by_mode)
          !          call get_vol_average (ntot(:,:,:,is), ntot(:,:,:,is), n2, n2_by_mode)
          !          call get_vol_average (tperp(:,:,:,is), tperp(:,:,:,is), T2, T2_by_mode)
          phase_tot = atan2(aimag(nTp),real(nTp))!/sqrt(n2*T2)
       end if
    end do
    
    deallocate (ntot, tperp)
    
  end subroutine get_cross_phase
end module diagnostics_turbulence
