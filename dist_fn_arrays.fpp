!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays. 

module dist_fn_arrays

  public :: g, gnew, gold, kx_shift, theta0_shift, source
  public :: gpnew, gpold!, ghnew, ghold
  public :: vpar, aj0, aj1, vparp, aj0p
  public :: apar_ext, kperp2, dkperp2dr
  public :: g_adjust
#ifdef LOWFLOW
  public :: hneoc, vparterm, wdfac, wstarfac
#endif


  ! dist fn
  complex, dimension (:,:,:,:), allocatable :: g, gnew, gold, source
  complex, dimension (:,:,:,:), allocatable :: gpnew, gpold!, ghnew, ghold
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, -g-layout-)

  real, dimension(:), allocatable :: kx_shift, theta0_shift
  ! (naky)

  real, dimension (:,:,:), allocatable :: vpar, vparp
  ! (-ntgrid:ntgrid, -nvgrid:nvgrid, nspec)

  real, dimension (:,:,:), allocatable :: aj0, aj1, aj0p
  ! (-ntgrid:ntgrid, ntheta0, -g-layout-)

  ! fieldeq
  complex, dimension (:,:,:), allocatable :: apar_ext
  real, dimension (:,:,:), allocatable :: kperp2, dkperp2dr
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

#ifdef LOWFLOW
  ! v-dependent factors in low-flow terms
  real, dimension (:,:,:), allocatable :: hneoc, vparterm, wdfac, wstarfac
  ! (-ntgrid:ntgrid,2, -g-layout-)
#endif

  private
contains

  subroutine g_adjust (g, phi, bpar, facphi, facbpar)
!CMR, 17/4/2012: 
! g_adjust transforms between representations of perturbed dist'n func'n.
!    <delta_f> = g_wesson J0(Z) - q phi/T F_m  where <> = gyroaverage
!        g_gs2 = g_wesson - q phi/T J0(Z) F_m - m v_||^2/T B_||/B J1(Z)/Z F_m
! For numerical convenience the GS2 variable g uses the form g_gs2.
! g_wesson (see Wesson's book, Tokamaks) is often a more convenient form:
!     e.g. for handling collisions, calculating v-space moments in real space.
!
! To transform gnew from g_gs2 to g_wesson form:
!    call g_adjust(gnew,phinew,bparnew,fphi,fbpar)
! or transform from gnew from g_wesson to g_gs2 form:
!    call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
!
    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid, vperp2, anon
    use gs2_layouts, only: g_lo, ik_idx, is_idx, imu_idx

    implicit none
    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar

    integer :: iglo, ig, ik, it, is, imu, iv, ntheta0
    complex :: adj

    ntheta0 = size(g,3)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
! BD:  bpar == delta B_parallel / B_0(theta) so no extra factor of 
! 1/bmag is needed here.
       do it = 1, ntheta0
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                adj = (2.0*vperp2(ig,imu)*aj1(ig,it,iglo) &
                     *bpar(ig,it,ik)*facbpar &
                     + spec(is)%z*phi(ig,it,ik)*aj0(ig,it,iglo) &
                     /spec(is)%temp*facphi) &
                     * anon(ig,iv,imu)
                g(ig,iv,it,iglo) = g(ig,iv,it,iglo) + adj
             end do
          end do
       end do
    end do
  end subroutine g_adjust

end module dist_fn_arrays
