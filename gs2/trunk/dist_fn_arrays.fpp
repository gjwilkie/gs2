!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays. 

module dist_fn_arrays

  public :: g, gnew, gold, kx_shift, theta0_shift, vpa, vpac
  public :: vperp2, vpar, ittp, aj0, aj1, aj2, aj0f, aj1f
  public :: apar_ext, kperp2, c_rate
  public :: g_adjust
#ifdef LOWFLOW
  public :: hneoc, vparterm, wdfac, wstarfac
#endif


  ! dist fn
  complex, dimension (:,:,:), allocatable :: g, gnew, gold
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension(:), allocatable :: kx_shift, theta0_shift
  ! (naky)

  real, dimension (:,:,:), allocatable :: vpa, vpac
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:), allocatable :: vperp2
  ! (-ntgrid:ntgrid, -g-layout-)

  real, dimension (:,:,:), allocatable :: vpar
  ! (-ntgrid:ntgrid,2, -g-layout-)

  integer, dimension (:), allocatable :: ittp
  ! (-ntgrid:ntgrid)

! DJA: 17/1/06, add variable aj2 to store J_2(x)
  real, dimension (:,:), allocatable :: aj0, aj1, aj2, aj0f, aj1f
  ! (-ntgrid:ntgrid, -g-layout-)

  ! fieldeq
  complex, dimension (:,:,:), allocatable :: apar_ext
  real, dimension (:,:,:), allocatable :: kperp2
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  ! collisional diagnostic of heating rate
  complex, dimension (:,:,:,:,:), allocatable :: c_rate
  ! (-ntgrid:ntgrid,ntheta0,naky,nspecies,2) replicated

#ifdef LOWFLOW
  ! v-dependent factors in low-flow terms
  real, dimension (:,:,:), allocatable :: hneoc, vparterm, wdfac, wstarfac
  ! (-ntgrid:ntgrid,2, -g-layout-)
#endif

  private
contains

  subroutine g_adjust (g, phi, bpar, facphi, facbpar)
    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: anon
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar

    integer :: iglo, ig, ik, it, ie, is
    complex :: adj

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
! BD:  bpar == delta B_parallel / B_0(theta) so no extra factor of 
! 1/bmag is needed here.
       do ig = -ntgrid, ntgrid
           adj = anon(ie)*2.0*vperp2(ig,iglo)*aj1(ig,iglo) &
                  *bpar(ig,it,ik)*facbpar &
               + spec(is)%z*anon(ie)*phi(ig,it,ik)*aj0(ig,iglo) &
                  /spec(is)%temp*facphi
          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do
  end subroutine g_adjust

end module dist_fn_arrays
