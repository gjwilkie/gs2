module vpamu_grids

  implicit none

  public :: init_vpamu_grids, finish_vpamu_grids
  public :: integrate_moment, integrate_species
  public :: integrate_mu
  public :: vpa, nvgrid, wgts_vpa, dvpa, vpac
  public :: mu, nmu, wgts_mu, vpa_imp
  public :: vperp2, energy, anon, anonc, anon_thetac
  public :: use_gaussian_quadrature
  
  integer :: nvgrid
  integer :: nmu
  real :: vpa_imp, mu_max, vpa_max
  logical :: use_gaussian_quadrature
  
  ! arrays that are filled in vpamu_grids
  real, dimension (:), allocatable :: mu, wgts_mu
  real, dimension (:), allocatable :: vpa, wgts_vpa, dvpa
  real, dimension (:,:), allocatable :: vpac

  ! vpa-mu related arrays that are declared here
  ! but allocated and filled elsewhere because they depend on theta, etc.
  real, dimension (:,:), allocatable :: vperp2
  real, dimension (:,:,:), allocatable :: energy, anon, anonc
  real, dimension (:,:), allocatable :: anon_thetac

contains

  subroutine init_vpamu_grids

    implicit none

    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_vpa_grid
    call init_mu_grid

  end subroutine init_vpamu_grids

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    namelist /vpamu_grids_parameters/ nvgrid, nmu, vpa_max, mu_max, vpa_imp, &
         use_gaussian_quadrature

    integer :: in_file
    logical :: exist

    if (proc0) then

       nvgrid = 8
       vpa_max = 3.0
       nmu = 16
       mu_max = 4.5
       vpa_imp = 0.6
       use_gaussian_quadrature = .true.

       in_file = input_unit_exist("vpamu_grids_parameters", exist)
       if (exist) read (unit=in_file, nml=vpamu_grids_parameters)

    end if

    call broadcast (nvgrid)
    call broadcast (vpa_max)
    call broadcast (nmu)
    call broadcast (mu_max)
    call broadcast (vpa_imp)
    call broadcast (use_gaussian_quadrature)

  end subroutine read_parameters

  subroutine init_vpa_grid

    use centering, only: get_cell_value

    implicit none

    integer :: iv, i, idx, iseg, nvpa_seg
    real :: del

    if (.not. allocated(vpa)) then
       ! vpa is the parallel velocity at grid points
       allocate (vpa(-nvgrid:nvgrid)) ; vpa = 0.0
       ! wgts_vpa are the integration weights assigned
       ! to the parallel velocity grid points
       allocate (wgts_vpa(-nvgrid:nvgrid)) ; wgts_vpa = 0.0
       ! vpac is the parallel velocity at cells
       ! the 2nd dimension keeps track of dvpa/dt < 0 (1) and dvpa/dt > 0 (2)
       ! which for lowest order GK equation corresponds to theta < 0 and theta > 0
       allocate (vpac(-nvgrid:nvgrid,2)) ; vpac = 0.0
       ! dvpa is the grid spacing in vpa
       allocate (dvpa(-nvgrid:nvgrid)) ; dvpa = 0.0
    end if

    ! velocity grid goes from -vpa_max to vpa_max
    ! with a point at vpa = 0

    ! obtain vpa grid for vpa >= 0
    do iv = 0, nvgrid
       vpa(iv) = real(iv)*vpa_max/nvgrid
    end do
    ! fill in vpa grid for vpa < 0
    vpa(-nvgrid:-1) = -vpa(nvgrid:1:-1)

    dvpa(-nvgrid:nvgrid-1) = (/ (vpa(i+1)-vpa(i), i=-nvgrid,nvgrid-1) /)
    ! dvpa(nvgrid) should never be needed, but give it a nonzero
    ! value so division by zero never occurs when doing matrix operations
    dvpa(nvgrid) = dvpa(nvgrid-1)

    ! get integration weights corresponding to vpa grid points
    ! for now use Simpson's rule; 
    ! i.e. subdivide grid into 3-point segments, with each segment spanning vpa_low to vpa_up
    ! then the contribution of each segment to the integral is
    ! (vpa_up - vpa_low) * (f1 + 4*f2 + f3) / 6
    ! inner boundary points are used in two segments, so they get double the weight
    nvpa_seg = (2*nvgrid+1)/2
    do iseg = 1, nvpa_seg
       idx = -nvgrid + (iseg-1)*2
       del = (dvpa(idx)+dvpa(idx+1))/6.
       wgts_vpa(idx) = wgts_vpa(idx) + del
       wgts_vpa(idx+1) = wgts_vpa(idx+1) + 4.*del
       wgts_vpa(idx+2) = wgts_vpa(idx+2) + del
    end do

    ! get vpa at cell centers
    ! example: if vpa_imp=1 (fully upwinded), vpac(iv,1)=vpa(iv+1)
    ! and vpac(iv,2)=vpa(iv)
    call get_cell_value (vpa_imp, vpa, vpac(:,1), -nvgrid)
    call get_cell_value (1.0-vpa_imp, vpa, vpac(:,2), -nvgrid)
    ! vpac(nvgrid) should not be needed but set it nonzero to avoid possible
    ! divide by zero
    vpac(nvgrid,:) = vpa(nvgrid)

  end subroutine init_vpa_grid

  subroutine integrate_mu (ig, g, total)

    use species, only: nspec
    use theta_grid, only: bmag

    implicit none

    integer, intent (in) :: ig
    real, dimension (:,:), intent (in) :: g
    real, dimension (:), intent (out) :: total

    integer :: is, imu

    total = 0.

    do is = 1, nspec
       ! sum over mu
       do imu = 1, nmu
          total(is) = total(is) + wgts_mu(imu)*bmag(ig)*g(imu,is)
       end do
    end do

  end subroutine integrate_mu

  ! integrate over v-space for each species
  subroutine integrate_moment (g, total, all)

    use mp, only: nproc, sum_reduce, sum_allreduce
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx
    use theta_grid, only: ntgrid, bmag

    implicit none

    integer :: iglo, iv, ik, is, imu, ig

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    integer, intent (in), optional :: all

    total = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)

       do ig = -ntgrid, ntgrid
          ! sum over parallel speeds
          do iv = -nvgrid, nvgrid
             total(ig,:,ik,is) = total(ig,:,ik,is) + &
                  wgts_mu(imu)*wgts_vpa(iv)*bmag(ig)*g(ig,iv,:,iglo)
          end do
       end do
    end do

    if (nproc > 1) then
       if (present(all)) then
          call sum_allreduce (total)
       else
          call sum_reduce (total, 0)
       end if
    end if

  end subroutine integrate_moment

  ! integrave over v-space and sum over species
  subroutine integrate_species (g, weights, total)

    use mp, only: sum_allreduce
    use gs2_layouts, only: g_lo, ik_idx, imu_idx, is_idx
    use theta_grid, only: ntgrid, bmag

    implicit none

    integer :: iglo, iv, ik, ig, is, imu

    complex, dimension (-ntgrid:,-nvgrid:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (-ntgrid:,:,:), intent (out) :: total

    total = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       imu = imu_idx(g_lo,iglo)

       do ig = -ntgrid, ntgrid
          ! sum over parallel speeds
          do iv = -nvgrid, nvgrid
             total(ig,:,ik) = total(ig,:,ik) + &
                  wgts_mu(imu)*wgts_vpa(iv)*bmag(ig)*g(ig,iv,:,iglo)*weights(is)
          end do
       end do
    end do

    call sum_allreduce (total)

  end subroutine integrate_species

  subroutine finish_vpa_grid

    implicit none

    if (allocated(vpa)) deallocate (vpa)
    if (allocated(wgts_vpa)) deallocate (wgts_vpa)
    if (allocated(dvpa)) deallocate (dvpa)
    if (allocated(vpac)) deallocate (vpac)

  end subroutine finish_vpa_grid

  subroutine init_mu_grid

    use constants, only: pi
    use gauss_quad, only: get_laguerre_grids
    
    implicit none

    integer :: imu, iseg, nmu_seg, idx, i
    real :: del
    real, dimension (:), allocatable :: dmu

    ! allocate arrays and initialize to zero
    if (.not. allocated(mu)) then
       allocate (mu(nmu)) ; mu = 0.0
       allocate (wgts_mu(nmu)) ; wgts_mu = 0.0
    end if

    allocate (dmu(nmu-1)) ; dmu = 0.0
    
    if (use_gaussian_quadrature) then
       ! use Gauss-Laguerre quadrature in mu
       call get_laguerre_grids (mu, wgts_mu)
       wgts_mu = wgts_mu*exp(mu)
    else

       ! construct mu grid
       ! first get sqrt(mu) grid to be equally spaced
       do imu = 1, nmu
          mu(imu) = real(imu-1)*sqrt(mu_max)/(nmu-1)
       end do
       dmu = (/ (mu(i+1)-mu(i), i=1,nmu-1) /)

       ! get integration weights corresponding to mu grid points
       ! for now use Simpson's rule; 
       ! i.e. subdivide grid into 3-point segments, with each segment spanning mu_low to mu_up
       ! then the contribution of each segment to the integral is
       ! (mu_up - mu_low) * (f1 + 4*f2 + f3) / 6
       ! inner boundary points are used in two segments, so they get double the weight

       ! first get number of segments with 3 points.  If nmu is even, there will be 
       ! an extra segment with only 2 points that we will treat separately
       nmu_seg = (nmu-1)/2
       do iseg = 1, nmu_seg
          idx = (iseg-1)*2 + 1
          del = (dmu(idx)+dmu(idx+1))/6.
          wgts_mu(idx) = wgts_mu(idx) + del*mu(idx)
          wgts_mu(idx+1) = wgts_mu(idx+1) + 4.*del*mu(idx+1)
          wgts_mu(idx+2) = wgts_mu(idx+2) + del*mu(idx+2)
       end do
       ! if there is an extra segment with only 2 points in it,
       ! assign weights using trapezoid rule
       if (mod(nmu,2)==0) wgts_mu(nmu-1:) = wgts_mu(nmu-1:) + dmu(nmu-1)*0.5*mu(nmu-1)

       ! now convert from sqrt(mu) to mu
       mu = mu**2
       dmu = (/ (mu(i+1)-mu(i), i=1,nmu-1) /)
    end if
       
    ! factor of 2./sqrt(pi) necessary to account for 2pi from 
    ! integration over gyro-angle and 1/pi^(3/2) normalization
    ! of velocity space Jacobian
    ! note that a factor of bmag is missing and will have to be
    ! applied when doing integrals
    wgts_mu = wgts_mu*4./sqrt(pi)

    deallocate (dmu)

  end subroutine init_mu_grid

  subroutine finish_mu_grid

    implicit none

    if (allocated(mu)) deallocate (mu)
    if (allocated(wgts_mu)) deallocate (wgts_mu)

  end subroutine finish_mu_grid

  subroutine finish_vpamu_grids

    implicit none
    
    call finish_vpa_grid
    call finish_mu_grid

  end subroutine finish_vpamu_grids

end module vpamu_grids
