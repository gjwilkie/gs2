module gs2_transforms

  use redistribute, only: redist_type

  implicit none

  public :: init_transforms, init_zf, kz_spectrum
  public :: init_x_transform
  public :: transform_x, transform_y, transform2
  public :: inverse_x, inverse_y, inverse2

  private

  interface transform_x
     module procedure transform_x5d
     module procedure transform_x3d
     module procedure transform_x1d
  end interface

  interface transform_y
     module procedure transform_y5d
     module procedure transform_y3d
     module procedure transform_y1d
  end interface

  interface transform2
     module procedure transform2_5d, transform2_5d_accel
     module procedure transform2_3d
  end interface

  interface inverse_x
     module procedure inverse_x5d
     module procedure inverse_x3d
     module procedure inverse_x1d
  end interface

  interface inverse_y
     module procedure inverse_y5d
     module procedure inverse_y3d
     module procedure inverse_y1d
  end interface

  interface inverse2
     module procedure inverse2_5d, inverse2_5d_accel
     module procedure inverse2_3d
  end interface

  ! redistribution

  type (redist_type) :: g2x, x2y

  ! fft

  complex, dimension (:,:), allocatable :: xxf_tmp, fft_tmp

  real, dimension (:), allocatable :: fft_wrk, xtable, ytable
  real :: xscale, yscale

! accel will be set to true if the v layout is used AND the number of
! PEs is such that each PE has a complete copy of the x,y space --
! in that case, no communication is needed to evaluate the nonlinear
! terms
  logical :: accel = .false.

  integer :: igmin_proc, igmax_proc
  integer, dimension (:), allocatable :: igproc, ia, iak
  logical, dimension (:), allocatable :: aidx  ! aidx == aliased index

! k_parallel filter items
  real, dimension (:), allocatable :: kpar  ! (2*ntgrid)
  real, dimension(:), allocatable :: zwork, ztable
  real :: zscale
  integer :: nztable, nzwork


contains

  subroutine init_transforms &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerated)
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    logical, intent (out) :: accelerated

    logical, save :: initialized = .false.
    character (1) :: char

    if (initialized) return
    initialized = .true.

    call init_3d_layouts

!    call pe_layout (char) 

!    if (char == 'v' .and. mod (negrid*nlambda*nspec, nproc) == 0) then
!       accel = .true.
!       call init_accel_transform_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
!    else
       call init_y_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
!    end if

    call init_y_fft (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
    accelerated = accel

  end subroutine init_transforms

  subroutine init_x_transform (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    integer :: nwork

    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_3d_layouts

    call init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    call init_x_fft (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nwork)
    if (allocated(fft_wrk)) deallocate (fft_wrk) 
    allocate (fft_wrk(nwork))
  end subroutine init_x_transform

  subroutine init_x_fft (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nwork)
    use gs2_layouts, only: xxf_lo
    use fft_work, only: ccfft_work0, ccfft_work1, ccfft_table0, ccfft_table1
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (out) :: nwork
    integer :: nx
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    nx = xxf_lo%nx
    nwork = ccfft_work0 + ccfft_work1*nx
    xscale = 1.0/real(nx)
    allocate (xtable(ccfft_table0 + ccfft_table1*nx))
    call ccfft (0, nx, 1.0, xtable, xtable, xtable, xtable, 0)
  end subroutine init_x_fft

  subroutine init_y_fft (ntgrid, naky, ntheta0, nlambda, negrid, nspec)
    use gs2_layouts, only: yxf_lo
    use fft_work, only: csfft_work0, csfft_work1, csfft_table0, csfft_table1
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer :: nwork
    integer :: ny
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_x_fft (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nwork)

    ny = yxf_lo%ny
    nwork = max(nwork,csfft_work0 + csfft_work1*ny)
    yscale = 1.0/real(ny)
    allocate (ytable(csfft_table0 + csfft_table1*ny))
    call scfft (0, ny, 1.0, ytable, ytable, ytable, ytable, 0)
    if (allocated(fft_wrk)) deallocate (fft_wrk) 
    allocate (fft_wrk(nwork))
  end subroutine init_y_fft

  subroutine init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    use gs2_layouts, only: init_x_transform_layouts
    use gs2_layouts, only: g_lo, xxf_lo, gidx2xxfidx, proc_id, idx_local
    use mp, only: iproc, nproc
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (2) :: to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx
    logical :: initialized = .false.

    integer :: iglo, isign, ig, it, ixxf
    integer :: n, ip

    if (initialized) return
    initialized = .true.

    call init_x_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
             if (idx_local(g_lo,iglo)) &
                  nn_from(proc_id(xxf_lo,ixxf)) = nn_from(proc_id(xxf_lo,ixxf)) + 1
             if (idx_local(xxf_lo,ixxf)) &
                nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
          end do
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do isign = 1, 2
          do ig = -ntgrid, ntgrid
             call gidx2xxfidx (ig, isign, iglo, g_lo, xxf_lo, it, ixxf)
             if (idx_local(g_lo,iglo)) then
                ip = proc_id(xxf_lo,ixxf)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = ig
                from_list(ip)%second(n) = isign
                from_list(ip)%third(n) = iglo
             end if
             if (idx_local(xxf_lo,ixxf)) then
                ip = proc_id(g_lo,iglo)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = it
                to_list(ip)%second(n) = ixxf
             end if
          end do
       end do
    end do

    from_low (1) = -ntgrid
    from_low (2) = 1
    from_low (3) = g_lo%llim_proc

    to_low = xxf_lo%llim_proc
    
    to_high(1) = xxf_lo%nx
    to_high(2) = xxf_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = 2
    from_high(3) = g_lo%ulim_alloc

    call init_redist (g2x, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_x_redist

  subroutine init_y_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    use gs2_layouts, only: init_y_transform_layouts
    use gs2_layouts, only: xxf_lo, yxf_lo, xxfidx2yxfidx, proc_id, idx_local
    use mp, only: iproc, nproc
    use redistribute, only: index_list_type, init_redist, delete_list
    implicit none
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (2) :: from_low, from_high, to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny

    integer :: it, ixxf, ik, iyxf
    integer :: n, ip
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)
    call init_y_transform_layouts &
         (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do ixxf = xxf_lo%llim_world, xxf_lo%ulim_world
       do it = 1, yxf_lo%nx
          call xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
          if (idx_local(xxf_lo,ixxf)) &
             nn_from(proc_id(yxf_lo,iyxf)) = nn_from(proc_id(yxf_lo,iyxf)) + 1
          if (idx_local(yxf_lo,iyxf)) &
             nn_to(proc_id(xxf_lo,ixxf)) = nn_to(proc_id(xxf_lo,ixxf)) + 1
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0
    do ixxf = xxf_lo%llim_world, xxf_lo%ulim_world
       do it = 1, yxf_lo%nx
          call xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
          if (idx_local(xxf_lo,ixxf)) then
             ip = proc_id(yxf_lo,iyxf)
             n = nn_from(ip) + 1
             nn_from(ip) = n
             from_list(ip)%first(n) = it
             from_list(ip)%second(n) = ixxf
          end if
          if (idx_local(yxf_lo,iyxf)) then
             ip = proc_id(xxf_lo,ixxf)
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_list(ip)%first(n) = ik
             to_list(ip)%second(n) = iyxf
          end if
       end do
    end do

    from_low(1) = 1
    from_low(2) = xxf_lo%llim_proc

    to_low = yxf_lo%llim_proc

    to_high(1) = yxf_lo%ny/2+1
    to_high(2) = yxf_lo%ulim_alloc

    from_high(1) = xxf_lo%nx
    from_high(2) = xxf_lo%ulim_alloc

    call init_redist (x2y, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)
    
  end subroutine init_y_redist

  subroutine transform_x5d (g, xxf)
    use gs2_layouts, only: xxf_lo, g_lo
    use mp, only: iproc, nproc, send, receive
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather
    implicit none
    complex, dimension (-xxf_lo%ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    integer :: i, idp, ipto, ipfrom, iadp

    call prof_entering ("transform_x5d", "gs2_transforms")

    ! intent statement in gather actually makes this next line non-standard: 
    xxf = 0.
    call gather (g2x, g, xxf)

    ! do ffts (no memory problem here)
    do i = xxf_lo%llim_proc, xxf_lo%ulim_proc
       call ccfft (1, xxf_lo%nx, 1.0, xxf(:,i), xxf(:,i), xtable, fft_wrk, 0)
    end do

    call prof_leaving ("transform_x5d", "gs2_transforms")
  end subroutine transform_x5d

  subroutine inverse_x5d (xxf, g)
    use gs2_layouts, only: xxf_lo, g_lo
    use mp, only: iproc, nproc, send, receive
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter
    implicit none
    complex, dimension (:,xxf_lo%llim_proc:), intent (in out) :: xxf
    complex, dimension (-xxf_lo%ntgrid:,:,g_lo%llim_proc:), intent (out) :: g
    integer :: i, idp, ipto, ipfrom, iadp

    call prof_entering ("inverse_x5d", "gs2_transforms")

    ! do ffts
    do i = xxf_lo%llim_proc, xxf_lo%ulim_proc
       call ccfft (-1, xxf_lo%nx, xscale, xxf(:,i), xxf(:,i), xtable, fft_wrk, 0)
    end do

! not necessary.  double-checked on 9.1.99
!    call reality(xxf)

    g = 0.
    call scatter (g2x, xxf, g)

    call prof_leaving ("inverse_x5d", "gs2_transforms")
  end subroutine inverse_x5d

  subroutine transform_y5d (xxf, yxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather
    implicit none
    complex, dimension (:,xxf_lo%llim_proc:), intent (in) :: xxf
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf
    integer :: i

    call prof_entering ("transform_y5d", "gs2_transforms")

    fft_tmp = 0.
    call gather (x2y, xxf, fft_tmp)

    ! normalize
    fft_tmp(2:yxf_lo%naky,:) = fft_tmp(2:yxf_lo%naky,:)/2.0

    ! do ffts (no memory problem)
    do i = yxf_lo%llim_proc, yxf_lo%ulim_proc
       call csfft (1, yxf_lo%ny, 1.0, fft_tmp(:,i), yxf(:,i), ytable, fft_wrk, 0)
    end do

    call prof_leaving ("transform_y5d", "gs2_transforms")
  end subroutine transform_y5d

  subroutine inverse_y5d (yxf, xxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in) :: yxf
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    integer :: i
!    complex, dimension (yxf_lo%ny/2+1,yxf_lo%llim_proc:yxf_lo%ulim_alloc) :: fft

    call prof_entering ("inverse_y5d", "gs2_transforms")

    ! do ffts
    do i = yxf_lo%llim_proc, yxf_lo%ulim_proc
       call scfft (-1, yxf_lo%ny, yscale, yxf(:,i), fft_tmp(:,i), ytable, fft_wrk, 0)
       fft_tmp(2:,i) = fft_tmp(2:,i)*2.0
    end do

    call scatter (x2y, fft_tmp, xxf)

    call prof_leaving ("inverse_y5d", "gs2_transforms")
  end subroutine inverse_y5d

  subroutine transform2_5d (g, yxf)
    use gs2_layouts, only: g_lo, xxf_lo, yxf_lo
    implicit none
    complex, dimension (:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf

    call init_xf
    call transform_x (g, xxf_tmp)
    call transform_y (xxf_tmp, yxf)

  end subroutine transform2_5d

  subroutine inverse2_5d (yxf, g)
    use gs2_layouts, only: g_lo, xxf_lo, yxf_lo
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in) :: yxf
    complex, dimension (:,:,g_lo%llim_proc:), intent (out) :: g

    xxf_tmp = 0.; g = 0.           ! dealias
    call inverse_y (yxf, xxf_tmp)
    call inverse_x (xxf_tmp, g)

  end subroutine inverse2_5d

  subroutine transform2_5d_accel (g, yxf, i)

!!!!!!!!!!!!!!!!! Not operational  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use gs2_layouts, only: g_lo, xxf_lo, yxf_lo
    implicit none
    complex, dimension (:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:,:,yxf_lo%llim_proc:), intent (out) :: yxf
    integer :: i

  end subroutine transform2_5d_accel

  subroutine inverse2_5d_accel (yxf, g, i)

!!!!!!!!!!!!!!!!! Not operational  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use gs2_layouts, only: g_lo, xxf_lo, yxf_lo
    implicit none
    real, dimension (:,:,yxf_lo%llim_proc:), intent (in) :: yxf
    complex, dimension (:,:,g_lo%llim_proc:), intent (out) :: g
    integer :: i

  end subroutine inverse2_5d_accel

  subroutine init_3d_layouts
    use theta_grid, only: ntgrid
    use mp, only: nproc, iproc
    implicit none
    integer :: ntot, nblock, ig, i_proc, n_proc
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    allocate (igproc(-ntgrid:ntgrid))
    ntot = 2*ntgrid+1
    nblock = (ntot+1)/nproc+1

    igmin_proc = ntgrid+1
    igmax_proc = -ntgrid-1
    i_proc = 0
    n_proc = 0
    do ig = -ntgrid, ntgrid
       igproc(ig) = i_proc
       if (i_proc == iproc) then
          if (igmin_proc > ig) igmin_proc = ig
          if (igmax_proc < ig) igmax_proc = ig
       end if
       n_proc = n_proc + 1
       if (n_proc >= nblock) then
          n_proc = 0
          i_proc = i_proc + 1
       end if
    end do
  end subroutine init_3d_layouts

  subroutine transform_x3d (phi, phixxf)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: xxf_lo, yxf_lo
    use mp, only: broadcast
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    complex, dimension (:,:,-ntgrid:), intent (out) :: phixxf
    integer :: ig, ik, nlo

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    if (igmin_proc <= igmax_proc) then
       do ig = igmin_proc, igmax_proc
          do ik = 1, naky
             phixxf(:,ik,ig) = 0.0
             phixxf(1:(ntheta0+1)/2,ik,ig) = phi(ig,1:(ntheta0+1)/2,ik)
             phixxf(nlo:xxf_lo%nx,ik,ig) &
                  = phi(ig,(ntheta0+1)/2+1:ntheta0,ik)
             call ccfft (1, xxf_lo%nx, 1.0, phixxf(:,ik,ig), phixxf(:,ik,ig), &
                  xtable, fft_wrk, 0)
          end do
       end do
    end if
    do ig = -ntgrid, ntgrid
       do ik = 1, naky
          call broadcast (phixxf(:,ik,ig), igproc(ig))
       end do
    end do
  end subroutine transform_x3d

  subroutine inverse_x3d (phixxf, phi)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: xxf_lo, yxf_lo
    implicit none
    complex, dimension (:,:,-ntgrid:), intent (in) :: phixxf
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi
    complex, dimension (xxf_lo%nx) :: f
    integer :: ig, ik, nlo


! new field array layout fixes (never observed?) bug here!

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    do ig = -ntgrid, ntgrid
       do ik = 1, naky
          call ccfft (-1, xxf_lo%nx, xscale, phixxf(:,ik,ig), f, &
               xtable, fft_wrk, 0)
          phi(ig,1:(ntheta0+1)/2,ik) = f(1:(ntheta0+1)/2)
          phi(ig,(ntheta0+1)/2+1:ntheta0,ik) = f(nlo:xxf_lo%nx)
       end do
    end do
  end subroutine inverse_x3d

  subroutine transform_y3d (phixxf, phixf)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: xxf_lo, yxf_lo
    use mp, only: broadcast
    implicit none
    complex, dimension (:,:,-ntgrid:), intent (in) :: phixxf
    real, dimension (:,:,-ntgrid:), intent (out) :: phixf
    complex, dimension (yxf_lo%ny+1) :: f
    integer :: ig, ix

    if (igmin_proc <= igmax_proc) then
       do ig = igmin_proc, igmax_proc
          do ix = 1, xxf_lo%nx
             f(1) = real(phixxf(ix,1,ig))
             f(2:naky) = phixxf(ix,2:naky,ig)/2.0
             f(naky+1:) = 0.0
             call csfft (1, yxf_lo%nx, 1.0, f, phixf(:,ix,ig), &
                  ytable, fft_wrk, 0)
          end do
       end do
    end if
    do ig = -ntgrid, ntgrid
       do ix = 1, xxf_lo%nx
          call broadcast (phixf(:,ix,ig), igproc(ig))
       end do
    end do
  end subroutine transform_y3d

  subroutine inverse_y3d (phixf, phixxf)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: xxf_lo, yxf_lo
    implicit none
    real, dimension (:,:,-ntgrid:), intent (in) :: phixf
    complex, dimension (:,:,-ntgrid:), intent (out) :: phixxf
    complex, dimension (yxf_lo%ny+1) :: f
    integer :: ig, ix

    do ig = -ntgrid, ntgrid
       do ix = 1, xxf_lo%nx
          call scfft (-1, yxf_lo%ny, yscale, phixf(:,ix,ig), f, &
               ytable, fft_wrk, 0)
          phixxf(ix,1,ig) = f(1)
          phixxf(ix,2:naky,ig) = f(2:naky)*2.0
       end do
    end do
  end subroutine inverse_y3d

  subroutine transform2_3d (phi, phixf)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky
    use gs2_layouts, only: xxf_lo, yxf_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,-ntgrid:), intent (out) :: phixf
    complex, dimension (xxf_lo%nx,naky,-ntgrid:ntgrid) :: phixxf

    call transform_x (phi, phixxf)
    call transform_y (phixxf, phixf)
  end subroutine transform2_3d

  subroutine inverse2_3d (phi, phixf)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky
    use gs2_layouts, only: xxf_lo, yxf_lo
    implicit none
    real, dimension (:,:,-ntgrid:), intent (in) :: phixf
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi
    complex, dimension (xxf_lo%nx,naky,-ntgrid:ntgrid) :: phixxf

    call inverse_y (phixf, phixxf)
    call inverse_x (phixxf, phi)
  end subroutine inverse2_3d

  subroutine transform_x1d (f, xf)
    use kt_grids, only: ntheta0
    use gs2_layouts, only: xxf_lo
    implicit none
    complex, dimension (:), intent (in) :: f
    complex, dimension (:), intent (out) :: xf
    integer :: nlo

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    xf = 0.0
    xf(1:(ntheta0+1)/2) = f(1:(ntheta0+1)/2)
    xf(nlo:xxf_lo%nx) = f((ntheta0+1)/2+1:ntheta0)
    call ccfft (1, xxf_lo%nx, 1.0, xf, xf, xtable, fft_wrk, 0)
  end subroutine transform_x1d

  subroutine inverse_x1d (xf, f)
    use kt_grids, only: ntheta0
    use gs2_layouts, only: xxf_lo
    implicit none
    complex, dimension (:), intent (in out) :: xf
    complex, dimension (:), intent (out) :: f
    integer :: nlo

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd

    call ccfft (-1, xxf_lo%nx, xscale, xf, xf, xtable, fft_wrk, 0)
    f(1:(ntheta0+1)/2) = xf(1:(ntheta0+1)/2)
    f((ntheta0+1)/2+1:ntheta0) = xf(nlo:xxf_lo%nx)
  end subroutine inverse_x1d

  subroutine transform_y1d (f, xf)
    use kt_grids, only: naky
    use gs2_layouts, only: yxf_lo
    implicit none
    complex, dimension (:), intent (in) :: f
    real, dimension (:), intent (out) :: xf
    complex, dimension (yxf_lo%ny+1) :: fx

    fx(1) = real(f(1))
    fx(2:naky) = f(2:naky)/2.0
    fx(naky+1:) = 0.0
    call csfft (1, yxf_lo%ny, 1.0, fx, xf, ytable, fft_wrk, 0)
  end subroutine transform_y1d

  subroutine inverse_y1d (xf, f)
    use kt_grids, only: naky
    use gs2_layouts, only: yxf_lo
    implicit none
    real, dimension (:), intent (in) :: xf
    complex, dimension (:), intent (out) :: f
    complex, dimension (yxf_lo%ny+1) :: fx

    call scfft (-1, yxf_lo%ny, yscale, xf, fx, ytable, fft_wrk, 0)
    f(1) = fx(1)
    f(2:naky) = fx(2:naky)*2.0
  end subroutine inverse_y1d

  subroutine reality (xxf)

    use kt_grids, only: ntheta0
    use gs2_layouts, only: xxf_lo, xxf_ky_is_zero
    implicit none
    complex, dimension(:,xxf_lo%llim_proc:), intent (in out) :: xxf
    integer :: it, idx

    do idx = xxf_lo%llim_proc, xxf_lo%ulim_proc
       if (xxf_ky_is_zero(xxf_lo, idx)) then
          do it = 1, ntheta0/2
             xxf(xxf_lo%nx-ntheta0/2+it,idx) = conjg(xxf((ntheta0+1)/2+1-it,idx))
!             xxf((ntheta0+1)/2+1-it,idx) = conjg(xxf(xxf_lo%nx-ntheta0/2+it,idx))
          enddo
       endif
    enddo
   
  end subroutine reality

  subroutine init_xf

    use gs2_layouts, only: xxf_lo, yxf_lo
    logical :: initialized = .false.
        
    if (initialized) return
    initialized = .true.

    allocate (xxf_tmp(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc))
    allocate (fft_tmp(yxf_lo%ny/2+1,yxf_lo%llim_proc:yxf_lo%ulim_alloc))


  end subroutine init_xf

  subroutine init_zf (ntgrid, nperiod)

    use fft_work, only: ccfft_work0, ccfft_work1, ccfft_table0, ccfft_table1
    implicit none
    integer, intent (in) :: ntgrid, nperiod
    complex, dimension(2*ntgrid) :: b
    real :: L_theta
    integer :: i
    logical :: done = .false.

    if (done) return
    done = .true.

!    allocate (kpar(2*ntgrid))
    
!    do i = 1, ntgrid
!       kpar(i) = real(i-1)
!       kpar(i+ntgrid) = real(i-ntgrid-1)
!    enddo

!    L_theta = real(2*nperiod-1)
!    kpar = kpar / L_theta * kfilter
    
    zscale=1.0/sqrt(real(2*ntgrid))
    
! FFT work array: 
    nzwork=ccfft_work0 + ccfft_work1*2*ntgrid
    
! trig table: 
    nztable = ccfft_table0 + ccfft_table1*2*ntgrid
    
    if (nzwork > 0 .and. nztable > 0) then
       allocate(zwork(nzwork), ztable(nztable)) 

       b = 0.
       call ccfft(0, 2*ntgrid, zscale, b, b, ztable, zwork, 0)
    else
       ! no fft will be done
    end if

  end subroutine init_zf
    
  subroutine kz_spectrum (an, an2, ntgrid, ntheta0, naky)

    complex, dimension (:,:,:) :: an, an2
    integer, intent (in) :: ntgrid, ntheta0, naky
    integer :: it, ik

    do ik = 1, naky
       do it = 1, ntheta0
          call ccfft (-1, 2*ntgrid, zscale, an(:,it,ik), an2(:,it,ik), ztable, zwork, 0)
       end do
    end do

    an2 = conjg(an2)*an2

  end subroutine kz_spectrum

!  subroutine par_filter(an)

!    use fft_work
!    use theta_grid, only: ntgrid
!    use kt_grids, only: naky, ntheta0
!    complex, dimension(:,:,:) :: an
!    complex, dimension(2*ntgrid) :: bn
!    integer :: it, ik

!    if (allocated(tablekp)) then
!       do ik = 1, naky
!          do it = 1, ntheta0
!             call ccfft(-1, 2*ntgrid, scale, an(:,it,ik), bn, tablekp, work, 0)
!             bn=bn/(1+kpar**4)
!             call ccfft( 1, 2*ntgrid, scale, bn, an(:,it,ik), tablekp, work, 0)
!          enddo
!       enddo
!       
!! take care of last grid point in theta; should not matter
!       an(2*ntgrid+1,:,:) = an(1,:,:)
!    else
!       ! no fft in this case
!    end if

!  end subroutine par_filter

end module gs2_transforms

