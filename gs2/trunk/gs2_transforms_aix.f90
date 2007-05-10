module gs2_transforms

  use redistribute, only: redist_type
  use fft_work, only: fft_type

  implicit none

  public :: init_transforms
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
     module procedure transform2_5d
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
     module procedure inverse2_5d
     module procedure inverse2_3d
  end interface

  ! redistribution

  type (redist_type), save :: g2x, x2y

  ! fft

  type (fft_type) :: xf_fft, xb_fft, yf_fft, yb_fft

  logical :: xfft_initted = .false.
  logical :: yfft_initted = .false.

  integer :: igmin_proc, igmax_proc
  integer, dimension (:), allocatable :: igproc

  complex, dimension (:, :), allocatable :: fft, xxf

contains

  subroutine init_transforms &
       (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accelerate)
    implicit none
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec
    integer, intent (in) :: nx, ny
    logical, intent (out) :: accelerate

    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_3d_layouts

    call init_y_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny)
    call init_y_fft 

    accelerate = .false.

  end subroutine init_transforms

  subroutine init_x_transform (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    use fft_work, only: init_dcft
    use gs2_layouts, only: xxf_lo
    integer, intent (in) :: ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx

    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_3d_layouts

    call init_x_redist (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx)

    call init_dcft (xf_fft, -1, xxf_lo%nx)
    call init_dcft (xb_fft,  1, xxf_lo%nx)
    xfft_initted = .true.

  end subroutine init_x_transform

  subroutine init_y_fft 

    use gs2_layouts, only: xxf_lo, yxf_lo
    use fft_work, only: init_dcrft, init_drcft, init_dcft
    implicit none

    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    allocate (fft(yxf_lo%ny/2+1, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
    allocate (xxf(xxf_lo%nx, xxf_lo%llim_proc:xxf_lo%ulim_alloc))

    if (.not. xfft_initted) then
       call init_dcft (xf_fft, -1, xxf_lo%nx)
       call init_dcft (xb_fft,  1, xxf_lo%nx)
    end if

    call init_dcrft (yf_fft, -1, yxf_lo%ny)
    call init_drcft (yb_fft,  1, yxf_lo%ny)

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
    real, dimension(:), allocatable :: aux2
    integer :: i, idp, ipto, ipfrom, iadp

    call prof_entering ("transform_x5d", "gs2_transforms")

    ! intent statement in gather actually makes this next line non-standard: 
    xxf = 0.
    call gather (g2x, g, xxf)

    ! do ffts
    allocate (aux2(xf_fft%naux2))
    do i = xxf_lo%llim_proc, xxf_lo%ulim_proc
       call dcft (0, xxf(:,i), 1, 0, xxf(:,i), 1, 0, xxf_lo%nx, 1, &
            xf_fft%is, xf_fft%scale, &
            xf_fft%aux1, xf_fft%naux1, &
            aux2, xf_fft%naux2)
    end do
    deallocate (aux2)

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
    real, dimension(:), allocatable :: aux2
    integer :: i, idp, ipto, ipfrom, iadp

    call prof_entering ("inverse_x5d", "gs2_transforms")

    ! do ffts
    allocate (aux2(xb_fft%naux2))
    do i = xxf_lo%llim_proc, xxf_lo%ulim_proc
       call dcft (0, xxf(:,i), 1, 0, xxf(:,i), 1, 0, xxf_lo%nx, 1, &
            xb_fft%is, xb_fft%scale, &
            xb_fft%aux1, xb_fft%naux1, &
            aux2, xb_fft%naux2)
    end do
    deallocate (aux2)

    call scatter (g2x, xxf, g)

    call prof_leaving ("inverse_x5d", "gs2_transforms")
  end subroutine inverse_x5d

  subroutine transform_y5d (xxf, yxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use mp, only: iproc, nproc, send, receive
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather
    implicit none
    complex, dimension (:,xxf_lo%llim_proc:), intent (in) :: xxf
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf
    integer :: i, idp, ipto, ipfrom, iadp
!    complex, dimension (yxf_lo%ny/2+1,yxf_lo%llim_proc:yxf_lo%ulim_alloc) :: fft
    real, dimension(:), allocatable :: aux2

    call prof_entering ("transform_y5d", "gs2_transforms")

    fft = 0.
    call gather (x2y, xxf, fft)

    ! normalize
    fft(2:yxf_lo%naky,:) = fft(2:yxf_lo%naky,:)/2.0

    ! do ffts
    allocate (aux2(yf_fft%naux2))
    do i = yxf_lo%llim_proc, yxf_lo%ulim_proc
       call dcrft(0, fft(:,i), 0, yxf(:,i), 0, yxf_lo%ny, 1, &
            yf_fft%is,   yf_fft%scale, &
            yf_fft%aux1, yf_fft%naux1, &
            aux2,        yf_fft%naux2)
    end do
    deallocate (aux2)
    call prof_leaving ("transform_y5d", "gs2_transforms")
  end subroutine transform_y5d

  subroutine inverse_y5d (yxf, xxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use mp, only: iproc, nproc, send, receive
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in out) :: yxf
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    integer :: i, idp, ipto, ipfrom, iadp
!    complex, dimension (yxf_lo%ny/2+1,yxf_lo%llim_proc:yxf_lo%ulim_alloc) :: fft
    real, dimension(:), allocatable :: aux2

    call prof_entering ("inverse_y5d", "gs2_transforms")

    ! do ffts
    allocate (aux2(yb_fft%naux2))
    do i = yxf_lo%llim_proc, yxf_lo%ulim_proc
       call drcft(0, yxf(:,i), 0, fft(:,i), 0, yxf_lo%ny, 1, &
            yb_fft%is,   yb_fft%scale, &
            yb_fft%aux1, yb_fft%naux1, &
            aux2,        yb_fft%naux2)
       fft(2:,i) = fft(2:,i)*2.0
    end do
    deallocate (aux2)

    call scatter (x2y, fft, xxf)

    call prof_leaving ("inverse_y5d", "gs2_transforms")
  end subroutine inverse_y5d

  subroutine transform2_5d (g, yxf)
    use gs2_layouts, only: g_lo, xxf_lo, yxf_lo
    implicit none
    complex, dimension (:,:,g_lo%llim_proc:), intent (in) :: g
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf

    call transform_x (g, xxf)
    call transform_y (xxf, yxf)
  end subroutine transform2_5d

  subroutine inverse2_5d (yxf, g)
    use gs2_layouts, only: g_lo, xxf_lo, yxf_lo
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in out) :: yxf
    complex, dimension (:,:,g_lo%llim_proc:), intent (out) :: g

    xxf = 0.; g = 0.
    call inverse_y (yxf, xxf)
    call inverse_x (xxf, g)
  end subroutine inverse2_5d

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
    real, dimension(:), allocatable :: aux2
    integer :: ig, ik, nlo

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    if (igmin_proc <= igmax_proc) then
       allocate (aux2(xf_fft%naux2))
       do ig = igmin_proc, igmax_proc
          do ik = 1, naky
             phixxf(:,ik,ig) = 0.0
             phixxf(1:(ntheta0+1)/2,ik,ig) = phi(ig,ik,1:(ntheta0+1)/2)
             phixxf(nlo:xxf_lo%nx,ik,ig) &
                  = phi(ig,ik,(ntheta0+1)/2+1:ntheta0)
             call dcft (0, phixxf(:,ik,ig), 1, 0, phixxf(:,ik,ig), 1, 0, &
                  xxf_lo%nx, 1, &
                  xf_fft%is, xf_fft%scale, &
                  xf_fft%aux1, xf_fft%naux1, &
                  aux2, xf_fft%naux2)
          end do
       end do
       deallocate (aux2)
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
    real, dimension(:), allocatable :: aux2

    allocate (aux2(xb_fft%naux2))

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    do ig = -ntgrid, ntgrid
       do ik = 1, naky
          call dcft (0, phixxf(:,ik,ig), 1, 0, phixxf(:,ik,ig), 1, 0, &
               xxf_lo%nx, 1, &
               xb_fft%is, xb_fft%scale, &
               xb_fft%aux1, xb_fft%naux1, &
               aux2, xb_fft%naux2)
          phi(ig,1:(ntheta0+1)/2,ik) = f(1:(ntheta0+1)/2)
          phi(ig,(ntheta0+1)/2+1:ntheta0,ik) = f(nlo:xxf_lo%nx)
       end do
    end do

    deallocate (aux2)
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
    real, dimension(:), allocatable :: aux2

    if (igmin_proc <= igmax_proc) then
       allocate (aux2(yf_fft%naux2))
       do ig = igmin_proc, igmax_proc
          do ix = 1, xxf_lo%nx
             f(1) = real(phixxf(ix,1,ig))
             f(2:naky) = phixxf(ix,2:naky,ig)/2.0
             f(naky+1:) = 0.0
             call dcrft(0, f, 0, phixf(:,ix,ig), 0, yxf_lo%ny, 1, &
                  yf_fft%is, yf_fft%scale, &
                  yf_fft%aux1, yf_fft%naux1, &
                  aux2, yf_fft%naux2)
          end do
       end do
       deallocate (aux2)
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
    real, dimension(:), allocatable :: aux2

    allocate (aux2(yb_fft%naux2))
    do ig = -ntgrid, ntgrid
       do ix = 1, xxf_lo%nx
       call drcft(0, phixf(:,ix,ig), 0, f, 0, yxf_lo%ny, 1, &
            yb_fft%is,   yb_fft%scale, &
            yb_fft%aux1, yb_fft%naux1, &
            aux2,        yb_fft%naux2)

          phixxf(ix,1,ig) = f(1)
          phixxf(ix,2:naky,ig) = f(2:naky)*2.0
       end do
    end do
    deallocate (aux2)
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
    real, dimension(:), allocatable :: aux2

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    xf = 0.0
    xf(1:(ntheta0+1)/2) = f(1:(ntheta0+1)/2)
    xf(nlo:xxf_lo%nx) = f((ntheta0+1)/2+1:ntheta0)

    allocate (aux2(xf_fft%naux2))
    call dcft (0, xf, 1, 0, xf, 1, 0, xxf_lo%nx, 1, xf_fft%is, xf_fft%scale, &
         xf_fft%aux1, xf_fft%naux1, aux2, xf_fft%naux2)
    deallocate (aux2)

  end subroutine transform_x1d

  subroutine inverse_x1d (xf, f)
    use kt_grids, only: ntheta0
    use gs2_layouts, only: xxf_lo
    implicit none
    complex, dimension (:), intent (in out) :: xf
    complex, dimension (:), intent (out) :: f
    integer :: nlo
    real, dimension(:), allocatable :: aux2

    allocate (aux2(xb_fft%naux2))

    nlo = (ntheta0+1)/2+1+xxf_lo%nadd
    call dcft (0, xf, 1, 0, xf, 1, 0, xxf_lo%nx, 1, &
         xb_fft%is, xb_fft%scale, xb_fft%aux1, xb_fft%naux1, aux2, xb_fft%naux2)
    f(1:(ntheta0+1)/2) = xf(1:(ntheta0+1)/2)
    f((ntheta0+1)/2+1:ntheta0) = xf(nlo:xxf_lo%nx)

    deallocate (aux2)

  end subroutine inverse_x1d

  subroutine transform_y1d (f, xf)
    use kt_grids, only: naky
    use gs2_layouts, only: yxf_lo
    implicit none
    complex, dimension (:), intent (in) :: f
    real, dimension (:), intent (out) :: xf
    complex, dimension (yxf_lo%ny+1) :: fx
    real, dimension(:), allocatable :: aux2

    allocate (aux2(yf_fft%naux2))

    fx(1) = real(f(1))
    fx(2:naky) = f(2:naky)/2.0
    fx(naky+1:) = 0.0
    call dcrft(0, fx, 0, xf, 0, yxf_lo%ny, 1, &
         yf_fft%is,   yf_fft%scale, &
         yf_fft%aux1, yf_fft%naux1, &
         aux2,        yf_fft%naux2)

    deallocate (aux2)

  end subroutine transform_y1d

  subroutine inverse_y1d (xf, f)
    use kt_grids, only: naky
    use gs2_layouts, only: yxf_lo
    implicit none
    real, dimension (:), intent (in) :: xf
    complex, dimension (:), intent (out) :: f
    complex, dimension (yxf_lo%ny+1) :: fx
    real, dimension(:), allocatable :: aux2

    allocate (aux2(yb_fft%naux2))
    call drcft(0, xf, 0, fx, 0, yxf_lo%ny, 1, &
         yb_fft%is,   yb_fft%scale, &
         yb_fft%aux1, yb_fft%naux1, &
         aux2,        yb_fft%naux2)
    f(1) = fx(1)
    f(2:naky) = fx(2:naky)*2.0
    deallocate (aux2)
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

end module gs2_transforms
