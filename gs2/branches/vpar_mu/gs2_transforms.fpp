! Modifications for using FFTW version 3:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2009 
!                                 on behalf of the HECToR project

# include "define.inc"

module gs2_transforms

  use redistribute, only: redist_type
  use fft_work, only: fft_type

  implicit none

  public :: init_transforms
  public :: init_zf, kz_spectrum
  public :: transform_x, transform_y, transform2
  public :: inverse_x, inverse_y, inverse2

  private

  interface transform_x
     module procedure transform_x5d
!     module procedure transform_x3d
!     module procedure transform_x1d
  end interface

  interface transform_y
     module procedure transform_y5d
!     module procedure transform_y3d
!     module procedure transform_y1d
  end interface

  interface transform2
     module procedure transform2_5d
     module procedure transform2_4d
     module procedure transform2_3d
     module procedure transform2_2d
  end interface

  interface inverse_x
     module procedure inverse_x5d
!     module procedure inverse_x3d
!     module procedure inverse_x1d
  end interface

  interface inverse_y
     module procedure inverse_y5d
!     module procedure inverse_y3d
!     module procedure inverse_y1d
  end interface

  interface inverse2
     module procedure inverse2_5d
     module procedure inverse2_3d
     module procedure inverse2_2d
  end interface

  ! redistribution

  type (redist_type), save :: g2x, x2y

  ! fft

  type (fft_type) :: xf_fft, xb_fft, yf_fft, yb_fft, zf_fft
  type (fft_type) :: xf3d_cr, xf3d_rc

  logical :: xfft_initted = .false.

  complex, save, dimension (:, :), allocatable :: fft, xxf

contains

  subroutine init_transforms &
       (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)

    use mp, only: nproc
    use gs2_layouts, only: init_gs2_layouts
    use gs2_layouts, only: init_y_transform_layouts
    use gs2_layouts, only: init_x_transform_layouts

    implicit none

    integer, intent (in) :: ntgrid, naky, ntheta0, nvgrid, nmu, nspec
    integer, intent (in) :: nx, ny
    
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_gs2_layouts

    call init_y_redist (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)

! need these for movies
    call init_y_transform_layouts (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)
    call init_x_transform_layouts (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx)

    call init_y_fft

  end subroutine init_transforms

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine init_y_fft

    use gs2_layouts, only: xxf_lo, yxf_lo, g_lo
    use fft_work, only: init_crfftw, init_rcfftw, init_ccfftw

    implicit none

    logical :: initialized = .false.
    integer :: nb_ffts

    if (initialized) return
    initialized = .true.

    ! non-accelerated
    allocate (fft(yxf_lo%ny/2+1, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
    if (.not.allocated(xxf))then
       allocate (xxf(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc))
    endif
    
    if (.not. xfft_initted) then
       
       !JH FFTW plan creation for transform x5d and inverse
# if FFT == _FFTW_
       call init_ccfftw (xf_fft,  1, xxf_lo%nx)
       call init_ccfftw (xb_fft, -1, xxf_lo%nx)
# elif FFT == _FFTW3_
       ! number of ffts to be calculated
       !JH 7th December 2011
       !JH xxf_lo%ulim_alloc is used here rather than xxf_lo%lulim_proc
       !JH because there are situations where xxf_lo%llim_proc is greater 
       !JH than xxf_lo%ulim_proc and that would create a negative number 
       !JH of FFTs to be calculated.  However, xxf_lo%ulim_alloc is set
       !JH to be xxf_lo%llim_proc in this situation, and that will give 
       !JH 1 FFT to be calculated which the code can correctly undertake.
       nb_ffts = xxf_lo%ulim_alloc - xxf_lo%llim_proc + 1
       
       call init_ccfftw (xf_fft,  1, xxf_lo%nx, nb_ffts, xxf)
       call init_ccfftw (xb_fft, -1,  xxf_lo%nx, nb_ffts, xxf)
       
# endif
       xfft_initted = .true.
    end if
    
    !JH FFTW plan creation for transform y5d and inverse
# if FFT == _FFTW_
    call init_crfftw (yf_fft,  1, yxf_lo%ny)
    call init_rcfftw (yb_fft, -1, yxf_lo%ny)
#elif FFT == _FFTW3_
    ! number of ffts to be calculated
    !JH 7th December 2011
    !JH yxf_lo%ulim_alloc is used here rather than yxf_lo%lulim_proc
    !JH because there are situations where yxf_lo%llim_proc is greater 
    !JH than yxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, yxf_lo%ulim_alloc is set
    !JH to be yxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    nb_ffts = yxf_lo%ulim_alloc - yxf_lo%llim_proc + 1
    
    call init_crfftw (yf_fft,  1, yxf_lo%ny, nb_ffts)
    call init_rcfftw (yb_fft, -1,  yxf_lo%ny, nb_ffts)
    
#endif

  end subroutine init_y_fft

  subroutine init_x_redist (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx)

    use gs2_layouts, only: init_x_transform_layouts
    use gs2_layouts, only: g_lo, xxf_lo, gidx2xxfidx, proc_id, idx_local
    use gs2_layouts, only: layout
    use mp, only: nproc
    use redistribute, only: index_list_type, init_redist, delete_list
    use redistribute, only: set_redist_character_type

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (4) :: from_low, from_high
    integer, dimension (2) :: to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx
    logical :: initialized = .false.

    integer :: iglo, iv, ig, it0, it, ixxf
    integer :: n, ip

    if (initialized) return
    initialized = .true.

    call init_x_transform_layouts &
         (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx)

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       do it0 = 1, ntheta0
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                call gidx2xxfidx (ig, iv, it0, iglo, g_lo, xxf_lo, it, ixxf)
                if (idx_local(g_lo,iglo)) &
                     nn_from(proc_id(xxf_lo,ixxf)) = nn_from(proc_id(xxf_lo,ixxf)) + 1
                if (idx_local(xxf_lo,ixxf)) &
                     nn_to(proc_id(g_lo,iglo)) = nn_to(proc_id(g_lo,iglo)) + 1
             end do
          end do
       end do
    end do

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
          allocate (from_list(ip)%fourth(nn_from(ip)))
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
       do it0 = 1, ntheta0
          do iv = -nvgrid, nvgrid
             do ig = -ntgrid, ntgrid
                call gidx2xxfidx (ig, iv, it0, iglo, g_lo, xxf_lo, it, ixxf)
                if (idx_local(g_lo,iglo)) then
                   ip = proc_id(xxf_lo,ixxf)
                   n = nn_from(ip) + 1
                   nn_from(ip) = n
                   from_list(ip)%first(n) = ig
                   from_list(ip)%second(n) = iv
                   from_list(ip)%third(n) = it0
                   from_list(ip)%fourth(n) = iglo
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
    end do

    from_low (1) = -ntgrid
    from_low (2) = -nvgrid
    from_low (3) = 1
    from_low (4) = g_lo%llim_proc

    to_low = xxf_lo%llim_proc
    
    to_high(1) = xxf_lo%nx
    to_high(2) = xxf_lo%ulim_alloc

    from_high(1) = ntgrid
    from_high(2) = nvgrid
    from_high(3) = ntheta0
    from_high(4) = g_lo%ulim_alloc

    call set_redist_character_type(g2x, 'g2x')

    call init_redist (g2x, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_x_redist

  subroutine init_y_redist (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)

    use gs2_layouts, only: init_y_transform_layouts
    use gs2_layouts, only: xxf_lo, yxf_lo, xxfidx2yxfidx, proc_id, idx_local
    use mp, only: nproc
    use redistribute, only: index_list_type, init_redist, delete_list
    use redistribute, only: set_redist_character_type

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (2) :: from_low, from_high, to_high
    integer :: to_low
    integer, intent (in) :: ntgrid, naky, ntheta0, nvgrid, nmu, nspec
    integer, intent (in) :: nx, ny

    integer :: it, ixxf, ik, iyxf
    integer :: n, ip
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_x_redist (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx)
    call init_y_transform_layouts &
         (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, nx, ny)

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
!
!CMR: loop over all xxf indices, find corresponding yxf indices
!
    do ixxf = xxf_lo%llim_world, xxf_lo%ulim_world
       do it = 1, yxf_lo%nx
!
!CMR obtain corresponding yxf indices
!
          call xxfidx2yxfidx (it, ixxf, xxf_lo, yxf_lo, ik, iyxf)
          if (idx_local(xxf_lo,ixxf)) then
!CMR: if xxf index local, set:
!         ip = corresponding yxf processor
!        from_list%first,second arrays = it,ixxf  (ie xxf indices)
!     later will send from_list to proc ip
             ip = proc_id(yxf_lo,iyxf)
             n = nn_from(ip) + 1
             nn_from(ip) = n
             from_list(ip)%first(n) = it
             from_list(ip)%second(n) = ixxf
          end if
          if (idx_local(yxf_lo,iyxf)) then
!CMR: if yxf index local, set ip to corresponding xxf processor
!     set to_list%first,second arrays = ik,iyxf  (ie yxf indices)
!     will receive to_list from ip
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

    call set_redist_character_type(x2y, 'x2y')

    call init_redist (x2y, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)
    
  end subroutine init_y_redist

  subroutine transform_x5d (g, xxf)

    use gs2_layouts, only: xxf_lo, g_lo, it_idx, ik_idx, imu_idx
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather
    use mp, only: mp_abort
    use theta_grid, only: theta
    use vpamu_grids, only: anon

    implicit none

    complex, dimension (-xxf_lo%ntgrid:,-xxf_lo%nvgrid:,:,g_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    complex, dimension(:), allocatable :: aux
    integer :: i!, ik, it, imu, iglo, iv

!    complex, dimension (-xxf_lo%ntgrid:xxf_lo%ntgrid, -xxf_lo%nvgrid:xxf_lo%nvgrid, g_lo%llim_proc:g_lo%ulim_alloc) :: g0

    call prof_entering ("transform_x5d", "gs2_transforms")

    ! intent statement in gather actually makes this next line non-standard: 
    xxf = 0.
!CMR, 7/3/2011: gather pulls appropriate pieces of g onto this processor for
!    local Fourier transform in x, and may also pad with zeros for dealiasing
    call gather (g2x, g, xxf)

!    call write_mpdistxxf (xxf, '.xxf')
!    call mp_abort ('transform_x5d')

    ! g0 = 0.0
    ! do iglo = g_lo%llim_proc, g_lo%ulim_proc
    !    ik = ik_idx(g_lo,iglo)
    !    it = it_idx(g_lo,iglo)
    !    imu = imu_idx(g_lo,iglo)
    !    if (ik /= 2 .or. it /= 3) cycle
    !    do iv = -g_lo%nvgrid, g_lo%nvgrid
    !       g0(:,iv,iglo) = exp(-theta**2)*anon(:,iv,imu)

    !    end do
    ! end do

!    call gather (g2x, g0, xxf)

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH xxf_lo%ulim_alloc is used here rather than xxf_lo%lulim_proc
    !JH because there are situations where xxf_lo%llim_proc is greater 
    !JH than xxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, xxf_lo%ulim_alloc is set
    !JH to be xxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = xxf_lo%ulim_alloc - xxf_lo%llim_proc + 1

    allocate (aux(xf_fft%n))
    call fftw_f77 (xf_fft%plan, i, xxf, 1, xxf_lo%nx, aux, 0, 0)
    deallocate (aux)
# elif FFT == _FFTW3_
    call dfftw_execute(xf_fft%plan)
# endif

!    call write_mpdistxxf (xxf, '.xxf')
!    call mp_abort ('transform_x5d')

    call prof_leaving ("transform_x5d", "gs2_transforms")
  end subroutine transform_x5d

  subroutine inverse_x5d (xxf, g)

    use gs2_layouts, only: xxf_lo, g_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter

    implicit none

    complex, dimension (:,xxf_lo%llim_proc:), intent (in out) :: xxf
    complex, dimension (-xxf_lo%ntgrid:,-xxf_lo%nvgrid:,:,g_lo%llim_proc:), intent (out) :: g
    complex, dimension(:), allocatable :: aux
    integer :: i

    call prof_entering ("inverse_x5d", "gs2_transforms")

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH xxf_lo%ulim_alloc is used here rather than xxf_lo%lulim_proc
    !JH because there are situations where xxf_lo%llim_proc is greater 
    !JH than xxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, xxf_lo%ulim_alloc is set
    !JH to be xxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = xxf_lo%ulim_alloc - xxf_lo%llim_proc + 1

    allocate (aux(xb_fft%n))
    call fftw_f77 (xb_fft%plan, i, xxf, 1, xxf_lo%nx, aux, 0, 0)
    deallocate (aux)
# elif FFT == _FFTW3_
    call dfftw_execute(xb_fft%plan)
# endif

    call scatter (g2x, xxf, g)

    call prof_leaving ("inverse_x5d", "gs2_transforms")
  end subroutine inverse_x5d

  subroutine transform_y5d (xxf, yxf)

    use gs2_layouts, only: xxf_lo, yxf_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: gather

    implicit none

    complex, dimension (:,xxf_lo%llim_proc:), intent (in) :: xxf
# ifdef FFT
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf
# else
    real, dimension (:,yxf_lo%llim_proc:) :: yxf
# endif
    integer :: i

    call prof_entering ("transform_y5d", "gs2_transforms")

    fft = 0.
    call gather (x2y, xxf, fft)

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH yxf_lo%ulim_alloc is used here rather than yxf_lo%lulim_proc
    !JH because there are situations where yxf_lo%llim_proc is greater 
    !JH than yxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, yxf_lo%ulim_alloc is set
    !JH to be yxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = yxf_lo%ulim_alloc - yxf_lo%llim_proc + 1
    call rfftwnd_f77_complex_to_real (yf_fft%plan, i, fft, 1, yxf_lo%ny/2+1, yxf, 1, yxf_lo%ny)
# elif FFT == _FFTW3_

    call dfftw_execute_dft_c2r (yf_fft%plan, fft, yxf)
# endif

    call prof_leaving ("transform_y5d", "gs2_transforms")
  end subroutine transform_y5d

  subroutine inverse_y5d (yxf, xxf)
    use gs2_layouts, only: xxf_lo, yxf_lo
    use prof, only: prof_entering, prof_leaving
    use redistribute, only: scatter
    implicit none
    real, dimension (:,yxf_lo%llim_proc:), intent (in out) :: yxf
    complex, dimension (:,xxf_lo%llim_proc:), intent (out) :: xxf
    integer :: i 

    call prof_entering ("inverse_y5d", "gs2_transforms")

    ! do ffts
# if FFT == _FFTW_
    !JH 7th December 2011
    !JH yxf_lo%ulim_alloc is used here rather than yxf_lo%lulim_proc
    !JH because there are situations where yxf_lo%llim_proc is greater 
    !JH than yxf_lo%ulim_proc and that would create a negative number 
    !JH of FFTs to be calculated.  However, yxf_lo%ulim_alloc is set
    !JH to be yxf_lo%llim_proc in this situation, and that will give 
    !JH 1 FFT to be calculated which the code can correctly undertake.
    i = yxf_lo%ulim_alloc - yxf_lo%llim_proc + 1
    call rfftwnd_f77_real_to_complex (yb_fft%plan, i, yxf, 1, yxf_lo%ny, fft, 1, yxf_lo%ny/2+1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_r2c (yb_fft%plan, yxf, fft)
# endif

    call scatter (x2y, fft, xxf)

    call prof_leaving ("inverse_y5d", "gs2_transforms")
  end subroutine inverse_y5d

  subroutine transform2_5d (g, yxf)

    use gs2_layouts, only: g_lo, yxf_lo, ik_idx, it_idx, imu_idx
    use theta_grid, only: theta
    use vpamu_grids, only: anon
    use mp, only: mp_abort

    implicit none

    complex, dimension (-g_lo%ntgrid:,-g_lo%nvgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    real, dimension (:,yxf_lo%llim_proc:), intent (out) :: yxf
    integer :: iglo!, it, ik, imu, iv

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       if (ik_idx(g_lo, iglo) == 1) cycle
       g(:,:,:,iglo) = g(:,:,:,iglo) / 2.0
    end do

    ! g = 0.0
    ! do iglo = g_lo%llim_proc, g_lo%ulim_proc
    !    ik = ik_idx(g_lo,iglo)
    !    it = it_idx(g_lo,iglo)
    !    imu = imu_idx(g_lo,iglo)
    !    if (ik /= 2 .or. it /= 3) cycle
    !    do iv = -g_lo%nvgrid, g_lo%nvgrid
    !       g(:,iv,iglo) = exp(-theta**2)*anon(:,iv,imu)
    !    end do
    ! end do

    call transform_x (g, xxf)
!    call write_mpdistxxf (xxf, '.xxf')
    call transform_y (xxf, yxf)

!    call write_mpdistyxf (yxf, '.yxf')

!    call mp_abort ('testing xxf layout')

  end subroutine transform2_5d

  subroutine inverse2_5d (yxf, g)

    use gs2_layouts, only: g_lo, yxf_lo, ik_idx

    implicit none

    real, dimension (:,yxf_lo%llim_proc:), intent (in out) :: yxf
    complex, dimension (-yxf_lo%ntgrid:,-yxf_lo%nvgrid:,:,g_lo%llim_proc:), intent (out) :: g
    integer :: iglo

    call inverse_y (yxf, xxf)
    call inverse_x (xxf, g)

    g = g * xb_fft%scale * yb_fft%scale

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       if (ik_idx(g_lo, iglo) == 1) cycle
       g(:,:,:,iglo) = g(:,:,:,iglo) * 2.0
    end do

  end subroutine inverse2_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_3d (nny_in, nnx_in, how_many_in)

    use fft_work, only: init_crfftw, init_rcfftw, delete_fft

    logical :: initialized = .false.
    integer :: nny_in, nnx_in, how_many_in
    integer, save :: nnx, nny, how_many

    if (initialized) then
       if (nnx /= nnx_in .or. nny /= nny_in) then
          call delete_fft(xf3d_cr)
          call delete_fft(xf3d_rc)
# if FFT == _FFTW3_
       elseif ( how_many /= how_many_in) then
          call delete_fft(xf3d_cr)
          call delete_fft(xf3d_rc)
# endif
       else
          return
       end if
    end if
    initialized = .true.
    nny = nny_in
    nnx = nnx_in
    how_many = how_many_in

#if FFT == _FFTW_    

    call init_crfftw (xf3d_cr,  1, nny, nnx)
    call init_rcfftw (xf3d_rc, -1, nny, nnx)

#elif FFT == _FFTW3_

    call init_crfftw (xf3d_cr,  1, nny, nnx, how_many)
    call init_rcfftw (xf3d_rc, -1, nny, nnx, how_many)

#endif

  end subroutine init_3d

  subroutine transform2_3d (phi, phixf, nny, nnx)

    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, nx, aky

    implicit none

    integer :: nnx, nny
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,-ntgrid:), intent (out) :: phixf  
    real, dimension (:,:,:), allocatable :: phix
    complex, dimension (:,:,:), allocatable :: aphi
    real :: fac
    integer :: ig, ik, it, i

! scale, dealias and transpose

    call init_3d (nny, nnx, 2*ntgrid+1)

    allocate (phix (-ntgrid:ntgrid, nny, nnx))
    allocate (aphi (-ntgrid:ntgrid, nny/2+1, nnx))
    aphi = 0.

    do ik=1,naky
       fac = 0.5
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it=1,(ntheta0+1)/2
          do ig=-ntgrid, ntgrid
             aphi(ig,ik,it) = phi(ig,it,ik)*fac
          end do
       end do
       do it=(ntheta0+1)/2+1,ntheta0
          do ig=-ntgrid, ntgrid
!CMR, 30/3/2010: bug fix to replace nx by nnx on next line
             aphi(ig,ik,it-ntheta0+nnx) = phi(ig,it,ik)*fac
          end do
       end do
    end do

! transform
    i = 2*ntgrid+1
# if FFT == _FFTW_
    call rfftwnd_f77_complex_to_real (xf3d_cr%plan, i, aphi, i, 1, phix, i, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_c2r (xf3d_cr%plan, aphi, phix)
# endif

    do it=1,nnx
       do ik=1,nny
          do ig=-ntgrid, ntgrid
             phixf (it,ik,ig) = phix (ig,ik,it)
          end do
       end do
    end do

    deallocate (aphi, phix)

  end subroutine transform2_3d
  
  subroutine inverse2_3d (phixf, phi, nny, nnx)
!CMR, 30/4/2010:  
! Fixed up previously buggy handling of dealiasing.
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky

    implicit none

    real, dimension (:,:,-ntgrid:):: phixf
    complex, dimension (-ntgrid:,:,:) :: phi
    integer :: nnx, nny
    complex, dimension (:,:,:), allocatable :: aphi
    real, dimension (:,:,:), allocatable :: phix
    real :: fac
    integer :: i, ik, it, ig

    allocate (aphi (-ntgrid:ntgrid, nny/2+1, nnx))
    allocate (phix (-ntgrid:ntgrid, nny, nnx))

    do it=1,nnx
       do ik=1,nny
          do ig=-ntgrid, ntgrid
             phix (ig,ik,it) = phixf (it,ik,ig)
          end do
       end do
    end do

! transform
    i = 2*ntgrid+1
# if FFT == _FFTW_
    call rfftwnd_f77_real_to_complex (xf3d_rc%plan, i, phix, i, 1, aphi, i, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_r2c (xf3d_rc%plan, phix, aphi)
# endif

! dealias and scale
    do it=1,(ntheta0+1)/2
       do ik=1,naky
          fac = 2.0
          if (aky(ik) < epsilon(0.0)) fac = 1.0
          do ig=-ntgrid, ntgrid
             phi (ig,it,ik) = aphi (ig,ik,it)*fac*xf3d_rc%scale
          end do
       end do
    end do

!CMR, 30/4/2010:  fixed up previously buggy handling of dealiasing
    do it=(ntheta0+1)/2+1,ntheta0
       do ik=1,naky
          fac = 2.0
          if (aky(ik) < epsilon(0.0)) fac = 1.0
          do ig=-ntgrid, ntgrid
             phi (ig,it,ik) = aphi (ig,ik,it-ntheta0+nnx)*fac*xf3d_rc%scale
          end do
       end do
    end do

    deallocate (aphi, phix)

  end subroutine inverse2_3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform2_2d (phi, phixf, nny, nnx)

    use fft_work, only: FFTW_BACKWARD, delete_fft, init_crfftw
    use kt_grids, only: naky, nakx => ntheta0, nx, aky, akx

    implicit none

    integer :: nnx, nny
    complex, intent (in) :: phi(:,:)
    real, intent (out) :: phixf  (:,:)
    real, allocatable :: phix(:,:)
    complex, allocatable :: aphi(:,:)
    real :: fac
    integer :: ik, it
    type (fft_type) :: xf2d

#if FFT == _FFTW_

    call init_crfftw (xf2d, FFTW_BACKWARD, nny, nnx)

#elif FFT == _FFTW3_

    call init_crfftw (xf2d, FFTW_BACKWARD, nny, nnx, 1)

#endif

    allocate (phix (nny, nnx))
    allocate (aphi (nny/2+1, nnx))
    phix(:,:)=0.; aphi(:,:)=cmplx(0.,0.)

! scale, dealias and transpose
    do ik=1,naky
       fac = 0.5
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it=1,(nakx+1)/2
          aphi(ik,it) = phi(it,ik)*fac
       end do
       do it=(nakx+1)/2+1,nakx
          aphi(ik,it-nakx+nx) = phi(it,ik)*fac
       end do
    end do

! transform
# if FFT == _FFTW_
    call rfftwnd_f77_complex_to_real (xf2d%plan, 1, aphi, 1, 1, phix, 1, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_c2r (xf2d%plan, aphi, phix)
# endif

    phixf(:,:)=transpose(phix(:,:))

    deallocate (aphi, phix)
!RN> this statement causes error for lahey with DEBUG. I do not know why
!    call delete_fft(xf2d)
  end subroutine transform2_2d


  subroutine inverse2_2d (phixf, phi, nny, nnx)

    use fft_work, only: FFTW_FORWARD, delete_fft, init_rcfftw
    use kt_grids, only: naky, nakx => ntheta0, aky

    implicit none

    real, intent(in) :: phixf(:,:)
    complex, intent(out) :: phi(:,:)
    integer :: nnx, nny
    complex, allocatable :: aphi(:,:)
    real, allocatable :: phix(:,:)
    real :: fac
    integer :: ik, it
    type (fft_type) :: xf2d

#if FFT == _FFTW_

    call init_rcfftw (xf2d, FFTW_FORWARD, nny, nnx)

#elif FFT == _FFTW3_

    call init_rcfftw (xf2d, FFTW_FORWARD, nny, nnx, 1)

#endif    

    allocate (aphi (nny/2+1, nnx))
    allocate (phix (nny, nnx))
    phix(:,:)=cmplx(0.,0.); aphi(:,:)=cmplx(0.,0.)

    phix(:,:)=transpose(phixf(:,:))

! transform
# if FFT == _FFTW_
    call rfftwnd_f77_real_to_complex (xf2d%plan, 1, phix, 1, 1, aphi, 1, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_r2c (xf2d%plan, phix, aphi)
# endif

! scale, dealias and transpose
    do it=1,nakx
       do ik=1,naky
          fac = 2.0
          if (aky(ik) < epsilon(0.0)) fac = 1.0
          phi(it,ik) = aphi(ik,it)*fac*xf2d%scale
       end do
    end do

    deallocate (aphi, phix)
!RN> this statement causes error for lahey with DEBUG. I do not know why
!    call delete_fft(xf2d)
  end subroutine inverse2_2d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform2_4d (den, phixf, nny, nnx)
!CMR, 30/3/2010: den input has 4th index being species, 
!     but transform only operates for species=1
!     anyone who uses this routine should be aware of/fix this!
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, nx, aky

    implicit none

    integer :: nnx, nny
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: den
    real, dimension (:,:,-ntgrid:), intent (out) :: phixf  
    real, dimension (:,:,:), allocatable :: phix
    complex, dimension (:,:,:), allocatable :: aphi
    real :: fac
    integer :: ig, ik, it, i

! scale, dealias and transpose

    call init_3d (nny, nnx, 2*ntgrid+1)

    allocate (phix (-ntgrid:ntgrid, nny, nnx))
    allocate (aphi (-ntgrid:ntgrid, nny/2+1, nnx))
    aphi = 0.

    do ik=1,naky
       fac = 0.5
       if (aky(ik) < epsilon(0.)) fac = 1.0
       do it=1,(ntheta0+1)/2
          do ig=-ntgrid, ntgrid
             aphi(ig,ik,it) = den(ig,it,ik,1)*fac
          end do
       end do
       do it=(ntheta0+1)/2+1,ntheta0
          do ig=-ntgrid, ntgrid
!CMR, 30/3/2010: bug fix to replace nx by nnx on next line
             aphi(ig,ik,it-ntheta0+nnx) = den(ig,it,ik,1)*fac
          end do
       end do
    end do

! transform
    i = 2*ntgrid+1
# if FFT == _FFTW_
    call rfftwnd_f77_complex_to_real (xf3d_cr%plan, i, aphi, i, 1, phix, i, 1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft_c2r (xf3d_cr%plan, aphi, phix)
# endif

    do it=1,nnx
       do ik=1,nny
          do ig=-ntgrid, ntgrid
             phixf (it,ik,ig) = phix (ig,ik,it)
          end do
       end do
    end do

    deallocate (aphi, phix)

  end subroutine transform2_4d
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_zf (ntgrid, howmany)

    use fft_work, only: init_z

    implicit none
    integer, intent (in) :: ntgrid, howmany
    logical :: done = .false.

    if (done) return
    done = .true.

    call init_z (zf_fft, 1, 2*ntgrid, howmany)
    
  end subroutine init_zf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  subroutine kz_spectrum (an, an2, ntheta0, naky)

    complex, dimension (:,:,:), intent(in)  :: an
    complex, dimension (:,:,:), intent(out) :: an2
    integer, intent (in) :: ntheta0, naky

# if FFT == _FFTW_    
    call fftw_f77 (zf_fft%plan, ntheta0*naky, an, 1, zf_fft%n+1, an2, 1, zf_fft%n+1)
# elif FFT == _FFTW3_
    call dfftw_execute_dft(zf_fft%plan, an, an2)
# endif
    an2 = conjg(an2)*an2

  end subroutine kz_spectrum

  subroutine write_mpdistxxf (dist, extension)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: xxf_lo, ig_idx, iv_idx, ik_idx, is_idx
    use gs2_layouts, only: imu_idx, idx_local, proc_id
    use gs2_time, only: code_time
    use theta_grid, only: ntgrid, bmag, theta, shat
    use vpamu_grids, only: vpa, nvgrid, mu
    use kt_grids, only: theta0, lx, ly
    use constants, only: pi

    implicit none
    
    complex, dimension (:,xxf_lo%llim_proc:), intent (in) :: dist
    character (*), intent (in) :: extension
    
    integer :: i, j, ik, is, imu, ig, iv
    integer, save :: unit
    logical, save :: done = .false.
    complex :: gtmp
    
    if (.not. done) then
       if (proc0) call open_output_file (unit, trim(extension))
       do j=xxf_lo%llim_world, xxf_lo%ulim_world
          ig = ig_idx(xxf_lo, j)
          iv = iv_idx(xxf_lo, j)
          ik = ik_idx(xxf_lo, j)
          is = is_idx(xxf_lo, j) ; if (is /= 1) cycle
          imu = imu_idx(xxf_lo, j)
          do i = 1, xxf_lo%nx
             if (idx_local (xxf_lo, ig, iv, ik, imu, is)) then
                if (proc0) then
                   gtmp = dist(i,j)
                else
                   call send (dist(i,j), 0)
                end if
             else if (proc0) then
                call receive (gtmp, proc_id(xxf_lo, j))
             end if
             if (proc0) then
                write (unit,'(a1,6e14.4,i4,2e14.4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
                     (lx/xxf_lo%nx)*(i-1), ik, real(gtmp), aimag(gtmp)
!                write (unit,'(a1,5e14.4,2i4,2e14.4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
!                     i, ik, real(gtmp), aimag(gtmp)
             end if
          end do
          if (proc0) then
             write (unit,*)
             write (unit,*)
          end if
       end do
       if (proc0) call close_output_file (unit)
    end if
    
  end subroutine write_mpdistxxf

!   subroutine write_mpdistyxf (dist, extension)

!     use mp, only: proc0, send, receive
!     use file_utils, only: open_output_file, close_output_file
!     use gs2_layouts, only: yxf_lo, ig_idx, iv_idx, it_idx, is_idx
!     use gs2_layouts, only: imu_idx, idx_local, proc_id
!     use gs2_time, only: code_time
!     use theta_grid, only: ntgrid, bmag, theta, shat
!     use vpamu_grids, only: vpa, nvgrid, mu
!     use kt_grids, only: theta0, lx, ly
!     use constants, only: pi

!     implicit none
    
!     real, dimension (:,yxf_lo%llim_proc:), intent (in) :: dist
!     character (*), intent (in) :: extension
    
!     integer :: i, j, it, is, imu, ig, iv
!     integer, save :: unit
!     logical, save :: done = .false.
!     real :: gtmp
    
!     if (.not. done) then
!        if (proc0) call open_output_file (unit, trim(extension))
!        do j=yxf_lo%llim_world, yxf_lo%ulim_world
!           ig = ig_idx(yxf_lo, j)
!           iv = iv_idx(yxf_lo, j)
!           it = it_idx(yxf_lo, j)
!           is = is_idx(yxf_lo, j) ; if (is /= 1) cycle
!           imu = imu_idx(yxf_lo, j)
!           do i = 1, yxf_lo%ny
!              if (idx_local (yxf_lo, i)) then
!                 if (proc0) then
!                    gtmp = dist(i,j)
!                 else
!                    call send (dist(i,j), 0)
!                 end if
!              else if (proc0) then
!                 call receive (gtmp, proc_id(yxf_lo, j))
!              end if
!              if (proc0) then
!                 write (unit,'(a1,8e14.4)') "", code_time, theta(ig), vpa(iv), mu(imu), bmag(ig), &
!                      (lx/yxf_lo%nx)*(it-1), (ly/yxf_lo%ny)*(i-1), gtmp
!              end if
!           end do
!           if (proc0) then
!              write (unit,*)
!              write (unit,*)
!           end if
!        end do
!        if (proc0) call close_output_file (unit)
!     end if
    
!   end subroutine write_mpdistyxf

!   subroutine write_mpdistfft (dist, extension)

!     use mp, only: proc0, send, receive
!     use file_utils, only: open_output_file, close_output_file
!     use gs2_layouts, only: yxf_lo, ig_idx, iv_idx, it_idx, is_idx
!     use gs2_layouts, only: imu_idx, idx_local, proc_id
!     use gs2_time, only: code_time
!     use theta_grid, only: ntgrid, bmag, theta, shat
!     use vpamu_grids, only: vpa, nvgrid, mu
!     use kt_grids, only: theta0, lx, ly
!     use constants, only: pi

!     implicit none
    
!     complex, dimension (:,yxf_lo%llim_proc:), intent (in) :: dist
!     character (*), intent (in) :: extension
    
!     integer :: i, j, it, is, imu, ig, iv
!     integer, save :: unit
!     logical, save :: done = .false.
!     complex :: gtmp
    
!     if (.not. done) then
!        if (proc0) call open_output_file (unit, trim(extension))
!        do j=yxf_lo%llim_world, yxf_lo%ulim_world
!           ig = ig_idx(yxf_lo, j)
!           iv = iv_idx(yxf_lo, j)
!           it = it_idx(yxf_lo, j)
!           is = is_idx(yxf_lo, j) ; if (is /= 1) cycle
!           imu = imu_idx(yxf_lo, j)
!           do i = 1, yxf_lo%ny/2+1
!              if (idx_local (yxf_lo, i)) then
!                 if (proc0) then
!                    gtmp = dist(i,j)
!                 else
!                    call send (dist(i,j), 0)
!                 end if
!              else if (proc0) then
!                 call receive (gtmp, proc_id(yxf_lo, j))
!              end if
!              if (proc0) then
!                 write (unit,'(a1,6e14.4,i4,2e14.4)') "", code_time, theta(ig), &
!                      vpa(iv), mu(imu), bmag(ig), &
!                      (lx/yxf_lo%nx)*(it-1), i, real(gtmp), aimag(gtmp)
!              end if
!           end do
!           if (proc0) then
!              write (unit,*)
!              write (unit,*)
!           end if
!        end do
!        if (proc0) call close_output_file (unit)
!     end if
    
!   end subroutine write_mpdistfft

end module gs2_transforms
